##############################################################################
# Classes for ingesting the downloaded HSC data.
##############################################################################

import os
import glob
from stage import baseStage
from astropy.table import Table
import h5py
import numpy as np


############################
#          STAGES          #
############################

class ingestBase(baseStage):
    '''
    Base class for loading FITS files and converting to HDF5.
    '''
    def __init__(
            self,
            config_file
    ):
        super().__init__(config_file)
        self.infiles = {}
        self.outfile = None
        self.columns = self.config.cols_include
        self.run_name = self.config.run_name

    def define_io(self):
        '''
        Identifies the input files and sets the output file location.

        Since this is the base class, this is just a placeholder and should be
        defined in each subclass individually.
        '''
        raise TypeError(
            'Base class cannot be executed. Please execute one of the '
            'subclasses instead.'
        )

    def _select_columns(self, t):
        '''
        Identifies columns for inclusion in the output file.

        Parameters
        ----------
        t: astropy.table.Table
            Input table from which columns are to be selected.
        '''
        # Check if list of columns explicitly stated in config
        if self.columns is None:
            # Exclude '_isnull' columns
            self.columns = [
                col for col in t.colnames if not col.endswith('_isnull')
            ]
        else:
            # Check each column is actually in the table
            self.columns = [
                col for col in self.columns if col in t.colnames
            ]

    def _apply_cuts(self, sample, t):
        '''
        Applies any additional cuts not applied at the query stages.

        Cuts are specified via the config file in the form of True/False
        flags or min/max thresholds on certain columns. Note that this
        requires that the relevant columns have been downloaded.

        Parameters
        ----------
        t: astropy.table.Table
            Table of data to which the cuts will be applied.

        sample: str
            Name of the sample for which the cuts are relevant. Can be
            'all' if the cuts are agnostic to the sample.

        Returns
        -------
        t_masked: astropy.table.Table
            Trimmed version of the table with the cuts applied.
        '''
        cf = self.config
        mask = np.ones(len(t), dtype=bool)

        # Select sources where these columns are False
        for col in cf.cuts.false_cols[sample]:
            for b in cf.cuts.false_cols[sample][col]:
                if col == 'merge_peak':
                    mask *= ~t[f'{col}_{b}']
                else:
                    mask *= ~t[f'{b}_{col}']

        # Select sources where these columns are True
        for col in cf.cuts.true_cols[sample]:
            for b in cf.cuts.true_cols[sample][col]:
                if col == 'merge_peak':
                    mask *= t[f'{col}_{b}']
                else:
                    mask *= t[f'{b}_{col}']

        # Select sources above these minimum thresholds
        for col in cf.cuts.min_cols[sample]:
            for b in cf.cuts.min_cols[sample][col]:
                mask *= (t[f'{col}_{b}']
                         >= cf.cuts.min_cols[sample][col][b])

        # Select sources below these maximum thresholds
        for col in cf.cuts.max_cols[sample]:
            for b in cf.cuts.max_cols[sample][col]:
                mask *= (t[f'{col}_{b}']
                         < cf.cuts.max_cols[sample][col][b])

        return t[mask]

    def write_hdf5(self):
        '''
        Constructs the output HDF5 file and saves it.
        '''
        # Verify output path exists
        path, _ = os.path.split(self.outfile)
        if not os.path.exists(path):
            os.system(f'mkdir -p {path}')

        with h5py.File(self.outfile, 'w') as f:
            # Create a Group for each sample
            for sample in self.infiles:
                infiles = self.infiles[sample]
                gp = f.require_group(sample)
                # Load first input file and use it to define the HDF5 datasets
                t = Table.read(infiles[0])
                t.sort('object_id')

                # Apply additional cuts (if any)
                if self.config.cuts is not None:
                    self._apply_cuts(t)

                # For initialising resizable datasets
                rowcount = len(t)
                shape = (rowcount,)
                maxshape = (None,)

                # Identify columns to include
                self._select_columns(t)
                for col in self.columns:
                    data = t[col]
                    dt = data.dtype
                    # Create resizable dataset
                    dset = gp.create_dataset(
                        col,
                        shape=shape,
                        maxshape=maxshape,
                        dtype=dt
                    )
                    dset[:] = data

                # Cycle through remaining tables and append data
                if len(infiles) < 2:
                    continue
                for infile in infiles[1:]:
                    # Load data
                    t = Table.read(infile)
                    t.sort('object_id')
                    rows_now = len(t)

                    shape = (rowcount,)
                    for col in self.columns:
                        data = t[col]
                        dt = data.dtype
                        dset = gp.require_dataset(col, shape=shape, dtype=dt)

                        # Resize the dataset
                        dset.resize(rowcount + rows_now, axis=0)

                        # Write the next chunk
                        dset[rowcount:] = data

                    # Update row count
                    rowcount += rows_now

    def run(self):
        '''
        Runs the successive stages to produce the output file.
        '''
        for fd in self.config.fields:
            self.field = fd
            self.define_io()
            self.write_hdf5()


class ingestFlags(ingestBase):
    '''
    Stage for ingesting the dust attenuation data.
    '''
    def define_io(self):
        '''
        Identifies the input files and sets the output file location.
        '''
        # Identify required input files
        self.infiles['/'] = sorted(
            glob.glob(
                self.config.paths.data +
                f'flags/{self.config.run_name}/*_{self.field}_*.fits'
            )
        )

        # Set output file
        self.outfile = os.path.join(
            self.config.paths.out,
            self.field,
            'flags.hdf5'
        )


class ingestDustAttenuation(ingestBase):
    '''
    Stage for ingesting the dust attenuation data.
    '''
    def define_io(self):
        '''
        Identifies the input files and sets the output file location.
        '''
        # Identify required input files
        self.infiles['/'] = sorted(
            glob.glob(
                self.config.paths.data +
                f'dust_attenuation/*_{self.field}_*.fits'
            )
        )

        # Set output file
        self.outfile = os.path.join(
            self.config.paths.out,
            self.field,
            'dust_attenuations.hdf5'
        )


class ingestStarsForDepth(ingestBase):
    '''
    Stage for ingesting the stellar data required for making depth maps.
    '''
    def define_io(self):
        '''
        Identifies the input files and sets the output file location.
        '''
        # Identify required input files
        self.infiles['/'] = sorted(
            glob.glob(
                self.config.paths.data +
                f'stars/for_depth_map/*_{self.field}_*_fluxerrs.fits'
            )
        )

        # Set output file
        self.outfile = os.path.join(
            self.config.paths.out,
            self.field,
            'star_catalogue_for_depth_map.hdf5'
        )


class ingestGalaxies(ingestBase):
    '''
    Stage for ingesting the data for the galaxy samples being analysed.
    '''
    def define_io(self):
        '''
        Identifies the input files and sets the output file location.

        Downloaded data for each sample are compiled into a single HDF5 file
        with one Group per sample.
        '''
        # Identify files corresponding to each sample
        for sample in self.config.samples:
            self.infiles[sample] = sorted(
                glob.glob(
                    self.config.paths.data +
                    f'galaxies/{self.run_name}/*_{self.field}_*_{sample}.fits'
                )
            )

        # Set output file
        self.outfile = os.path.join(
            self.config.paths.out,
            self.field,
            'galaxy_catalogue.hdf5'
        )


class ingestLikeStars(ingestBase):
    '''
    Stage for ingesting the data for analogous stellar samples.
    '''
    def define_io(self):
        '''
        Identifies the input files and sets the output file location.
        '''
        # Identify files corresponding to each sample
        self.infiles['/'] = sorted(
            glob.glob(
                self.config.paths.data +
                f'stars/{self.run_name}/*_{self.field}_*.fits'
            )
        )

        # Set output file
        self.outfile = os.path.join(
            self.config.paths.out,
            self.field,
            'star_catalogue.hdf5'
        )


class ingestCosmosHSC(ingestBase):
    '''
    Stage for ingesting the data for analogous COSMOS galaxies.
    '''
    def define_io(self):
        '''
        Identifies the input files and sets the output file location.
        '''
        # Identify files corresponding to each sample
        self.infiles['/'] = sorted(
            glob.glob(
                self.config.paths.data +
                f'galaxies/{self.run_name}/*_cosmos.fits'
            )
        )

        # Set output file
        self.outfile = os.path.join(
            self.config.paths.out,
            self.field,
            'dir_catalogue.hdf5'
        )
