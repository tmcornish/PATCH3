##############################################################################
# Classes for ingesting the downloaded HSC data.
##############################################################################

import os
import glob
from stage import baseStage
from astropy.table import Table
import h5py


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

    def write_hdf5(self):
        '''
        Constructs the output HDF5 file and saves it.
        '''
        with h5py.File(self.outfile, 'w') as f:
            # Create a Group for each sample
            for sample in self.infiles:
                infiles = self.infiles[sample]
                gp = f.create_group(sample)
                # Load first input file and use it to define the HDF5 datasets
                t = Table.read(infiles[0])
                t.sort('object_id')

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
                    return
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
                    f'galaxies/*_{self.field}_*_{sample}.fits'
                )
            )

        # Set output file
        self.outfile = os.path.join(
            self.config.paths.out,
            self.field,
            'galaxy_catalogue.hdf5'
        )


class ingestStars(ingestBase):
    '''
    Stage for ingesting the data for analogous stellar samples.

    TODO: Join this to ingestGalaxies somehow, rather than having two separate
    but very similar stages?
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
                    f'stars/*_{self.field}_*_{sample}.fits'
                )
            )

        # Set output file
        self.outfile = os.path.join(
            self.config.paths.out,
            self.field,
            'star_catalogue.hdf5'
        )
