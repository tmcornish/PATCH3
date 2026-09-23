##############################################################################
# Classes for constructing maps from the catalogue data.
##############################################################################

from stage import baseStage
import healpy as hp
import healsparse as hsp
import numpy as np


class mapperBase(baseStage):
    '''
    Base class for constructing maps from catalogue data.
    '''
    def __init__(
            self,
            config_file
    ):
        super().__init__(config_file)
        self.catalogue = None
        self.nside = self.config.nside
        self.nside_cover = self.config.nside_coverage
        self.maps = {}
        self.vpix = None
        self.field = None
        self.combine_fields = self.config.combine_fields

    def get_data_and_pix_ids(self):
        '''
        Retrieves the position info and required data for mapping.

        Here this is just a placeholder; it will be defined for each subclass.
        '''
        pass

    def _compute_statistic(
        self,
        pix_data,
        data=None,
        operation='binary',
        pix_comp=None,
        Nmin=0
    ):
        '''
        Creates a map showing the chosen summary statistic for some data.

        Parameters
        ----------
        pix_data: numpy.ndarray[int]
            Array of pixel IDs corresponding to the data being summarised.

        data: numpy.ndarray (None)
            Array of values corresponding to each pixel ID. Can be None if
            `operation` is 'count' or 'binary'.

        operation: str ('binary')
            Summary statistic to be mapped. Can be one of the following:
            - 'binary' (pixels are 1 if occupied and 0 otherwise)
            - 'count' (pixels display the number of sources within)
            - 'mean' (pixels display the mean of the data within)
            - 'std' (pixels display the standard deviation of the data within)

        pix_comp: numpy.ndarray[int] (None)
            Pixels in which the statistic is to be computed. If None, will
            compute the statistic for all pixels containing at least `Nmin`
            sources.

        Nmin: int (0)
            Minimum number of sources required for the summary statistic to be
            computed. Only relevant if `pix_comp` is not None, in which case
            pixels with fewer than `Nmin` sources have their values derived by
            interpolating neighbouring pixels.

        Returns
        -------
        '''
        # Check operation is compatible if data is None
        if data is None and operation not in ['binary', 'count']:
            raise ValueError(
                '`operation` must be either "binary" or "count" if `data` is '
                'set to `None`.'
            )

        # Define dtype of map based on summary statistic
        if operation in ['binary', 'count']:
            dt = np.int8
        else:
            dt = np.float64

        # Construct empty map
        m = hsp.HealSparseMap.make_empty(
            self.nside_cover,
            self.nside,
            dt
        )

        # Use pix_data as pix_comp if None provided
        if pix_comp is None:
            pix_comp = pix_data

        # Count sources in each pixel
        pmax = np.max(pix_comp)
        N = np.bincount(pix_data, minlength=pmax)[pix_comp]
        occupied = N > 0
        Nocc = N[occupied]

        # If operation is 'count', fill map with counts
        if operation == 'count':
            m[pix_comp] = N

        # If 'binary', identify which pixels actually contain sources
        elif operation == 'binary':
            m[pix_comp[~occupied]] = 0
            m[pix_comp[occupied]] = 1

        # Otherwise, additional steps are required
        else:
            # Calculate sum of the data in each pixel
            dsum = np.bincount(
                pix_data,
                weights=data,
                minlength=pmax
            )[pix_comp]
            # Calculate mean in occupied pixels
            dmean = dsum[occupied] / Nocc

            # Populate map if operation is 'mean'
            if operation == 'mean':
                m[pix_comp[occupied]] = dmean

            # Otherwise, move on to computing the standard deviation
            else:
                # Compute sum of the data squared in each pixel
                dsqsum = np.bincount(
                    pix_data,
                    weights=data*data,
                    minlength=pmax
                )[pix_comp]
                # Compute mean of squared data
                dsqmean = dsqsum[occupied] / Nocc
                # Compute variance
                dvar = dsqmean + (Nocc - 2) * dmean * dmean
                # Set negative variance (from rounding errors) to 0
                dvar[dvar < 0.] = 0.
                # Square root to get standard deviation
                dstd = np.sqrt(dvar)
                m[pix_comp[occupied]] = dstd

        # Identify pixels with counts below the threshold
        pix_few = pix_comp[N < Nmin]
        if len(pix_few) > 0:
            # Get the coordinates of these pixels
            ra_few, dec_few = hp.pix2ang(
                self.nside,
                pix_few,
                nest=True,
                lonlat=True
            )
            # Use nearest-neighbour interpolation to set pixel values
            vals_interp = m.interpolate_pos(
                ra_few,
                dec_few,
                lonlat=True,
                allow_partial=True
            )
            m[pix_few] = vals_interp

        return m

    def build_maps(self):
        '''
        Constructs the maps required of this stage and stores them.

        Here this is just a placeholder; it will be defined for each subclass.
        '''
        pass

    def write_maps(self):
        '''
        Writes the maps to their target files.
        '''
        for name, map_main in self.maps.items():
            # If a dictionary of maps, save to single file using recarray
            if isinstance(map_main, dict):
                # Retrieve map names and dtypes
                dts = [
                    (label, m.dtype) for label, m in map_main.items()
                ]
                primary = dts[0][0]
                # Construct empty recarray
                map_out = hsp.HealSparseMap.make_empty(
                    self.nside_cover,
                    self.nside,
                    dtype=dts,
                    primary=primary
                )
                # Initialise pixels and fill them in
                map_out.update_values_pix(
                    self.vpix,
                    np.rec.fromarrays(
                        [m[self.vpix] for m in map_main.values()],
                        dtype=dts
                    )
                )
            # Otherwise, output is a single map
            else:
                map_out = map_main

            # Save to file
            path_out = self.config.paths.out + self.field
            outfile = f'{path_out}/{name}_nside{self.nside}.hsp'
            map_out.write(outfile, clobber=True)

    def combine_maps(self):
        '''
        Combines maps from each individual HSC field and writes to files.
        '''
        for name in self.maps:
            # Load the maps from each field
            maps = [
                hsp.HealSparseMap.read(
                    f'{self.config.paths.out}{fd}/'
                    f'{name}_nside{self.nside}'
                ) for fd in self.config.fields
            ]
            # Check if the maps are recarrays of multiple maps
            dt = maps[0].dtype
            ndt = len(dt)
            if ndt > 0:
                # Cycle through individual maps and combine from each field
                unions_indiv = [
                    hsp.operations.max_union(
                        [m[m.dtype.names[i]] for m in maps]
                    ) for i in range(ndt)
                ]
                # Initialise empty map recarray
                union = hsp.HealSparseMap.make_empty(
                    self.nside_cover,
                    self.nside,
                    dtype=dt,
                    primary=dt.names[0]
                )
                # Fill in values
                vpix_all = unions_indiv[0].valid_pixels
                union.update_values_pix(
                    vpix_all,
                    np.rec.fromarrays(
                        [m[self.vpix] for m in unions_indiv],
                        dtype=dt
                    )
                )
            else:
                union = hsp.operations.max_union(maps)

            # Write to file
            path_out = self.config.paths.out + 'combined/'
            outfile = f'{path_out}/{name}_nside{self.nside}.hsp'
            union.write(outfile, clobber=True)

    def run(self):
        '''
        Runs the successive stages to produce the output maps.
        '''
        for fd in self.config.fields:
            self.field = fd
            self.get_data_and_pix_ids()
            self.build_maps()
            self.write_maps()

            if self.combine_fields:
                self.combine_maps()
