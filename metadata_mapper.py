##############################################################################
# Classes for creating various maps necessary for the analysis.
##############################################################################

import os
from stage import baseStage
from decasu.configuration import Configuration
from decasu.multi_healpix_mapper import MultiHealpixMapper


############################
#          STAGES          #
############################

class decasuMapperBase(baseStage):
    '''
    Base class for producing maps with decasu; should not be used dirctly.
    '''
    def __init__(
        self,
        config_file,
    ):
        super().__init__(config_file)
        # Default parameters for decasu
        self.decasu_params = {
            'outbase': 'decasu',
            'map_types': {'coverage': ['sum']},
            'nside': 32768,
            'use_two_amps': False,
            'arcsec_per_pix': 0.17,
            'maglim_aperture': 4.08,
            'maglim_nsig': 5.0,
            'zp_global': 27.0,
            'magzp_field': "zeropt",
            'exp_field': "visit",
            'ccd_field': "ccd",
            'band_field': "filter",
            'mjd_field': "mjd",
            'skyvar_field': "skylevel",
            'fwhm_field': "seeing",
            'bad_amps': {},
            'bad_ccds': [9],
            'extra_fields': {},
            'band_replacement': {'i2': 'i', 'r2': 'r'},
            'use_wcs': False,
            'ra_corner_fields': ['llcra', 'lrcra', 'urcra', 'ulcra'],
            'dec_corner_fields': ['llcdec', 'lrcdec', 'urcdec', 'ulcdec'],
            'latitude': 19.8230,
            'longitude': -155.4694,
            'elevation': 4205.0,
        }

    def update_config(self):
        '''
        Updates the decasu parameters using the config file.
        '''
        for p in self.config.decasu_params:
            self.decasu_params[p] = self.config.decasu_params[p]

    def build_config(self):
        '''
        Builds a Configuration object for running decasu.
        '''
        self.decasu_conf = Configuration(**self.decasu_params)

    def select_metadata(self):
        '''
        Identifies the metadata files to use for creating the requested maps.
        '''
        import glob

        path_meta = self.config.paths.data + 'metadata/'
        self.infiles = []
        self.bands = []
        self.fields = []

        # Cycle through fields
        for fd in self.config.fields:
            # Were metadata downloaded for each badn separately?
            if self.config.split_by_band:
                files_fd = sorted(glob.glob(path_meta + f'*_{fd}_?.fits'))
                # Only include those for requested bands
                files_fd[:] = [
                    f for f in files_fd
                    for b in self.config.bands_to_run
                    if f'_{b}' in f
                ]
                self.infiles.extend(files_fd)
                self.bands.extend(list(self.config.bands_to_run))
            else:
                files_fd = glob.glob(path_meta + f'*_{fd}.fits')
                self.infiles.append(files_fd[0])
                self.bands.append(','.join(list(self.config.bands_to_run)))

            self.fields.extend([fd] * len(files_fd))

        # If running with MPI, split the tasks and seelct for this rank
        ntasks = self.config.qsub.resources.select
        if ntasks > 1:
            from mpi4py import MPI
            rank = MPI.COMM_WORLD.Get_rank()

            self.infiles = self.split_list(self.infiles, ntasks)[rank]
            self.bands = self.split_list(self.bands, ntasks)[rank]

    def get_ncpus(self):
        '''
        Determines number of CPUs to use based on platform.
        '''
        if self.config.platform == 'local':
            import multiprocessing as mp
            # TODO: add config option to set ncpu manually?
            ncpus = mp.cpu_count() - 1
        else:
            ncpus = self.config.qsub.resources.ncpus
        return ncpus

    def run(self):
        '''
        Runs decasu to produce the requested maps.
        '''
        self.update_config()
        self.build_config()
        self.select_metadata()
        ncpus = self.get_ncpus()

        for infile, b, fd in zip(self.infiles, self.bands, self.fields):
            path_out = self.config.paths.out + fd + 'maps/'\
                + self.config.subdir

            if not os.path.exists(path_out):
                os.system(f'mkdir -p {path_out}')

            # Set up decasu mapper
            mapper = MultiHealpixMapper(
                self.decasu_conf,
                path_out,
                ncores=ncpus
            )

            # Run the mapper
            mapper(infile, bands=b, clear_intermediate_files=True)


class coverageMapper(decasuMapperBase):
    '''
    Stage for producing coverage maps to define survey geometry.
    '''
    def __init__(self, config_file):
        super().__init__(config_file)


class surveyPropertyMapper(decasuMapperBase):
    '''
    Stage for producing maps of various survey properties.
    '''
    def __init__(self, config_file):
        super().__init__(config_file)
