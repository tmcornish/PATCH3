###############################################################
# - Creates maps of various systematics using frame metadata.
###############################################################

# TODO: Try to determine a better way of configuring for different machines.

import os
import sys
from configuration import PipelineConfig as PC
from decasu.multi_healpix_mapper import MultiHealpixMapper
from decasu.configuration import Configuration
from output_utils import colour_string

# SETTINGS #
config_file = sys.argv[1]
cf = PC(config_file, stage='makeMapsFromMetadata')

#######################################################
#                  START OF SCRIPT                    #
#######################################################

# Retrieve the path of the Decasu config file
configfile = cf.auxfiles.decasu_config
# Create Configuration object for Decasu
CONF = Configuration.load_yaml(configfile)

# Overwrite the output file basenames and NSIDE to match pipeline config
CONF.outbase = f'decasu_nside{cf.nside_hi}'
CONF.nside = cf.nside_hi
# Do the same with replacements for band names
for k in cf.bands.altnames:
    for j in cf.bands.altnames[k]:
        CONF.band_replacement[j] = k

# Rest of script now depends on whether it is being run locally or on glamdring
if cf.platform == 'local':
    # Number of cores to use
    import multiprocessing as mp
    ncores = min(mp.cpu_count()-1, cf.ncores)
    # Cycle through the fields being analysed
    for fd in cf.fields:
        print(colour_string(fd.upper(), 'orange'))
        # Directory containing outputs for this field
        OUT = cf.paths.out + fd + '/'
        PATH_SYST = OUT + 'systmaps/'

        # Set up Decasu mapper
        mapper = MultiHealpixMapper(CONF, PATH_SYST, ncores=ncores)

        # Check if metadata was split by filter
        if cf.split_by_band:
            for b in cf.bands:
                infile = f'{cf.path.data}'\
                         f'{cf.data_files.metadata[:-5]}_{fd}_{b}.fits'
                # Run the Decasu mapper for each band individually
                mapper(infile, bands=b, clear_intermediate_files=True)
        else:
            infile = f'{cf.path.data}{cf.data_files.metadata[:-5]}_{fd}.fits'
            bands = ','.join(cf.bands)
            # Run the Decasu mapper for all bands
            mapper(infile, bands=b, clear_intermediate_files=True)

else:
    # Retrieve band for which the script is being run from the system arguments
    b = sys.argv[2]
    # Environment variable specifying no. of CPUs differs on NERSC vs glamdring
    if cf.platform == 'nersc':
        import multiprocessing as mp
        ncores = mp.cpu_count() - 1
    else:
        ncores = int(os.getenv('SLURM_NTASKS_PER_NODE'))
    # Cycle through the fields being analysed
    for fd in cf.fields:
        print(colour_string(fd.upper(), 'orange'))
        # Directory containing outputs for this field
        OUT = cf.paths.out + fd + '/'
        PATH_SYST = OUT + 'systmaps/'

        # Set up Decasu mapper
        mapper = MultiHealpixMapper(CONF, PATH_SYST, ncores=ncores)

        infile = f'{cf.path.data}{cf.data_files.metadata[:-5]}_{fd}_{b}.fits'
        mapper(infile, bands=b, clear_intermediate_files=True)
