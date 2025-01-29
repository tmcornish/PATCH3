#####################################################################################################
# - Creates maps of various systematics using frame metadata.
#####################################################################################################

#TODO: Try to determine a better way of configuring for different machines.

import os
import sys
from configuration import PipelineConfig as PC
from decasu.multi_healpix_mapper import MultiHealpixMapper
from decasu.configuration import Configuration
from output_utils import colour_string

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='makeMapsFromMetadata')



#######################################################
###############    START OF SCRIPT    #################
#######################################################

#retrieve the path of the Decasu config file
configfile = cf.auxfiles.decasu_config
#create Configuration object for Decasu
CONF = Configuration.load_yaml(configfile)


#the rest of the script now depends on whether it is being run locally or on glamdring
if cf.platform == 'local':
	#number of cores to use
	import multiprocessing as mp
	ncores = min(mp.cpu_count()-1, cf.ncores)
	#cycle through the fields being analysed
	for fd in cf.fields:
		print(colour_string(fd.upper(), 'orange'))
		#directory containing outputs for this field
		OUT = cf.paths.out + fd + '/'
		PATH_SYST = OUT + 'systmaps/'

		#set up Decasu mapper
		mapper = MultiHealpixMapper(CONF, PATH_SYST, ncores=ncores)

		#check if metadata was split by filter
		if cf.split_by_band:
			for b in cf.bands:
				infile = f'{cf.path.data}{cf.data_files.metadata[:-5]}_{fd}_{b}.fits'
				#run the Decasu mapper for each band individually
				mapper(infile, bands=b, clear_intermediate_files=True)
		else:
			infile = f'{cf.path.data}{cf.data_files.metadata[:-5]}_{fd}.fits'
			bands = ','.join(cf.bands)
			#run the Decasu mapper for all bands
			mapper(infile, bands=b, clear_intermediate_files=True)

else:
	#retrieve the band for which the script is being run from the system arguments
	b = sys.argv[2]
	#environment variable specifying number of CPUs will differ on NERSC and glamdring
	if cf.platform == 'nersc':
		import multiprocessing as mp
		ncores = mp.cpu_count() - 1
	else:
		ncores = int(os.getenv('SLURM_NTASKS_PER_NODE'))
	#cycle through the fields being analysed
	for fd in cf.fields:
		print(colour_string(fd.upper(), 'orange'))
		#directory containing outputs for this field
		OUT = cf.paths.out + fd + '/'
		PATH_SYST = OUT + 'systmaps/'

		#set up Decasu mapper
		mapper = MultiHealpixMapper(CONF, PATH_SYST, ncores=ncores)

		infile = f'{cf.path.data}{cf.data_files.metadata[:-5]}_{fd}_{b}.fits'
		mapper(infile, bands=b, clear_intermediate_files=True)
		
	