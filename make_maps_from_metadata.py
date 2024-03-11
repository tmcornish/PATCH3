#####################################################################################################
# - Creates maps of various systematics using frame metadata.
#####################################################################################################

#TODO: write function to create decasu config file from pipeline config info
#TODO: Integrate properly into pipeline (currently written to be run via SLURM bash script for multiprocessing)

import os, sys
import config
from decasu.multi_healpix_mapper import MultiHealpixMapper
from decasu.configuration import Configuration
from astropy.table import Table
from output_utils import colour_string
import glob

### SETTINGS ###
cf = config.makeMapsFromMetadata



#######################################################
###############    START OF SCRIPT    #################
#######################################################

#retrieve the path of the Decasu config file
configfile = cf.configfile
#create Configuration object for Decasu
CONF = Configuration.load_yaml(configfile)

#get the extension of the metadata file and its character length
ext = cf.metafile.split('.')[-1]
lenext = len(ext)

#the rest of the script now depends on whether it is being run locally or on glamdring
if cf.LOCAL:
	#number of cores to use
	import multiprocessing as mp
	ncores = min(mp.cpu_count()-1, cf.ncores)
	#cycle through the fields being analysed
	for fd in cf.get_global_fields():
		print(colour_string(fd.upper(), 'orange'))
		#directory containing outputs for this field
		OUT = cf.PATH_OUT + fd + '/'
		PATH_SYST = OUT + 'systmaps/'

		#set up Decasu mapper
		mapper = MultiHealpixMapper(CONF, PATH_SYST, ncores=ncores)

		#check if metadata was split by filter
		if cf.split_by_band:
			for b in cf.bands:
				infile = f'{OUT}{cf.metasplit[:-(lenext+1)]}_{b}.fits'
				#run the Decasu mapper for each band individually
				mapper(infile, bands=b, clear_intermediate_files=True)
		else:
			infile = f'{OUT}{cf.metasplit}'
			bands = ','.join(cf.bands)
			#run the Decasu mapper for all bands
			mapper(infile, bands=b, clear_intermediate_files=True)

else:
	#retrieve the band for which the script is being run from the system arguments
	b = sys.argv[1]
	#get the number of cores assigned from SLURM configuration
	ncores = int(os.getenv('SLURM_NTASKS_PER_NODE'))
	#cycle through the fields being analysed
	for fd in cf.get_global_fields():
		print(colour_string(fd.upper(), 'orange'))
		#directory containing outputs for this field
		OUT = cf.PATH_OUT + fd + '/'
		PATH_SYST = OUT + 'systmaps/'

		#set up Decasu mapper
		mapper = MultiHealpixMapper(CONF, PATH_SYST, ncores=ncores)

		infile = f'{OUT}{cf.metasplit[:-(lenext+1)]}_{b}.fits'
		mapper(infile, bands=b, clear_intermediate_files=True)
		
	