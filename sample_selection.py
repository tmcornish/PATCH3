#####################################################################################################
# - Identifies and flags galaxies in the cleaned catalogues as belonging to the desired
#   samples for clustering analysis.
#####################################################################################################

import sys
import h5py
import numpy as np
from configuration import PipelineConfig as PC
from output_utils import colour_string

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='sampleSelection')

#######################################################
###############    START OF SCRIPT    #################
#######################################################

#cycle through each of the fields
for fd in cf.fields:
	print(colour_string(fd.upper(), 'orange'))
	#output directory for this field
	OUT = cf.paths.out + fd
	#load the fully cleaned galaxy catalogue for this field
	with h5py.File(f'{OUT}/{cf.cats.main}', 'a') as cat:
		gp = cat['photometry']
		#define masks for each subsample as defined in the config file
		sample_masks = cf.get_samples(gp)
		#remove galaxies with secondary solutions at high redshift if told in config
		if cf.remove_pz_outliers:
			outlier_mask = ((gp['pz_err95_max_dnnz'][:] - gp['pz_err95_min_dnnz'][:] < 2.7)
						* (gp['pz_err95_max_mizu'][:] - gp['pz_err95_min_mizu'][:] < 2.7))
		else:
			outlier_mask = np.ones_like(sample_masks[0], dtype=bool)
		
		#open the summary file to write outputs as they are determined
		with open(f'{OUT}/{cf.sample_summary_file}', 'w') as outfile:
			outfile.write('Sample\tCounts\n')
			#cycle through the samples
			for sm, key in zip(sample_masks, cf.samples):
				sm *= outlier_mask
				#create a dataset for each of these masks
				d = gp.require_dataset(key, sm.shape, dtype=bool)
				d[:] = sm
				outfile.write(f'{key}\t{sm.sum()}\n')


	