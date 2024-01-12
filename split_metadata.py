#####################################################################################################
# - Creates maps of galaxy counts and galaxy density in each tomographic bin.
#####################################################################################################

import os
import config
import healpy as hp
import healsparse as hsp
import numpy as np
from astropy.table import Table


### SETTINGS ###
cf = config.splitMetadata


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#get the global fields being analysed
fields_global = cf.get_global_fields()

#load the metadata
t = Table.read(cf.metafile)

#cycle through the global fields
for g in fields_global:
	#directory in which the metadata for this field is to be stored
	PATH_OUT = f'{cf.PATH_OUT}{g}/'
	#filename for the metadata in this field
	fname = f'{PATH_OUT}{cf.metasplit}'
	#only bother with the next steps if the target file doesn't exist
	if os.path.exists(fname):
		print(f'Metadata file for {g} already exists; skipping...')
	else:
		#make the host directory if it doesn't exist
		if not os.path.exists(PATH_OUT):
			os.system(f'mkdir -p {PATH_OUT}')
		#define a mask for selecting data within the bounds of the field
		ra_min, ra_max, dec_min, dec_max = cf.bounds[g]
		if ra_max < ra_min:
			mask = ((t['ra2000'] >= ra_min) | (t['ra2000'] <= ra_max)) * (t['dec2000'] >= dec_min) * (t['dec2000'] <= dec_max)
		else:
			mask = (t['ra2000'] >= ra_min) * (t['ra2000'] <= ra_max) * (t['dec2000'] >= dec_min) * (t['dec2000'] <= dec_max)
		#save the masked table to the destination file
		t[mask].write(fname, format='fits')



