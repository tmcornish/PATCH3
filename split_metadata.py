#####################################################################################################
# - Creates maps of galaxy counts and galaxy density in each tomographic bin.
#####################################################################################################

import os
import config
import healpy as hp
import healsparse as hsp
import numpy as np
from astropy.table import Table
from output_utils import colour_string


### SETTINGS ###
cf = config.splitMetadata


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#get the global fields being analysed
fields_global = cf.get_global_fields()

#load the metadata
t = Table.read(cf.metafile)
#get the extension of the metadata file and its character length
ext = cf.metafile.split('.')[-1]
lenext = len(ext)

#cycle through the global fields
for g in fields_global:
	print(colour_string(g, 'orange'))
	#directory in which the metadata for this field is to be stored
	PATH_OUT = f'{cf.PATH_OUT}{g}/'
	#make the host directory if it doesn't exist
	if not os.path.exists(PATH_OUT):
		os.system(f'mkdir -p {PATH_OUT}')

	#define a mask for selecting data within the bounds of the field
	ra_min, ra_max, dec_min, dec_max = cf.bounds[g]
	if ra_max < ra_min:
		coord_mask = ((t['ra2000'] >= ra_min) | (t['ra2000'] <= ra_max)) * (t['dec2000'] >= dec_min) * (t['dec2000'] <= dec_max)
	else:
		coord_mask = (t['ra2000'] >= ra_min) * (t['ra2000'] <= ra_max) * (t['dec2000'] >= dec_min) * (t['dec2000'] <= dec_max)
	
	#determine whether to split the metadata by band (filter)
	if cf.split_by_band:
		for b in cf.bands:
			fname = f'{PATH_OUT}{cf.metasplit[:-(lenext+1)]}_{b}.fits'

			#only bother with the next steps if the target file doesn't exist
			if os.path.exists(fname):
				print(f'Metadata file for {g} in band {b} already exists; skipping...')
				continue
			
			print(b)
			#mask to select rows using the current filter
			band_mask = (t['filter'] == b)
			#see if other columns might also contain data for this band
			if b in cf.bands_alt:
				for b_alt in cf.bands_alt[b]:
					band_mask |= (t['filter'] == b_alt)
			#combine with the field mask
			mask = coord_mask * band_mask
		
			#save the masked table to the destination file
			t[mask].write(fname, format='fits')

	else:
		fname = f'{PATH_OUT}{cf.metasplit}'
		#only bother with the next steps if the target file doesn't exist
		if os.path.exists(fname):
			print(f'Metadata file for {g} already exists; skipping...')
		else:
			t[coord_mask].write(fname, format='fits')


