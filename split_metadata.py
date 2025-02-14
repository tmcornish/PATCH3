#####################################################################################################
# - Creates maps of galaxy counts and galaxy density in each tomographic bin.
#####################################################################################################

import os
import sys
from configuration import PipelineConfig as PC
from astropy.table import Table
from output_utils import colour_string


### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='splitMetadata')


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#load the metadata
t = Table.read(cf.data_files.metadata)

#output directory (same as the directory containing the downloaded data)
PATH_OUT = cf.paths.data
#cycle through the global fields
for fd in cf.fields:
	print(colour_string(fd, 'orange'))

	#define a mask for selecting data within the bounds of the field
	ra_min, ra_max, dec_min, dec_max = cf.get_field_boundaries(fd)
	if ra_max < ra_min:
		coord_mask = ((t['ra2000'] >= ra_min) | (t['ra2000'] <= ra_max)) * (t['dec2000'] >= dec_min) * (t['dec2000'] <= dec_max)
	else:
		coord_mask = (t['ra2000'] >= ra_min) * (t['ra2000'] <= ra_max) * (t['dec2000'] >= dec_min) * (t['dec2000'] <= dec_max)
	
	#determine whether to split the metadata by band (filter)
	if cf.split_by_band:
		for b in cf.bands.all:
			fname = f'{cf.data_files.metadata[:-5]}_{fd}_{b}.fits'

			#only bother with the next steps if the target file doesn't exist
			if os.path.exists(fname):
				print(f'Metadata file for {fd} in band {b} already exists; skipping...')
				continue
			
			print(b)
			#mask to select rows using the current filter
			band_mask = (t['filter'] == b)
			#see if other columns might also contain data for this band
			if b in cf.bands.altnames:
				for b_alt in cf.bands.altnames[b]:
					band_mask |= (t['filter'] == b_alt)
			#combine with the field mask
			mask = coord_mask * band_mask
		
			#save the masked table to the destination file
			t[mask].write(fname, format='fits')

	else:
		fname = f'{cf.data_files.metadata[:-5]}_{fd}.fits'
		#only bother with the next steps if the target file doesn't exist
		if os.path.exists(fname):
			print(f'Metadata file for {fd} already exists; skipping...')
		else:
			t[coord_mask].write(fname, format='fits')


