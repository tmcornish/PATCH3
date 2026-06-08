#####################################################################################################
# - Creates maps of galaxy counts and galaxy density in each tomographic bin.
#####################################################################################################

import os
import config
import healpy as hp
import healsparse as hsp
import numpy as np
import h5py
from astropy.table import Table
from output_utils import colour_string


### SETTINGS ###
cf = config.splitByPixel



#######################################################
###############    START OF SCRIPT    #################
#######################################################

#dictonary of the 'global' fields being analysed and the sub-fields belonging to them
f_in_g = cf.fields_in_global()

#cycle through each global field
for g in f_in_g.keys():
	print(colour_string(g))
	#load the metadata for this field
	meta = Table.read(f'{cf.PATH_OUT}{g}/{cf.metasplit}', format='fits')
	#define polygons for each frame described in the metadata
	ra_corners = np.array([meta['llcra'], meta['lrcra'], meta['urcra'], meta['ulrca']]).T
	dec_corners = np.array([meta['llcdec'], meta['lrcdec'], meta['urcdec'], meta['ulrcdec']]).T
	poly = [hsp.Polygon(ra=ra_corners[i], dec=dec_corners[i], value=1) for i in range(len(meta))]
	#identify the pixels in which each frame resides
	print('Identifying pixels for each metadata frame...')
	vpix_meta = [p.get_map(nside_coverage=cf.nside_cover, nside_sparse=cf.nside_hi, dtype=np.int8).degrade(cf.nside_cover).valid_pixels for p in poly]
	vpix_u = np.unique([v for vp in vpix_meta for v in vp])
	#cycle through each pixel
	print('Splitting metadata...')
	for v in vpix_u:
		#identifier for this pixel
		PID = f'{cf.nside_cover}_{v}'
		print(colour_string(PID, 'purple'))

		#output directory for (meta)data belonging to this pixel
		OUT = f'{cf.PATH_OUT}{g}/pix{PID}/'
		#if directory doesn't exist, make it
		if not os.path.exists(OUT):
			os.system(f'mkdir -p {OUT}')

		#create mask to identify frames overlapping the current pixel
		pixmask = [True if v in vp else False for vp in vpix_meta]
		#write those data to a new file
		meta[pixmask].write(f'{OUT}{cf.metasplit}', overwrite=True)


	#get the list of sub-fields belonging to this field
	fields_now = f_in_g[g]




'''
with h5py.File('test_pix8_167.hdf5', 'a') as hf:
	_ = hf.require_group('photometry')
	N = int(mask.sum())
	for i in range(4):
		for col in cat['photometry'].keys():
			col_data = cat[f'photometry/{col}'][mask]
			if col in hf['photometry'].keys():
				dset = hf[f'photometry/{col}']
				dset.resize((len(dset)+N,))
				dset[-N:] = col_data
			else:
				#continue
				dset = hf.create_dataset(f'photometry/{col}', shape=(N,), data=cat[f'photometry/{col}'][mask], maxshape=(None,), dtype=t[col].dtype)
			del col_data
'''