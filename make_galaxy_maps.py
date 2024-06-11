#####################################################################################################
# - Creates maps of galaxy counts and galaxy density in each tomographic bin.
#####################################################################################################

import os
import config
import healpy as hp
import healsparse as hsp
import numpy as np
from output_utils import colour_string
from map_utils import *
import h5py


### SETTINGS ###
cf = config.makeGalaxyMaps


###################
#### FUNCTIONS ####
###################

def makeNgalMaps(cat, mask, group=''):
	'''
	Creates galaxy count maps for each tomographic bin.

	Parameters
	----------
	cat: h5py.File object
		Catalogue containing (at least): RAs, Decs, and dust attenuation values at each of these
		coordinates in each of the specified bands.

	group: str
		Group within which the relevant data are expected to reside.

	mask: HealSparseMap or None
		If provided, will use the mask to determine the positions of valid pixels rather than
		the RAs and Decs from the catalogue.

	Returns
	-------
	ngal_maps: recarray
		recarray for which each entry is a galaxy counts map for a given redshift bins.
	'''

	print('Creating galaxy count maps...')
	#initialise a recarray to contain the maps for each band
	labels = [f'ngal_{i}' for i in range(len(cf.zbins)-1)]
	ngal_maps, _, _ = initialiseRecMap(cf.nside_lo, cf.nside_hi, labels, pixels=mask.vpix_nest, dtypes='f8')

	#cycle through the tomographic bins and count the galaxies in each pixel
	for b in range(len(cf.zbins)-1):
		zmin, zmax = cf.zbins[b:b+2]
		zmask = (cat[f'{group}/{cf.zcol}'][:] >= zmin) * (cat[f'{group}/{cf.zcol}'][:] < zmax)
		_, ngal = countsInPixels(cat[f'{group}/ra'][zmask], cat[f'{group}/dec'][zmask], cf.nside_lo, cf.nside_hi, mask.vpix_nest, return_vals=True)

		ngal_maps[labels[b]][mask.vpix_nest] = np.asarray(ngal, dtype='f8')

	return ngal_maps


def makeDensityMaps(ngal_maps, mask):
	'''
	Uses galaxy counts maps along with a mask to create maps of delta_g (galaxy overdensity) in
	each tomographic bin. Assumes one mask for all bins.

	Parameters
	----------
	ngal_maps: recarray
		recarray for which each entry is a galaxy counts map in each of the specified bands.

	mask: HealSparseMap
		Survey mask, in which each pixel is assumed to effectively be the detection fraction.

	Returns
	-------
	deltag_maps: recarray
		recarray for which each entry is a galaxy overdensity map for a given redshift bin.
	'''

	print('Creating delta_g maps...')
	#identify pixels in the mask above the specified threshold
	### NOTE: due to the way in which the maps and masks are constructed, they should all
	### have the same valid_pixels

	'''vpix_mask, ra, dec = mask.valid_pixels_pos(return_pixels=True)
	keep = mask[vpix_mask] > cf.weight_thresh
	#pix_remove = vpix_mask[mask[vpix_mask] <= cf.weight_thresh]
	ra = ra[keep]
	dec = dec[keep]
	pix_keep = vpix_mask[keep]'''

	#initialise a recarray to contain the maps for each band
	labels = [f'delta_{i}' for i in range(len(cf.zbins)-1)]
	#make a copy of the inputted recarray
	deltag_maps, *_ = initialiseRecMap(cf.nside_lo, cf.nside_hi, labels, pixels=mask.vpix_nest, dtypes='f8')
	for key in deltag_maps.dtype.names:
		#calculate the mean galaxy counts and mean weight for pixels being kept
		#print(np.array(list(set(pix_keep)-set(ngal_maps[key].valid_pixels))))
		ngal_key = key.replace('delta', 'ngal')
		#deltag_maps[key][pix_remove] = 0.
		mu_n = np.mean(ngal_maps[ngal_key][mask.vpix_nest])
		mu_w = np.mean(mask.mask[mask.vpix])
		#update the density map with delta_g values
		deltag_maps[key][mask.vpix_nest] = (ngal_maps[ngal_key][mask.vpix_nest] / (mask.mask[mask.vpix] * mu_n / mu_w)) - 1.

	return deltag_maps

		




#######################################################
###############    START OF SCRIPT    #################
#######################################################

#cycle through each of the fields
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))
	#output directory for this field
	OUT = cf.PATH_OUT + fd
	#load the fully cleaned galaxy catalogue for this field
	cat_main = h5py.File(f'{OUT}/{cf.cat_main}', 'r')
	#load the survey mask
	survey_mask = MaskData(f'{OUT}/{cf.survey_mask}', 
							mask_thresh=cf.weight_thresh,
							smooth=True,
							fwhm_arcmin=60,
							smoothed_thresh=0.4,
							smooth_file=f'{OUT}/{cf.survey_mask[:-4]}_smoothed.hsp')
	print(f'{OUT}/{cf.survey_mask[:-4]}_smoothed.hsp')

	#make the galaxy count maps in each redsift bin and store in a single recarray
	ngal_maps = makeNgalMaps(cat_main, mask=survey_mask, group='photometry')
	#write to a file
	ngal_maps.write(f'{OUT}/{cf.ngal_maps}', clobber=True)

	#make the delta_g maps in each redsift bin and store in a single recarray
	deltag_maps = makeDensityMaps(ngal_maps, survey_mask)
	#write to a file
	deltag_maps.write(f'{OUT}/{cf.deltag_maps}', clobber=True)