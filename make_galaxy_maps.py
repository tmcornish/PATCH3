#####################################################################################################
# - Creates maps of galaxy counts and galaxy density in each tomographic bin.
#####################################################################################################

import sys
from configuration import PipelineConfig as PC
import healsparse as hsp
import numpy as np
from output_utils import colour_string
import map_utils as mu
import h5py

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='makeGalaxyMaps')


###################
#### FUNCTIONS ####
###################

def makeNgalMaps(cat, footprint):
	'''
	Creates galaxy count maps for each tomographic bin.

	Parameters
	----------
	cat: h5py.File object
		Catalogue containing (at least): RAs, Decs, and dust attenuation values at each of these
		coordinates in each of the specified bands.

	group: str
		Group within which the relevant data are expected to reside.

	footprint: HealSparseMap
		HealSparse map identifying the pixels in which sources reside.

	Returns
	-------
	ngal_maps: recarray
		recarray for which each entry is a galaxy counts map for a given redshift bins.
	'''

	print('Creating galaxy count maps...')
	#initialise a recarray to contain the maps for each band
	labels = [k for k in cf.samples]
	ngal_maps, _, _ = mu.initialiseRecMap(cf.nside_lo, cf.nside_hi, labels, pixels=footprint.valid_pixels, dtypes='f8')

	#cycle through the sample definitions and count the galaxies in each pixel
	for k in cf.samples:
		sel = cat[k][:]
		_, ngal = mu.countsInPixels(cat[f'ra'][sel], cat[f'dec'][sel], cf.nside_lo, cf.nside_hi, footprint.valid_pixels, return_vals=True)
		ngal_maps[k][footprint.valid_pixels] = np.asarray(ngal, dtype='f8')

	return ngal_maps


def makeDensityMaps(ngal_maps, mask):
	'''
	Uses galaxy counts maps along with a mask to create maps of delta_g (galaxy overdensity) in
	each tomographic bin. Assumes one mask for all bins.

	Parameters
	----------
	ngal_maps: recarray
		recarray for which each entry is a galaxy counts map in each of the specified bands.

	mask: MaskData object
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

	#initialise a recarray to contain the maps for each band
	labels = [k for k in cf.samples]
	#make a copy of the inputted recarray
	deltag_maps, *_ = mu.initialiseRecMap(cf.nside_lo, cf.nside_hi, labels, pixels=mask.vpix_nest, dtypes='f8')
	for k in cf.samples:
		#calculate the mean galaxy counts and mean weight for pixels being kept
		mu_n = np.mean(ngal_maps[k][mask.vpix_nest])
		mu_w = np.mean(mask.mask_hsp[mask.vpix_nest])
		#update the density map with delta_g values
		deltag_maps[k][mask.vpix_nest] = (ngal_maps[k][mask.vpix_nest] / (mask.mask_hsp[mask.vpix_nest] * mu_n / mu_w)) - 1.

	return deltag_maps

		




#######################################################
###############    START OF SCRIPT    #################
#######################################################

#cycle through each of the fields
for fd in cf.fields:
	print(colour_string(fd.upper(), 'orange'))
	#output directory for this field
	OUT = cf.paths.out + fd
	#load the fully cleaned galaxy catalogue for this field
	cat_main = h5py.File(f'{OUT}/{cf.cats.main}', 'r')['photometry']
	#load the survey footprint and mask
	footprint = hsp.HealSparseMap.read(f'{OUT}/footprint_{cf.nside_hi}.hsp')
	survey_mask = mu.MaskData(f'{OUT}/{cf.maps.survey_mask}')

	#make the galaxy count maps in each redsift bin and store in a single recarray
	ngal_maps = makeNgalMaps(cat_main, footprint)
	#write to a file
	ngal_maps.write(f'{OUT}/{cf.maps.ngal_maps}', clobber=True)

	#make the delta_g maps in each redsift bin and store in a single recarray
	deltag_maps = makeDensityMaps(ngal_maps, survey_mask)
	#write to a file
	deltag_maps.write(f'{OUT}/{cf.maps.deltag_maps}', clobber=True)