#####################################################################################################
# - Creates the following maps and masks using HealSparse:
#	- dust attenuation in each band
#	- star counts
#	- binary bright object mask
#	- masked fraction map
#	- depth map
#####################################################################################################

import os
import config
import healpy as hp
import healsparse as hsp
from astropy.table import Table
import numpy as np
from output_utils import colour_string
from map_utils import *
import h5py
from functools import reduce

### SETTINGS ###
cf = config.makeMaps


###################
#### FUNCTIONS ####
###################


def makeDustMap(cat, group=''):
	'''
	Creates dust maps for each of the specified bands.

	Paameters
	---------
	cat: h5py.File object
		Catalogue containing (at least): RAs, Decs, and dust attenuation values at each of these
		coordinates in each of the specified bands.

	group: str
		Group within which the relevant data are expected to reside.

	Returns
	-------
	dust_maps: recarray
		recarray for which each entry is a dust attenuation map in each of the specified bands.
	'''

	print('Creating dust maps...')
	#initialise a recarray to contain the maps for each band
	dust_maps, px_data, px_data_u = initialiseRecMap(cf.nside_lo, cf.nside_hi, cat[f'{group}/ra'][:], cat[f'{group}/dec'][:], cf.bands, dtypes='f8', primary=cf.band)

	#cycle through the bands and calculate the mean dust attenuation in each pixel
	for b in cf.bands:
		a_means, _ = pixelMeanStd(cat[f'{group}/a_{b}'][:], px_data)
		dust_maps[b][px_data_u] = a_means

	return dust_maps


def makeBOMask(cat, group=''):
	'''
	Creates bright object mask.

	Paameters
	---------
	cat: h5py.File object
		Catalogue containing (at least): RAs, Decs.

	group: str
		Group within which the relevant data are expected to reside.

	Returns
	-------
	bo_mask: HealSparseMap
		HealSparse map containing 0s at the masked positions and 1 everywhere else.
	'''

	print('Creating bright object mask...')
	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'{group}/ra'][:]
	dec = cat[f'{group}/dec'][:]
	#get the columns containing bright-object flags
	flags = [cat[f'{group}/{flag_col}'][:] for flag_col in cf.bo_flags]

	bo_mask = createMask(ra, dec, flags, cf.nside_lo, cf.nside_hi)

	return bo_mask


def makeMaskedFrac(cat, group=''):
	'''
	Creates a map showing the fraction of each pixel that is maske by the bright object criteria.

	Paameters
	---------
	cat: h5py.File object
		Catalogue containing (at least): RAs, Decs.

	group: str
		Group within which the relevant data are expected to reside.

	Returns
	-------
	mf_map: HealSparseMap
		HealSparse map containing the fraction of each pixel that is masked.
	'''

	print('Creating masked fraction map...')
	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'{group}/ra'][:]
	dec = cat[f'{group}/dec'][:]
	#get the columns containing bright-object flags
	flags = [cat[f'{group}/{flag_col}'][:] for flag_col in cf.bo_flags]
	#multiply all the masks together
	mult = lambda x,y : x*y
	flagged = reduce(mult, flags)

	#create counts map of all sources
	Ntotal_map = pixelCountsFromCoords(ra, dec, cf.nside_lo, cf.nside_mask)
	vpix_total = Ntotal_map.valid_pixels
	#create counts map of flagged sources
	Nflagged_map = pixelCountsFromCoords(ra[flagged], dec[flagged], cf.nside_lo, cf.nside_mask)
	vpix_flagged = Nflagged_map.valid_pixels
	#identify pixels that only contain unflagged sources
	vpix_unflagged_only = np.array(list(set(vpix_total) - set(vpix_flagged)))
	Nflagged_map[vpix_unflagged_only] = np.zeros(len(vpix_unflagged_only), dtype=np.int32)
	#calculate the fraction of masked sources in each pixel
	mf_map = hsp.HealSparseMap.make_empty(cf.nside_lo, cf.nside_mask, dtype=np.float64)
	mf_map[vpix_total] = Nflagged_map[vpix_total] / Ntotal_map[vpix_total]
	#degrade the resolution to the resolution of the bright object mask
	mf_map = mf_map.degrade(cf.nside_hi, reduction='mean')

	return mf_map


def makeStarMap(cat, group=''):
	'''
	Creates map showing the number of stars in each pixel.

	Paameters
	---------
	cat: h5py.File object
		Catalogue containing (at least): RAs and Decs of each star detected in the field.

	group: str
		Group within which the relevant data are expected to reside.

	Returns
	-------
	star_map: HealSparseMap
		HealSparse map containing the number of stars in each pixel.
	'''
	print('Creating star counts map...')
	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'{group}/ra'][:]
	dec = cat[f'{group}/dec'][:]

	#count the stars in each pixel
	star_map = pixelCountsFromCoords(ra, dec, cf.nside_lo, cf.nside_hi)

	return star_map


def makeDepthMap(cat, group='', stars_only=True):
	'''
	Creates map showing the number of stars in each pixel.

	Paameters
	---------
	cat: h5py.File object
		Catalogue containing (at least): RAs and Decs of each star detected in the field.

	group: str
		Group within which the relevant data are expected to reside.

	stars_only: bool
		Whether to just use stars for the calculation of the depth.

	Returns
	-------
	depth_map: HealSparseMap
		HealSparse map containing the average N-sigma depth in each pixel, where N is the
		SNR threshold of the primary band (set in the config file).
	'''
	print('Creating depth map...')
	'''
	#begin by calculating SNR of each source in the primary band
	snrs = cat[f'{group}/{cf.band}_cmodel_flux'][:] / cat[f'{group}/{cf.band}_cmodel_fluxerr'][:]
	'''
	#if told to use stars only, use flags to identify which sources are stars
	if stars_only:
		try:
			star_mask = cat[f'{group}/is_star'][:]
		except KeyError:
			print(colour_string('Error: '), f'No dataset "{group}/is_star" found. Using all sources instead.')
			star_mask = np.full(len(cat[f'{group}/ra']), True)
	else:
		star_mask = np.full(len(cat[f'{group}/ra']), True)

	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'{group}/ra'][star_mask]
	dec = cat[f'{group}/dec'][star_mask]


	#retrieve the flux error in the primary band for each source
	fluxerr = cat[f'{group}/{cf.band}_cmodel_fluxerr'][star_mask]
	#retrieve the SNR threshold
	snr_thresh = int(cf.sn_pri)

	#create a map containing the mean fluxerr multiplied by the SNR threshold
	depth_map, _ = createMeanStdMap(ra, dec, snr_thresh * fluxerr, cf.nside_lo, cf.nside_hi)
	vpix = depth_map.valid_pixels
	#fluxes and associated errors are given in nJy - convert to AB mags
	depth_map[vpix] = -2.5 * np.log10(depth_map[vpix] * 10. ** (-9.)) + 8.9

	return depth_map


def makeSurveyMask(mf_map):
	'''
	Creates a binary mask from the masked fraction map by applying a threshold for the
	masked fraction.

	Parameters
	----------
	mf_map: HealSparseMap
		Masked fraction map.

	Returns
	-------
	mask: HealSparseMap
		Binary mask.
	'''
	print('Creating survey mask...')
	#make an empty map with the same properties as the masked fraction map
	mask = hsp.HealSparseMap.make_empty(mf_map.nside_coverage, mf_map.nside_sparse, dtype=np.float64)
	#identify valid pixels
	vpix = mf_map.valid_pixels
	'''
	#identify pixels above and below the masked fraction threshold
	above = mf_map[vpix] >= thresh
	below = mf_map[vpix] < thresh
	#populate the binary mask with 0s and 1s accordingly
	mask[vpix[above]] = 0.
	mask[vpix[below]] = 1.
	'''
	#mask[vpix] = 1. - mf_map[vpix]

	#fill the survey mask with the masked fraction at the relevant pixels
	mask[vpix] = mf_map[vpix]

	return mask


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#cycle through each of the fields
for fd in cf.fields:
	#output directory for this field
	OUT = cf.PATH_OUT + fd
	#load the basic and fully cleaned galaxy/star catalogues for this field
	cat_basic = h5py.File(f'{OUT}/{cf.cat_basic}', 'r')
	cat_main = h5py.File(f'{OUT}/{cf.cat_main}', 'r')
	cat_stars = h5py.File(f'{OUT}/{cf.cat_stars}', 'r')

	#make the dust maps in each band and store in a single recarray
	dust_maps = makeDustMap(cat_basic, group='photometry')
	#write to a file
	dust_maps.write(f'{OUT}/{cf.dustmaps}', clobber=True)

	#make the bright object mask
	bo_mask = makeBOMask(cat_basic, group='photometry')
	#write to a file
	bo_mask.write(f'{OUT}/{cf.bo_mask}', clobber=True)
	healsparseToHDF(bo_mask, f'{OUT}/{cf.bo_mask[:-4]}.hdf5', group='maps/mask')

	#make the masked fraction map
	mf_map = makeMaskedFrac(cat_basic, group='photometry')
	#write to a file
	mf_map.write(f'{OUT}/{cf.masked_frac}', clobber=True)

	#make a survey mask by applying a masked fraction threshold
	survey_mask = makeSurveyMask(mf_map, thresh=0.5)
	#write to a file
	survey_mask.write(f'{OUT}/{cf.survey_mask}', clobber=True)
	healsparseToHDF(survey_mask, f'{OUT}/{cf.survey_mask[:-4]}.hdf5', group='maps/mask')

	#make the star counts map
	star_map = makeStarMap(cat_stars, group='photometry')
	#write to a file
	star_map.write(f'{OUT}/{cf.star_map}', clobber=True)

	#make the depth map
	depth_map = makeDepthMap(cat_basic, group='photometry', stars_only=True)
	#write to a file
	depth_map.write(f'{OUT}/{cf.depth_map}', clobber=True)