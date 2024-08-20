#####################################################################################################
# - Creates the following maps and masks using HealSparse:
#	- dust attenuation in each band
#	- star counts
#	- binary bright object mask
#	- masked fraction map
#	- depth map
#	- survey mask
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
cf = config.makeMapsFromCat
npix = hp.nside2npix(cf.nside_hi)


###################
#### FUNCTIONS ####
###################


def makeFootprint(cat, nside, keep=None):
	'''
	Defines the survey footprint as any pixel in a map of resolution NSIDE within which
	there are sources. Returns a boolean HealSparse map.

	Parameters
	----------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs, Decs.
	
	nside: int
		Determines the resolution of the footprint. (Left as an argument rather
		than using the value in the config file for flexibility.)
		
	keep: array or None
		Boolean array with length equal to the input catalogue, specifying
		which sources to keep. If set to None, will keep all sources.

	Returns
	-------
	footprint: HealSparseMap
		HealSparse map containing True at all pixels in which sources exist.
	'''

	print(f'Determining survey footprint at NSIDE={nside}...')

	#see if sources are to be removed when considering the footprint
	if keep is None:
		keep = np.ones_like(cat[f'ra'][:], dtype=bool)
	#get the pixel IDs corresponding to each source
	ipix_all = hp.ang2pix(nside, cat[f'ra'][keep], cat[f'dec'][keep], lonlat=True, nest=True)
	#identify pixels where sources exist
	nall = np.bincount(ipix_all, minlength=npix)
	footprint = nall > 0
	#set up empty HealSparse boolean map
	footprint_hsp = hsp.HealSparseMap.make_empty(cf.nside_lo, nside, bool)
	#fill the occupied pixels with True
	footprint_hsp[np.where(footprint)[0]] = True

	return footprint_hsp


def makeDustMap(cat, band='i'):
	'''
	Creates dust maps for each of the specified bands.

	Parameters
	---------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs, Decs, and dust attenuation values at each of these
		coordinates in each of the specified bands.

	Returns
	-------
	dust_maps: recarray
		recarray for which each entry is a dust attenuation map in each of the specified bands.
	'''

	print(f'Creating dust map (band {band})...')

	#cycle through the bands and calculate the mean dust attenuation in each pixel
	dust_map, _ = createMeanStdMap(cat[f'ra'][:], cat[f'dec'][:], cat[f'a_{band}'][:], cf.nside_lo, cf.nside_hi)

	return dust_map


def makeBOMask(cat):
	'''
	Creates bright object mask.

	Parameters
	---------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs, Decs.

	Returns
	-------
	bo_mask: HealSparseMap
		HealSparse map containing 0s at the masked positions and 1 everywhere else.
	'''

	print('Creating bright object mask...')
	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'ra'][:]
	dec = cat[f'dec'][:]
	#get the columns containing bright-object flags
	flags = [cat[f'{flag_col}'][:] for flag_col in cf.bo_flags]

	bo_mask = createMask(ra, dec, flags, cf.nside_lo, cf.nside_hi)

	return bo_mask


def makeMaskedFrac(cat, nside):
	'''
	Creates a map showing the fraction of each pixel that is maske by the bright object criteria.

	Parameters
	---------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs, Decs.
	
	nside: int
		Determines the resolution of the map. (Left as an argument rather
		than using the value in the config file for flexibility.)
		
	Returns
	-------
	mf_map: HealSparseMap
		HealSparse map containing the fraction of each pixel that is masked.
	'''

	print('Creating masked fraction map...')
	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'ra'][:]
	dec = cat[f'dec'][:]
	#get the columns containing bright-object flags
	flags = [cat[f'{flag_col}'][:] for flag_col in cf.bo_flags]
	#add all the masks together
	add = lambda x,y : x+y
	flagged = reduce(add, flags)

	#get the valid pixels at the desired resolution
	vpix = makeFootprint(cat, nside).valid_pixels

	#create counts map of all sources
	Ntotal_map = countsInPixels(ra, dec, cf.nside_lo, nside, vpix)
	#create counts map of flagged sources
	Nflagged_map = countsInPixels(ra[flagged], dec[flagged], cf.nside_lo, nside, vpix)
	#calculate the fraction of masked sources in each pixel
	mf_map = hsp.HealSparseMap.make_empty(cf.nside_lo, nside, dtype=np.float64)
	mf_map[vpix] = Nflagged_map[vpix] / Ntotal_map[vpix]

	return mf_map


def makeStarMap(cat):
	'''
	Creates map showing the number of stars in each pixel.

	Parameters
	---------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs and Decs of each star detected in the field.

	Returns
	-------
	star_map: HealSparseMap
		HealSparse map containing the number of stars in each pixel.
	'''
	print('Creating star counts map...')
	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'ra'][:]
	dec = cat[f'dec'][:]

	#count the stars in each pixel
	star_map = pixelCountsFromCoords(ra, dec, cf.nside_lo, cf.nside_hi)

	return star_map


def makeDepthMap(cat, stars_only=True):
	'''
	Creates map showing the number of stars in each pixel.

	Parameters
	---------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs and Decs of each star detected in the field.

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
	snrs = cat[f'{cf.band}_cmodel_flux'][:] / cat[f'{cf.band}_cmodel_fluxerr'][:]
	'''
	#if told to use stars only, use flags to identify which sources are stars
	if stars_only:
		try:
			star_mask = cat[f'is_star'][:]
		except KeyError:
			print(colour_string('Error: '), f'No dataset "is_star" found. Using all sources instead.')
			star_mask = np.full(len(cat[f'ra']), True)
	else:
		star_mask = np.full(len(cat[f'ra']), True)

	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'ra'][star_mask]
	dec = cat[f'dec'][star_mask]


	#retrieve the flux error in the primary band for each source
	fluxerr = cat[f'{cf.band}_cmodel_fluxerr'][star_mask]
	#retrieve the SNR threshold
	snr_thresh = int(cf.sn_pri)

	#create a map containing the mean fluxerr multiplied by the SNR threshold
	depth_map, _ = createMeanStdMap(ra, dec, snr_thresh * fluxerr, cf.nside_lo, cf.nside_hi)
	vpix = depth_map.valid_pixels
	#fluxes and associated errors are given in nJy - convert to AB mags
	depth_map[vpix] = -2.5 * np.log10(depth_map[vpix] * 10. ** (-9.)) + 8.9

	return depth_map


def makeSurveyMask(cat, depth_map=None):
	'''
	Creates a binary mask from the masked fraction map by applying a threshold for the
	masked fraction.

	Parameters
	----------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs and Decs of each star detected in the field.

	depth_map: HealSparseMap or None
		(Optional) Map of the survey depth. If provided, will additionally set to 0 the
		weight of any pixels below the depth threshold specified in the config file. 

	Returns
	-------
	mask: HealSparseMap
		Binary mask.
	'''
	print('Creating survey mask...')

	#create the masked fraction map using flags from the catalogue
	if cf.highres_first:
		# if initially making the mask at high resolution the process becomes 
		# somewhat convoluted, as one needs to identify pixels with 0 galaxies
		# that also lie within the footprint...
		nside = cf.nside_mask
		# first the footprint must be defined at the resolution of the maps, 
		# then upgraded to high resolution (HealSparse format)
		fp = makeFootprint(cat, cf.nside_hi).upgrade(nside)
		# then the masked fraction must be computed at high resolution (HealSparse format)
		mf_map = makeMaskedFrac(cat, nside)
		# initially define the mask in HealPy format at high resolution
		mask = np.zeros(hp.nside2npix(nside))
		# fill valid pixels using the masked fraction map
		mask[mf_map.valid_pixels] = 1. - mf_map[mf_map.valid_pixels]
		# identify valid pixels in the footprint
		vpix = fp.valid_pixels
		# construct a mask for identifying pixels with sources in the footprint
		fp_bin = np.zeros(hp.nside2npix(nside), dtype=bool)
		fp_bin[vpix] = True
		# set all other pixels in the mask equal to UNSEEN
		mask[~fp_bin] = hp.UNSEEN
		#convert this to a HealSparse map and degrade to the final resolution
		mask = hsp.HealSparseMap(nside_coverage=32, healpix_map=mask, nest=True)
		mask = mask.degrade(cf.nside_hi)
	else:
		#otherwise, define the mask simply as 1-masked_fraction
		mask = makeMaskedFrac(cat, cf.nside_hi)
		vpix = mask.valid_pixels
		mask[vpix] = 1. - mask[vpix]

	#mask pixels below the depth threshold if a depth map is provided
	if depth_map is not None:
		vpix_dm = depth_map.valid_pixels
		vpix_shallow = vpix_dm[depth_map[vpix_dm] < cf.depth_cut]
		mask[vpix_shallow] = 0.

	return mask


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#cycle through each of the fields
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))
	#output directory for this field
	OUT = cf.PATH_OUT + fd
	#path for systematics maps (create if it doesn't exist already)
	PATH_SYST = OUT + '/systmaps'
	if not os.path.exists(PATH_SYST):
		os.system(f'mkdir -p {PATH_SYST}')

	#load the basic and fully cleaned galaxy/star catalogues for this field
	cat_basic = h5py.File(f'{OUT}/{cf.cat_basic}', 'r')['photometry']
	cat_main = h5py.File(f'{OUT}/{cf.cat_main}', 'r')['photometry']
	cat_stars = h5py.File(f'{OUT}/{cf.cat_stars}', 'r')['photometry']

	#make the footprint for the current field
	footprint = makeFootprint(cat_basic, cf.nside_hi, keep=None)
	#write to a file
	footprint.write(f'{OUT}/{cf.footprint}', clobber=True)
	#retrieve the IDs of the occupied pixels
	vpix = footprint.valid_pixels

	for b,dm in zip(cf.bands, cf.dustmaps):
		#make the dust maps in the current band
		dust_map = makeDustMap(cat_basic, band=b)
		#write to a file
		dust_map.write(f'{PATH_SYST}/{dm}', clobber=True)

	#make the bright object mask
	bo_mask = makeBOMask(cat_basic)
	#write to a file
	bo_mask.write(f'{OUT}/{cf.bo_mask}', clobber=True)
	healsparseToHDF(bo_mask, f'{OUT}/{cf.bo_mask[:-4]}.hdf5', group='maps/mask')

	'''#make the masked fraction map
	mf_map = makeMaskedFrac(cat_basic)
	#write to a file
	mf_map.write(f'{OUT}/{cf.masked_frac}', clobber=True)'''

	#make the star counts map
	star_map = makeStarMap(cat_stars)
	#write to a file
	star_map.write(f'{PATH_SYST}/{cf.star_map}', clobber=True)

	#make the depth map
	depth_map = makeDepthMap(cat_basic, stars_only=True)
	#write to a file
	depth_map.write(f'{OUT}/{cf.depth_map}', clobber=True)

	#make a survey mask by applying a masked fraction threshold
	survey_mask = makeSurveyMask(cat_basic, depth_map=depth_map)
	#write to a file
	survey_mask.write(f'{OUT}/{cf.survey_mask}', clobber=True)
	#calculate the area above the mask threshold and the fractional sky coverage
	A_unmasked, f_sky = maskAreaSkyCoverage(survey_mask, thresh=cf.weight_thresh)
	mask_meta = {'area' : A_unmasked, 'f_sky': f_sky}
	healsparseToHDF(survey_mask, f'{OUT}/{cf.survey_mask[:-4]}.hdf5', group='maps/mask', metadata=mask_meta)