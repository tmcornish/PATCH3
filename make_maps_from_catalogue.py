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
import flags as fl

### SETTINGS ###
cf = config.makeMapsFromCat
npix = hp.nside2npix(cf.nside_hi)


###################
#### FUNCTIONS ####
###################


def makeFootprint(cat, nside):
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

	Returns
	-------
	footprint: HealSparseMap
		HealSparse map containing True at all pixels in which sources exist.
	'''

	print(f'Determining survey footprint at NSIDE={nside}...')

	#get the pixel IDs corresponding to each source
	ipix_all = hp.ang2pix(nside, cat[f'ra'][:], cat[f'dec'][:], lonlat=True, nest=True)
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
	flags = fl.get_flags(cf.band, cf.bands, types=['brightstar'], incl_channelstop=cf.incl_channelstop)
	flags = [cat[flag][:] for flag in flags] 

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
	ra = cat['ra'][:]
	dec = cat['dec'][:]
	#get the columns containing the flags to be incorporated in the mask
	flags = fl.get_flags(cf.band, cf.bands, types=cf.flags_to_mask, incl_channelstop=cf.incl_channelstop)
	#combine the flags
	flagged = fl.combine_flags(cat, flags)

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


def makeStarMap(cat, footprint):
	'''
	Creates map showing the number of stars in each pixel.

	Parameters
	---------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs and Decs of each star detected in the field.

	footprint: HealSparseMap or None
		If provided, must be a Healpsarse boolean map identifying pixels belonging to the
		survey footprint. If not provided, will make from scratch using the catalogue.
		
	Returns
	-------
	star_map: HealSparseMap
		HealSparse map containing the number of stars in each pixel.
	'''
	print('Creating star counts map...')
	#get the RAs and Decs of all sources in the catalogue
	ra = cat[f'ra'][:]
	dec = cat[f'dec'][:]

	#identify pixels in the survey footprint
	vpix = footprint.valid_pixels
	#count the stars in each pixel
	star_map = countsInPixels(ra, dec, cf.nside_lo, cf.nside_hi, vpix)
	#pixelCountsFromCoords(ra, dec, cf.nside_lo, cf.nside_hi)

	return star_map


def makeDepthMap(cat, stars_only=True, min_sources=0, footprint=None):
	'''
	Creates map showing the number of stars in each pixel.

	Parameters
	---------
	cat: h5py Dataset or Group
		Catalogue containing (at least): RAs and Decs of each star detected in the field.

	stars_only: bool
		Whether to just use stars for the calculation of the depth.
	
	min_sources: int
		Minimum number of sources required in a pixel for depth to be calculated.	
	
	footprint: HealSparseMap or None
		If provided, must be a Healpsarse boolean map identifying pixels belonging to the
		survey footprint. If not provided, will only use pixels in which there are sources
		in the 

	Returns
	-------
	depth_map: HealSparseMap
		HealSparse map containing the average N-sigma depth in each pixel, where N is the
		SNR threshold of the primary band (set in the config file).
	'''
	print('Creating depth map...')

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

	#if a minimum number of sources is required, use interpolation to fill in the
	#values for pixels with fewer sources
	if min_sources > 0:
		if footprint:
			vpix_fp = footprint.valid_pixels
			counts = countsInPixels(ra, dec, cf.nside_lo, cf.nside_hi, vpix_fp)
		else:
			counts = pixelCountsFromCoords(ra, dec, cf.nside_lo, cf.nside_hi)
			vpix_fp = counts.valid_pixels 
		#identify valid pixels with fewer sources than the limit 
		pix_few = vpix_fp[counts[vpix_fp] < min_sources]
		#get the coordinates of these pixels
		ra_few, dec_few = hp.pix2ang(cf.nside_hi, pix_few, nest=True, lonlat=True)
		#get the values of these pixels through nearest-neighbour interpolation
		pix_new_vals = depth_map.interpolate_pos(ra_few, dec_few, lonlat=True, allow_partial=True)
		#update these pixels in the map
		depth_map[pix_few] = pix_new_vals

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
		mask = makeMaskedFrac(cat, nside)
		vpix = mask.valid_pixels
		# fill valid pixels using the masked fraction map
		mask[vpix] = 1. - mask[vpix]
		# identify valid pixels in the footprint
		vpix_fp = fp.valid_pixels
		# identify which of these pixels were not valid in the masked fraction map
		vpix_empty = list(set(vpix_fp) - set(vpix))
		#set these pixels equal to 0
		mask[vpix_empty] = 0.
		#degrade the mask to the final resolution
		mask = mask.degrade(cf.nside_hi)
		vpix = mask.valid_pixels
	else:
		#otherwise, define the mask simply as 1-masked_fraction
		mask = makeMaskedFrac(cat, cf.nside_hi)
		vpix = mask.valid_pixels
		mask[vpix] = 1. - mask[vpix]
	
	#set any pixels below the mask threshold equal to UNSEEN
	mask[vpix[mask[vpix] < cf.weight_thresh]] = hp.UNSEEN

	#mask pixels below the depth threshold if a depth map is provided
	if depth_map is not None:
		vpix_dm = depth_map.valid_pixels
		vpix_shallow = vpix_dm[depth_map[vpix_dm] < cf.depth_cut]
		mask[vpix_shallow] = hp.UNSEEN

		#create a full-sky binary version of the depth map
		depth_bin = depth_map[:] >= cf.depth_cut
		#smooth it with a Guassian kernel
		depth_bin = hp.smoothing(depth_bin, np.radians(cf.r_smooth), nest=True)
		#identify valid pixels below 0.5 in the smoothed version
		vpix_shallow = np.where(depth_bin[vpix] < 0.7)[0]
		mask[vpix[vpix_shallow]] = hp.UNSEEN

	if cf.use_nexp_maps:
		#set up a list of maps to be combined into a binary map
		nexp_bin = []
		#load the N_exp maps for each band
		for b in cf.bands:
			nexp = hsp.HealSparseMap.read(f'{PATH_SYST}/decasu_{cf.nside_hi}_{b}_nexp_sum.hsp')
			#set all pixels with n_exp > 1 equal to 1
			nexp[nexp.valid_pixels] = 1
			nexp_bin.append(nexp)
		#identify all pixels with exposures in all frames
		nexp_bin = hsp.operations.and_union(nexp_bin)
		#convert this to a healpix map and set all UNSEEN to 0
		nexp_bin = nexp_bin.generate_healpix_map(nest=True)
		nexp_bin[nexp_bin < 0] = 0
		#smooth with a Gaussian kernel
		nexp_bin = hp.smoothing(nexp_bin, np.radians(cf.r_smooth), nest=True)
		#identify valid pixels below 0.6
		vpix_noexp = np.where(nexp_bin[vpix] < 0.7)[0]
		mask[vpix[vpix_noexp]] = hp.UNSEEN



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
	footprint = makeFootprint(cat_basic, cf.nside_hi)
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

	#make the masked fraction map
	mf_map = makeMaskedFrac(cat_basic, cf.nside_hi)
	#write to a file
	mf_map.write(f'{OUT}/{cf.masked_frac}', clobber=True)

	#make the star counts map
	star_map = makeStarMap(cat_stars, footprint)
	#write to a file
	star_map.write(f'{PATH_SYST}/{cf.star_map}', clobber=True)

	#make the depth map
	depth_map = makeDepthMap(cat_basic, stars_only=cf.stars_for_depth, min_sources=cf.min_sources, footprint=footprint)
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