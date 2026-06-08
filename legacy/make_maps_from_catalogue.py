##############################################################################
# - Creates the following maps and masks using HealSparse:
#    - dust attenuation in each band
#    - star counts
#    - binary bright object mask
#    - masked fraction map
#    - depth map
#    - survey mask
##############################################################################

import os
import sys
from configuration import PipelineConfig as PC
import healpy as hp
import healsparse as hsp
import numpy as np
from output_utils import colour_string
import map_utils as mu
import h5py

# SETTINGS #
config_file = sys.argv[1]
cf = PC(config_file, stage='makeMapsFromCat')
npix = hp.nside2npix(cf.nside_hi)

###################
#    FUNCTIONS    #
###################


def makeFootprint(cat, nside):
    '''
    Defines the survey footprint as any pixel in a map of resolution NSIDE
    within which there are sources. Returns a boolean HealSparse map.

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

    # Get the pixel IDs corresponding to each source
    ipix_all = hp.ang2pix(nside, cat['ra'][:], cat['dec'][:],
                          lonlat=True, nest=True)
    # Identify pixels where sources exist
    nall = np.bincount(ipix_all, minlength=npix)
    footprint = nall > 0
    # Set up empty HealSparse boolean map
    footprint_hsp = hsp.HealSparseMap.make_empty(cf.nside_lo, nside, bool)
    # Fill the occupied pixels with True
    footprint_hsp[np.where(footprint)[0]] = True

    return footprint_hsp


def makeDustMap(cat, band='i'):
    '''
    Creates dust maps for each of the specified bands.

    Parameters
    ---------
    cat: h5py Dataset or Group
        Catalogue containing (at least): RAs, Decs, and dust attenuation
        values at each of these coordinates in each of the specified bands.

    Returns
    -------
    dust_maps: recarray
        recarray for which each entry is a dust attenuation map in each of
        the specified bands.
    '''

    print(f'Creating dust map (band {band})...')

    # Cycle through bands and calculate mean dust attenuation in each pixel
    dust_map, _ = mu.createMeanStdMap(cat['ra'][:], cat['dec'][:],
                                      cat[f'a_{band}'][:], cf.nside_lo,
                                      cf.nside_hi)

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
        HealSparse map containing 0s at the masked positions and 1
        everywhere else.
    '''

    print('Creating bright object mask...')
    # Get the RAs and Decs of all sources in the catalogue
    ra = cat['ra'][:]
    dec = cat['dec'][:]
    # Get the columns containing bright-object flags
    flags = cf.flags.brightstar
    flags = [cat[flag][:] for flag in flags]

    bo_mask = mu.createMask(ra, dec, flags, cf.nside_lo, cf.nside_hi)

    return bo_mask


def makeMaskedFrac(cat, nside):
    '''
    Creates a map showing the fraction of each pixel that is maske by the
    bright object criteria.

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
    # Get the RAs and Decs of all sources in the catalogue
    ra = cat['ra'][:]
    dec = cat['dec'][:]
    # Get the columns containing the flags to be incorporated in the mask
    flags = []
    for fl in cf.flags_to_mask:
        flags.extend(cf.flags[fl])
    # Combine the flags (OR operation)
    flagged = cf.combine_flags(cat, flags)

    # Get the valid pixels at the desired resolution
    vpix = makeFootprint(cat, nside).valid_pixels

    # Create counts map of all sources
    Ntotal_map = mu.countsInPixels(ra, dec, cf.nside_lo, nside, vpix)
    # Create counts map of flagged sources
    Nflagged_map = mu.countsInPixels(ra[flagged], dec[flagged], cf.nside_lo,
                                     nside, vpix)
    # Calculate the fraction of masked sources in each pixel
    mf_map = hsp.HealSparseMap.make_empty(cf.nside_lo, nside, dtype=np.float64)
    mf_map[vpix] = Nflagged_map[vpix] / Ntotal_map[vpix]

    return mf_map


def makeStarMap(cat, footprint):
    '''
    Creates map showing the number of stars in each pixel.

    Parameters
    ---------
    cat: h5py Dataset or Group
        Catalogue containing (at least): RAs and Decs of each star detected in
        the field.

    footprint: HealSparseMap or None
        If provided, must be a Healpsarse boolean map identifying pixels
        belonging to the survey footprint. If not provided, will make from
        scratch using the catalogue.

    Returns
    -------
    star_map: HealSparseMap
        HealSparse map containing the number of stars in each pixel.
    '''
    print('Creating star counts map...')
    # Get the RAs and Decs of all sources in the catalogue
    ra = cat['ra'][:]
    dec = cat['dec'][:]

    # Identify pixels in the survey footprint
    vpix = footprint.valid_pixels
    # Count the stars in each pixel
    star_map = mu.countsInPixels(ra, dec, cf.nside_lo, cf.nside_hi, vpix)

    return star_map


def makeDepthMap(cat, stars_only=True, min_sources=0, footprint=None):
    '''
    Creates map showing the number of stars in each pixel.

    Parameters
    ---------
    cat: h5py Dataset or Group
        Catalogue containing (at least): RAs and Decs of each star detected in
        the field.

    stars_only: bool
        Whether to just use stars for the calculation of the depth.

    min_sources: int
        Minimum number of sources required in a pixel for depth to be
        calculated.

    footprint: HealSparseMap or None
        If provided, must be a Healpsarse boolean map identifying pixels
        belonging to the survey footprint. If not provided, will only use
        pixels in which there are sources.

    Returns
    -------
    depth_map: HealSparseMap
        HealSparse map containing the average N-sigma depth in each pixel,
        where N is the SNR threshold of the primary band (set in the config
        file).
    '''
    print('Creating depth map...')

    # If told to use stars only, use flags to identify which sources are stars
    if stars_only:
        try:
            star_mask = cat['is_star'][:]
        except KeyError:
            print(colour_string('Error: '),
                  'No dataset "is_star" found. Using all sources instead.')
            star_mask = np.full(len(cat['ra']), True)
    else:
        star_mask = np.full(len(cat['ra']), True)

    # Get the RAs and Decs of all sources in the catalogue
    ra = cat['ra'][star_mask]
    dec = cat['dec'][star_mask]

    # Retrieve the flux error in the primary band for each source
    fluxerr = cat[f'{cf.bands.primary}_cmodel_fluxerr'][star_mask]
    # Retrieve the SNR threshold
    snr_thresh = int(cf.sn_pri)

    # Create a map containing the mean fluxerr multiplied by the SNR threshold
    depth_map, _ = mu.createMeanStdMap(ra, dec, snr_thresh * fluxerr,
                                       cf.nside_lo, cf.nside_hi)

    # If a minimum number of sources is required, use interpolation to fill in
    # the values for pixels with fewer sources
    if min_sources > 0:
        if footprint:
            vpix_fp = footprint.valid_pixels
            counts = mu.countsInPixels(ra, dec, cf.nside_lo,
                                       cf.nside_hi, vpix_fp)
        else:
            counts = mu.pixelCountsFromCoords(ra, dec,
                                              cf.nside_lo, cf.nside_hi)
            vpix_fp = counts.valid_pixels
        # Identify valid pixels with fewer sources than the limit
        pix_few = vpix_fp[counts[vpix_fp] < min_sources]
        # Get the coordinates of these pixels
        ra_few, dec_few = hp.pix2ang(cf.nside_hi, pix_few,
                                     nest=True, lonlat=True)
        # Get values of these pixels through nearest-neighbour interpolation
        pix_new_vals = depth_map.interpolate_pos(ra_few, dec_few, lonlat=True,
                                                 allow_partial=True)
        # Update these pixels in the map
        depth_map[pix_few] = pix_new_vals

    vpix = depth_map.valid_pixels
    # Fluxes and associated errors are given in nJy - convert to AB mags
    depth_map[vpix] = -2.5 * np.log10(depth_map[vpix] * 10. ** (-9.)) + 8.9

    return depth_map


def makeSurveyMask(cat, depth_map=None):
    '''
    Creates a binary mask from the masked fraction map by applying a threshold
    for the masked fraction.

    Parameters
    ----------
    cat: h5py Dataset or Group
        Catalogue containing (at least): RAs and Decs of each star detected in
        the field.

    depth_map: HealSparseMap or None
        (Optional) Map of the survey depth. If provided, will additionally set
        to 0 the weight of any pixels below the depth threshold specified in
        the config file.

    Returns
    -------
    mask: HealSparseMap
        Binary mask.
    '''
    print('Creating survey mask...')

    # Create the masked fraction map using flags from the catalogue
    if cf.highres_first:
        # If initially making the mask at high resolution the process becomes
        # somewhat convoluted, as one needs to identify pixels with 0 galaxies
        # that also lie within the footprint...
        nside = cf.nside_mask
        # First the footprint must be defined at the resolution of the maps,
        # then upgraded to high resolution (HealSparse format)
        fp = makeFootprint(cat, cf.nside_hi).upgrade(nside)
        # Then the masked fraction must be computed at high resolution
        # (HealSparse format)
        mask = makeMaskedFrac(cat, nside)
        vpix = mask.valid_pixels
        # Fill valid pixels using the masked fraction map
        mask[vpix] = 1. - mask[vpix]
        # Identify valid pixels in the footprint
        vpix_fp = fp.valid_pixels
        # Identify which of these pixels weren't valid in masked fraction map
        vpix_empty = list(set(vpix_fp) - set(vpix))
        # Set these pixels equal to 0
        mask[vpix_empty] = 0.
        # Degrade the mask to the final resolution
        mask = mask.degrade(cf.nside_hi)
        vpix = mask.valid_pixels
    else:
        # Otherwise, define the mask simply as 1-masked_fraction
        mask = makeMaskedFrac(cat, cf.nside_hi)
        vpix = mask.valid_pixels
        mask[vpix] = 1. - mask[vpix]

    # Set any pixels below the mask threshold equal to UNSEEN
    mask[vpix[mask[vpix] < cf.weight_thresh]] = hp.UNSEEN

    # Mask pixels below the depth threshold if a depth map is provided
    if depth_map is not None:
        vpix_dm = depth_map.valid_pixels
        vpix_shallow = vpix_dm[depth_map[vpix_dm] < cf.depth_cut]
        mask[vpix_shallow] = hp.UNSEEN

        # Create a full-sky binary version of the depth map
        depth_bin = depth_map[:] >= cf.depth_cut
        # Smooth it with a Guassian kernel
        depth_bin = hp.smoothing(depth_bin, np.radians(cf.r_smooth), nest=True)
        # Identify valid pixels below 0.5 in the smoothed version
        vpix_shallow = np.where(depth_bin[vpix] < 0.7)[0]
        mask[vpix[vpix_shallow]] = hp.UNSEEN

    if cf.use_nexp_maps:
        # Set up a list of maps to be combined into a binary map
        nexp_bin = []
        # Load the N_exp maps for each band
        for b in cf.bands.all:
            nexp = hsp.HealSparseMap.read(f'{PATH_SYST}/decasu_nside'
                                          f'{cf.nside_hi}_{b}_nexp_sum.hsp')
            # Set all pixels with n_exp > 1 equal to 1
            nexp[nexp.valid_pixels] = 1
            nexp_bin.append(nexp)
        # Identify all pixels with exposures in all frames
        nexp_bin = hsp.operations.and_union(nexp_bin)
        # Convert this to a healpix map and set all UNSEEN to 0
        nexp_bin = nexp_bin.generate_healpix_map(nest=True)
        nexp_bin[nexp_bin < 0] = 0
        # Smooth with a Gaussian kernel
        nexp_bin = hp.smoothing(nexp_bin, np.radians(cf.r_smooth), nest=True)
        # Identify valid pixels below 0.6
        vpix_noexp = np.where(nexp_bin[vpix] < 0.7)[0]
        mask[vpix[vpix_noexp]] = hp.UNSEEN

    return mask


#######################################################
#                  START OF SCRIPT                    #
#######################################################

# Cycle through each of the fields
for fd in cf.fields:
    print(colour_string(fd.upper(), 'orange'))
    # Output directory for this field
    OUT = cf.paths.out + fd
    # Path for systematics maps (create if it doesn't exist already)
    PATH_SYST = OUT + '/systmaps'
    if not os.path.exists(PATH_SYST):
        os.system(f'mkdir -p {PATH_SYST}')

    # Load the basic and fully cleaned galaxy/star catalogues for this field
    cat_basic = h5py.File(f'{OUT}/{cf.cats.basic}', 'r')['photometry']
    cat_main = h5py.File(f'{OUT}/{cf.cats.main}', 'r')['photometry']
    cat_stars = h5py.File(f'{OUT}/{cf.cats.stars}', 'r')['photometry']

    # Make the footprint for the current field
    footprint = makeFootprint(cat_basic, cf.nside_hi)
    # Write to a file
    footprint.write(f'{OUT}/{cf.maps.footprint}', clobber=True)
    # Retrieve the IDs of the occupied pixels
    vpix = footprint.valid_pixels

    for b, dm in zip(cf.bands.all, cf.maps.dustmaps):
        # Make the dust maps in the current band
        dust_map = makeDustMap(cat_basic, band=b)
        # Write to a file
        dust_map.write(f'{PATH_SYST}/{dm}', clobber=True)

    # Make the bright object mask
    bo_mask = makeBOMask(cat_basic)
    # Write to a file
    bo_mask.write(f'{OUT}/{cf.maps.bo_mask}', clobber=True)

    # Make the masked fraction map
    mf_map = makeMaskedFrac(cat_basic, cf.nside_hi)
    # Write to a file
    mf_map.write(f'{OUT}/{cf.maps.masked_frac}', clobber=True)

    # Make the star counts map
    star_map = makeStarMap(cat_stars, footprint)
    # Write to a file
    star_map.write(f'{PATH_SYST}/{cf.maps.star_map}', clobber=True)

    # Make the depth map
    depth_map = makeDepthMap(cat_basic, stars_only=cf.stars_for_depth,
                             min_sources=cf.min_sources, footprint=footprint)
    # Write to a file
    depth_map.write(f'{OUT}/{cf.maps.depth_map}', clobber=True)

    # Make a survey mask by applying a masked fraction threshold
    survey_mask = makeSurveyMask(cat_basic, depth_map=depth_map)
    # Write to a file
    survey_mask.write(f'{OUT}/{cf.maps.survey_mask}', clobber=True)
