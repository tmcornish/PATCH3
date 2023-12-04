####################################################################################################
# Functions for the creation of HealPIX/HealSparse maps.
####################################################################################################

import numpy as np
import healpy as hp
import healsparse as hsp


def initialiseRecMap(nside_cover, nside_sparse, ra, dec, labels, dtypes='f8', primary=None, return_pix=True, return_unique=True):
	'''
	Initialises a RecArray of HealSparse maps. Useful e.g. when a single quantity is
	measured for multiple bands. 

	Parameters
	----------
	nside_cover: int
		Resolution of the wider HealSparse map where no data exist.

	nside_sparse: int
		Resolution of the regions of the HealSparse map in which data exist.

	ra: array-like
		RA coordinates at which data exist for the map.

	dec: array-like
		Dec. coordinates at which data exist for the map.

	labels: list
		List of strings to use as labels for each map.

	dtypes: str or type or list
		Describes the type of data with which each map will be filled. If a list,
		must have the same length of labels; otherwise, the same dtype will be used
		for all maps.

	primary: str or None
		Label of the HealSparse map to be used as the primary map.
		If None, uses the first entry of labels.

	return_pix: bool
		Return the pixel positions corresponding to the provided coordinates.

	return_unique: bool
		Return the unique pixel IDs for which data exist.

	Returns
	-------
	all_maps: recarray
		Recarray in which each entry is a HealSparse map, initialised with zeros at the
		pixel positions in which data exist, based on the coordinates provided.

	px_data: array
		Pixel positions corresponding to the provided coordinates.

	px_data_u: array
		Unique pixel IDs for which data exist.
	'''

	#exported HealSparse map will contain dust maps in each band; set up the relevant array
	if hasattr(dtypes, '__len__') and not isinstance(dtypes, str):
		dtype_maps = [(l, d) for l,d in zip(labels, dtypes)]
	else:
		dtype_maps = [(l, dtypes) for l in labels]
	if primary is None:
		primary = l[0]
	all_maps = hsp.HealSparseMap.make_empty(nside_cover, nside_sparse, dtype=dtype_maps, primary=primary)

	#identify pixels in the HealSparse map that contain sources from the catalogue
	px_data = hp.ang2pix(nside_sparse, np.radians(90.-dec), np.radians(ra), nest=True)

	#determine the unique pixels
	px_data_u = np.unique(px_data)

	#initially fill these pixels with zeroes so that HealSparse can update their values later
	all_maps.update_values_pix(px_data_u, np.zeros(len(px_data_u), dtype=dtype_maps))

	return all_maps, px_data, px_data_u


def pixelCountsFromCoords(ra, dec, nside_cover, nside_sparse):
	'''
	Given sets of coordinates (RA and Dec.), counts the number of objects in each pixel
	of a HealPIX map with the specified NSIDE.

	Parameters
	----------
	ra: array-like
		RA coordinates of each object.

	dec: array-like
		Dec. coordinates of each object.

	nside_cover: int
		Resolution of the wider HealSparse map where no data exist.

	nside_sparse: int
		Resolution of the regions of the HealSparse map in which data exist.

	Returns
	-------
	counts_map: HealSparseMap
		HealSparse map containing the number of sources in each pixel.
	'''
	#initialise a HealSparse integer map
	counts_map = hsp.HealSparseMap.make_empty(nside_cover, nside_sparse, np.int32)
	#convert the provided coordinates into pixel IDs
	px_data = hp.ang2pix(nside_sparse, np.radians(90.-dec), np.radians(ra), nest=True)
	#get the unique pixel IDs
	px_data_u = np.unique(px_data)

	#count the number of sources in each pixel (going from 0 to max(px_data))
	N = np.bincount(px_data).astype(np.int32)
	#identify all pixels containing at least one source
	good = N > 0
	#fill the map at these positions with the number of sources in the pixel
	counts_map[px_data_u] = N[good]

	return counts_map



def pixelMeanStd(quant, pix=None, remove_zeros=True):
	'''
	Calculate the mean and standard deviation of a given quantity at each pixel in a map.

	Parameters
	----------
	quant: array-like
		Values of the quantity being pixelised.

	pix: array-like or None
		Pixels corresponding to the coordinates at which the quantity is measured. 
		If None, RAs and Decs must be provided.

	remove_zeros: bool
		If True, only returns results for pixels containing data.
	
	Returns
	-------
	qmean: array-like
		Mean value of the quantity at each pixel.

	qstd: array-like
		Standard deviation of the quantity at each pixel.
	'''

	#count the number of values of quant associated with each pixel
	N = np.bincount(pix)
	#calculate the sum and the sum of the squares of the quantity in each pixel
	qsum = np.bincount(pix, weights=quant)
	qsqsum = np.bincount(pix, weights=quant**2.)

	if remove_zeros:
		#identify pixels in which at least one value of the quantity exists
		good = (N > 0)
		'''
		P = np.arange(pix.max()+1)
		pix_good = P[good]
		#delete P to conserve memory
		del P
		'''
		#calculate the mean and standard deviation of the quantity at each 'good' pixel
		qmean = qsum[good] / N[good]
		qmeansq = qmean ** 2.
		qvar = (qsqsum[good] / N[good]) - qmeansq
		#identify pixels for which the variance is negative and whose absolute value is <1e-10
		zero_var = np.isclose(np.zeros(len(qvar)), qvar, atol=1e-6*qmeansq)
		qvar[zero_var] = 0.
		qstd = np.sqrt(qvar)
	else:
		#suppress numpy DivideByZero warnings
		with np.errstate(divide='ignore', invalid='ignore'):
			#calculate the mean and standard deviation of the quantity at all pixels
			qmean = qsum / N
			qstd = np.sqrt((qsqsum / N) - qmean ** 2.)

	return qmean, qstd



def createMask(ra, dec, flags, nside_cover, nside_sparse):
	'''
	Creates a mask using a set of flags defined at the provided coordinates.

	Parameters
	----------
	ra: array-like
		RAs at which the flags are provided.

	dec: array-like
		Decs at which the flags are provided.

	flags: array-like or list of array-likes
		Boolean arrays containing the flag at each position. Pixels containing any
		True flags will be masked.

	'''

	#begin by counting sources in each pixel, since pixels with zero sources will be masked
	counts_map = pixelCountsFromCoords(ra, dec, nside_cover, nside_sparse)

	#initialise a mask with the same NSIDE parameters as the counts map
	mask = hsp.HealSparseMap.make_empty(nside_cover, nside_sparse, np.int32)
	#fill any occupied pixels with a 1
	mask[counts_map.valid_pixels] = 1

	#convert the provided RAs and Decs to pixel IDs in the new-resolution mask
	px_data = hp.ang2pix(nside_sparse, np.radians(90.-dec), np.radians(ra), nest=True)
	
	if type(flags) != list:
		flags = [flags]
	#cycle through the flags
	for flag in flags:
		#identify where the flag=True
		px_to_mask = px_data[flag]
		px_u = np.unique(px_to_mask)
		mask[px_u] = 0

	return mask


