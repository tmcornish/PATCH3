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
		primary = labels[0]
	all_maps = hsp.HealSparseMap.make_empty(nside_cover, nside_sparse, dtype=dtype_maps, primary=primary)

	#identify pixels in the HealSparse map that contain sources from the catalogue
	px_data = hp.ang2pix(nside_sparse, np.radians(90.-dec), np.radians(ra), nest=True)

	#determine the unique pixels
	px_data_u = np.unique(px_data)

	#initially fill these pixels with zeroes so that HealSparse can update their values later
	all_maps.update_values_pix(px_data_u, np.zeros(len(px_data_u), dtype=dtype_maps))

	return all_maps, px_data, px_data_u


def pixelCountsFromCoords(ra, dec, nside_cover, nside_sparse, return_pix_and_vals=False):
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

	return_pix_and_vals: bool
		If True, also returns the IDs and corresponding values for occupied pixels.

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

	if return_pix_and_vals:
		return counts_map, px_data_u, N[good]
	else:
		return counts_map



def pixelMeanStd(quant, pix, remove_zeros=True):
	'''
	Calculate the mean and standard deviation of a given quantity at each pixel in a map.

	Parameters
	----------
	quant: array-like
		Values of the quantity being pixelised.

	pix: array-like or None
		Pixels corresponding to the coordinates at which the quantity is measured.

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
		#calculate the mean and standard deviation of the quantity at each 'good' pixel
		qmean = qsum[good] / N[good]
		qmeansq = qmean ** 2.
		qvar = (qsqsum[good] - (2 * qmean * qsum[good])) / N[good] + qmeansq
		#identify pixels for which the variance is < a millionth of the squared mean (arises from rounding errors)
		zero_var = np.isclose(np.zeros(len(qvar)), qvar, atol=1e-6*qmeansq)
		qvar[zero_var] = 0.
		qstd = np.sqrt(qvar)
	else:
		#suppress numpy DivideByZero warnings
		with np.errstate(divide='ignore', invalid='ignore'):
			#calculate the mean and standard deviation of the quantity at all pixels
			qmean = qsum / N
			qstd = np.sqrt((qsqsum - (2 * qmean * qsum)) / N + (qmean ** 2.))

	return qmean, qstd


def createMeanStdMap(ra, dec, quant, nside_cover, nside_sparse):
	'''
	Creates maps containing the mean and standard deviation of a given quantity in each pixel.

	Parameters
	----------
	ra: array-like
		RAs at which the flags are provided.

	dec: array-like
		Decs at which the flags are provided.

	quant: array-like or list of array-likes
		Values of the desired quantity at each position defined by the provided RAs and Decs.

	nside_cover: int
		NSIDE parameter definining the low-resolution regions of the map (where no data exist).

	nside_sparse: int
		NSIDE parameter definining the high-resolution regions of the map (where data exist).

	'''

	#set up two maps with the desired resolutions (one for mean and one for std)
	mean_map = hsp.HealSparseMap.make_empty(nside_cover, nside_sparse, np.float64)
	std_map = hsp.HealSparseMap.make_empty(nside_cover, nside_sparse, np.float64)

	#convert the provided coordinates into pixel coordinates in the high-resolution map
	px_data = hp.ang2pix(nside_sparse, np.radians(90.-dec), np.radians(ra), nest=True)
	px_data_u = np.unique(px_data)

	#calculate the mean and std of the quantity at each pixel and populate the maps
	mean_map[px_data_u], std_map[px_data_u] = pixelMeanStd(quant, px_data, remove_zeros=True)
	
	return mean_map, std_map


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

	nside_cover: int
		NSIDE parameter definining the low-resolution regions of the map (where no data exist).

	nside_sparse: int
		NSIDE parameter definining the high-resolution regions of the map (where data exist).

	Returns
	-------
	mask: HealSparseMap
		Map containing 0s at masked positions and 1s elsewhere.
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


def maskAreaSkyCoverage(mask, thresh=0.):
	'''
	Takes a HealSparseMap mask and calculates the unmasked area and the fractional sky coverage. 
	Assumes that pixels with smaller values are `more heavily masked'. 

	Parameters
	----------
	mask: HealSparseMap
		The mask for which the unmasked area and fractional sky coverage will be calculated.

	thresh: float
		Threshold below which pixels will be classed as masked.

	Returns
	-------
	A: float
		Unmasked area (in square degrees).

	f_sky: float
		Fractional sky coverage.
	'''
	
	#calculate the area of one pixel in square degrees
	A_pix = hp.nside2pixarea(mask.nside_sparse, degrees=True)
	#determine the number of pixels above the mask threshold
	vpix = mask.valid_pixels
	N_hits = mask[vpix][mask[vpix] > thresh].sum()
	#calculate the total unmasked area and fractional sky coverage
	A = N_hits * A_pix
	f_sky = A / (4. * np.pi * (180. / np.pi) ** 2.)

	return A, f_sky


def healsparseToHDF(hsp_map, fname, pix_scheme='ring', group='', metadata=None):
	'''
	Takes a HealSparse map and reformats it to an HDF5 dataset containing the pixel IDs
	and values.

	Parameters
	----------
	hsp_map: HealSparseMap
		Map to be converted from HealSparseMap to HDF5 format.

	fname: str
		Filename for the output HDF5 file.

	pix_scheme: str
		The desired pixelisation scheme of the map. Either 'ring' or 'nest'. Input map
		is assumed to be 'nest' as HealSparse maps have this scheme by definition.

	group: str
		Group within which the relevant data are expected to reside.

	metadata: dict or None
		Dictionary of additional metadata to be stored for the specified group.
	'''

	import h5py

	#get the valid pixels and correpsonding values from the map
	vpix = hsp_map.valid_pixels
	values = hsp_map[vpix]
	#if output is to be in ring scheme, convert pixel IDs to ring equivalent
	if pix_scheme.lower() == 'ring':
		vpix = hp.nest2ring(hsp_map.nside_sparse, vpix)
	elif pix_scheme.lower() == 'nest':
		pass
	else:
		raise ValueError('pix_scheme must be either "nest" or "ring".')

	with h5py.File(fname, 'w') as hf:
		#create the relevant group if it doesn't already exist
		if len(group) != 0:
			g = hf.create_group(group)
		else:
			g = hf
		#add some metadata to describe the pixelisation
		g.attrs['pixelization'] = 'healpix'
		g.attrs['nside'] = hsp_map.nside_sparse
		#add any additional metadata
		if metadata is not None:
			for md in metadata:
				#g.attrs = dict(g.attrs, **metadata)
				g.attrs[md] = metadata[md]
		_ = hf.create_dataset(f'{group}/pixel', data=vpix, dtype=vpix.dtype)
		_ = hf.create_dataset(f'{group}/value', data=values, dtype=values.dtype)



def healsparseToFITS(hsp_map, fname, nest=False):
	'''
	Takes a HealSparse map and reformats it to healpy FITS map.

	Parameters
	----------
	hsp_map: HealSparseMap
		Map to be converted from HealSparseMap to FITS format.

	fname: str
		Filename for the output HDF5 file.

	nest: bool
		Whether the output map should be in 'nest' format. Input map is assumed
		to be 'nest' as HealSparse maps have this scheme by definition.
	'''

	#convert the map to healpy format
	hp_map = hsp_map.generate_healpix_map(nest=nest)
	#write to file
	hp.write_map(fname, hp_map, nest=nest, column_names=['VALUE'], overwrite=True)



class MaskData:
	'''
	Convenience class for containing important information about the survey mask. Upon
	initialisation only requires the name of the HealSparse file.
	'''

	def __init__(self, hsp_file, mask_thresh=0.):
		#load the HealSparse map and convert it to full sky
		self.mask = hsp.HealSparseMap.read(hsp_file).generate_healpix_map(nest=False)
		#apply the mask threshold
		self.mask[self.mask <= mask_thresh] = 0.
		#get the IDs of all pixels above the threshold
		self.vpix = np.argwhere(self.mask > 0.).flatten()

		#compute the sum, mean, and mean squared of the mask
		self.sum = np.sum(self.mask)
		self.mean = np.mean(self.mask)
		self.meansq = self.mean ** 2.

