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
from map_utils import initialiseRecMap, pixelMeanStd, pixelCountsFromCoords, createMask
import h5py


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




#######################################################
###############    START OF SCRIPT    #################
#######################################################

#cycle through each of the fields
for fd in cf.fields:
	#output directory for this field
	OUT = cf.PATH_OUT + fd
	#load the basic and fully cleaned galaxy catalogues for this field
	cat_basic = h5py.File(f'{OUT}/{cf.cat_basic}', 'r')
	cat_main = h5py.File(f'{OUT}/{cf.cat_main}', 'r')

	#make the dust maps in each band and store in a single recarray
	dust_maps = makeDustMap(cat_basic, group='photometry')
	#write to a file
	dust_maps.write(f'{OUT}/{cf.dustmaps}', clobber=True)

	#make the bright object mask
	bo_mask = makeBOMask(cat_basic, group='photometry')
	#write to a file
	bo_mask.write(f'{OUT}/{cf.bo_mask}', clobber=True)
