#####################################################################################################
# Converts relevant catalogues and maps into the correct formats for TXPipe.
# Currently only set up to convert the systematics maps, as these are the only ones not created
# in TXPipe.
#####################################################################################################

import os
import config
import h5py
import healsparse as hsp
import numpy as np
from output_utils import colour_string
from map_utils import *
import h5py
import glob


### SETTINGS ###
cf = config.makeTXPipeInputs


###################
#### FUNCTIONS ####
###################

def make_tomography_cat(zbin_flags, fname):
	'''
	Makes an HDF file containing the information required by TXPipe's TXLensMaps stage.
	NOTE: This information includes data described as 'lens_weight', which is not relevant
	for this study. It is therefore assumed here that this quantity is 1 for all sources.

	Parameters
	----------
	zbin_flags: array-like
		Flags identifying the tomographic bin to which each galaxy belongs.

	fname: str
		Filename to be given to the output file. If the file already exists, will try
		to append to existing data in the file.
	'''

	#create an array to contain the counts in each bin
	counts = np.zeros(cf.nbins, dtype='i8')

	#cycle through the redshift bins
	for i in range(cf.nbins):
		zmask = zbin_flags == i
		counts[i] = zmask.sum()

	#create an array for containing the total counts in all bins
	counts_2d = np.array([counts.sum()])

	#create a column for 'lens_weight' and simply assign a value of 1 for every source
	lens_weight = np.ones(len(zbin_flags), dtype='f8')

	#compile the relevant data into a list and assign them names
	data = [zbin_flags, counts, counts_2d, lens_weight]
	names = ['bin', 'counts', 'counts_2d', 'lens_weight']

	#write to the output file
	with h5py.File(fname, 'w') as hf:
		g = hf.create_group('tomography')
		g.attrs['nbin'] = cf.nbins
		for d, n in zip(data, names):
			g.create_dataset(n, data=d, shape=(len(d),), dtype=d.dtype)
				

#######################################################
###############    START OF SCRIPT    #################
#######################################################


#cycle throught the fields
for fd in cf.get_global_fields():
	#path to directory containing the outputs for this field
	OUT = cf.PATH_OUT + fd + '/'
	#systematics maps directory
	PATH_SYST = OUT + 'systmaps/'
	#path in which the reformatted systematics maps will be placed
	PATH_TX = PATH_SYST + 'for_txpipe/'
	if not os.path.exists(PATH_TX):
		os.system(f'mkdir -p {PATH_TX}')

	#systematics maps
	maps = sorted(glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp'))
	#cycle through the maps
	for m in maps:
		#load the HealSparse map
		hsp_map = hsp.HealSparseMap.read(m)
		#retrieve the basename of the file
		bn = os.path.basename(m)
		#see if the input map contains a recarray
		if len(hsp_map.dtype) > 0:
			for n in hsp_map.dtype.names:
				fname = f'{PATH_TX}{bn[:-4]}_{n}.fits'
				healsparseToFITS(hsp_map[n], fname, nest=False)
		else:
			fname = f'{PATH_TX}{bn[:-4]}.fits'
			healsparseToFITS(hsp_map, fname, nest=False)

	#load the cleaned galaxy catalogue
	galcat = f'{OUT}clean_catalogue.hdf5'
	with h5py.File(galcat, 'r') as hf:
		#retrieve the tomographic bin labels
		zbin_flags = galcat[f'photometry/zbin'][:]
	#save relevant data to a tomographic catalogue
	tomo_out = PATH_TX + cf.cat_tomo
	make_tomography_cat(zbin_flags, tomo_out)