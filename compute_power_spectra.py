#####################################################################################################
# - Uses NaMaster to compute power spectra from the galaxy delta_g maps, deprojecting any 
#   systematics templates in the process.
#####################################################################################################

import os
import config
import healpy as hp
import healsparse as hsp
import numpy as np
from output_utils import colour_string
from map_utils import *
import h5py
import pymaster as nmt
from matplotlib import pyplot as plt
import plot_utils as pu



### SETTINGS ###
cf = config.computePowerSpectra
plt.style.use(pu.styledict)


###################
#### FUNCTIONS ####
###################


def hspToFullSky(hsp_map, normalise=False):
	'''
	Converts a HealSparse map into a full-sky realisation that NaMaster can use.

	Parameters
	----------
	hsp_map: HealSparseMap
		The map being converted.

	normalise: bool
		If True, normalises the data by performing the following operation to each valid pixel:
		value --> (value / <value>) - 1

	Returns
	-------
	fs_map: array
		Full-sky realisation containing values at every pixel (hp.UNSEEN in unoccupied pixels).
	'''

	#initialise an empty full-sky map (NOTE: can take a lot of memory for high nside_sparse)
	npix = hp.nside2npix(hsp_map.nside_sparse)
	fs_map = np.full(npix, 0, dtype=hsp_map.dtype)
	#find the valid pixels in the map and convert from NEST to RING ordering
	vpix = hsp_map.valid_pixels
	vpix_ring = hp.nest2ring(cf.nside_hi, vpix)

	#fill the relevant pixels in the full-sky map
	fs_map[vpix_ring] = hsp_map[vpix]
	#if told to normalise, divide by the mean and subtract 1
	if normalise:
		mu = np.mean(hsp_map[vpix])
		fs_map /= mu
		fs_map -= 1

	return fs_map

	


def load_maps(map_paths, normalise=False):
	'''
	Loads any number of provided HealSparse maps and returns their pixel values in healPIX 
	RING ordering.

	TODO: add functionality for reading maps in hdf5 format?

	Parameters
	----------
	map_paths: list
		List of paths to the maps being read.

	normalise: bool
		If True, normalises the data by performing the following operation to each valid pixel:
		value --> (value / <value>) - 1

	Returns
	-------
	maps: list
		List containing data from each map.
	'''
	#determine how many pixels are in a full-sky map at the resolution of the maps
	npix = hp.nside2npix(cf.nside_hi)

	maps = []

	#cycle through the list of paths
	for mp in map_paths:

		#use healSparse to read the maps
		map_data = hsp.HealSparseMap.read(mp)
		#see if multiple maps are stored in this file
		nmaps = len(map_data.dtype)
		if nmaps > 0:
			#append the data from each map to the list
			maps.extend([hspToFullSky(map_data[n], normalise) for n in map_data.dtype.names])
		else:
			#append the data to the list
			maps.append(hspToFullSky(map_data, normalise))

	return maps





#######################################################
###############    START OF SCRIPT    #################
#######################################################

#retrieve bandpower edges from config
bpw_edges = np.array(cf.bpw_edges)
#lower and upper edges
bpw_edges_lo = bpw_edges[:-1]
bpw_edges_hi = bpw_edges[1:].copy()
#add 1 to the lower edges so that each multipole is only included in one bandpower
bpw_edges_lo += 1

#create pymaster NmtBin object using these bandpower objects
b = nmt.NmtBin.from_edges(bpw_edges_lo, bpw_edges_hi)
#set up a pymaster Workspace object
w = nmt.NmtWorkspace()

#cycle through the fields being analysed (TODO: later change to global fields)
for fd in cf.fields:
	#path to the directory containing the maps
	PATH_MAPS = f'{cf.PATH_OUT}{fd}/'
	#load the delta_g maps
	deltag_maps = hsp.HealSparseMap.read(PATH_MAPS+cf.deltag_maps)

	#load the systematics maps and convert to full-sky realisations
	systmaps = load_maps([PATH_MAPS + s for s in cf.systs], normalise=True)
	#reshape the resultant list to have dimensions (nsyst, 1, npix)
	nsyst = len(systmaps)
	npix = len(systmaps[0])
	systmaps = np.array(systmaps).reshape([nsyst, 1, npix])
	#load the survey mask and convert to full-sky realisation
	mask, = load_maps([PATH_MAPS+cf.survey_mask])

	#cycle through the redshift bins (TODO: edit to allow for all pairings of z bins)
	for i in [0]:#range(len(cf.zbins)-1):
		dg_full = hspToFullSky(deltag_maps[f'delta_{i}'])
		#create an NmtField object
		f = nmt.NmtField(mask, [dg_full], templates=systmaps)
		#compute the mode coupling matrix
		w.compute_coupling_matrix(f, f, b)
		#compute full estimate of the power spectrum
		#TODO: bias deprojection, noise power spectrum calculation
		c = nmt.workspaces.compute_full_master(f, f, b, workspace=w)

		f, ax = plt.subplots()
		ax.set_ylabel(r'$\ell$')
		ax.set_ylabel(r'$C_{\ell}$')

		ax.plot(b.get_effective_ells(), c[0], label=f'Bin {i}')
		ax.legend()

		ax.set_xscale('log')
		ax.set_yscale('log')
		plt.tight_layout()
		plt.savefig(f'{PATH_MAPS}power_spectra_test.png', dpi=300)



