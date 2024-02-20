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
from matplotlib.lines import Line2D
import plot_utils as pu
import itertools


import faulthandler
faulthandler.enable()
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

'''
#retrieve bandpower edges from config
bpw_edges = np.array(cf.bpw_edges)
#lower and upper edges
bpw_edges_lo = bpw_edges[:-1]
bpw_edges_hi = bpw_edges[1:].copy()
#add 1 to the lower edges so that each multipole is only included in one bandpower
bpw_edges_lo += 1

#create pymaster NmtBin object using these bandpower objects
b = nmt.NmtBin.from_edges(bpw_edges_lo, bpw_edges_hi)
'''
#create pymaster NmtBin object using resolution of the maps
b = nmt.NmtBin.from_nside_linear(cf.nside_hi, 100)
ell_effs = b.get_effective_ells()
#use this to define the x-limits of the figures
xmin = ell_effs.min() / 1.5
xmax = ell_effs.max() * 1.2

#set up a pymaster Workspace object
w = nmt.NmtWorkspace()

#retrieve the number of redshift bins
nbins = len(cf.zbins) - 1
#use this to define the size of the power spectra figures
xsize = nbins * 4
ysize = nbins * 3.5

#maximum ell allowed by the resolution
ell_max = 3 * cf.nside_hi

#also get all possible pairings of bins
l = list(range(nbins))
pairings = [i for i in itertools.product(l,l) if tuple(reversed(i)) >= i]

#create a variable assignment that will later be occupied
cw = None


#cycle through the fields being analysed (TODO: later change to global fields)
for fd in cf.get_global_fields():
	#set up a figure for the power spectra from each redshift bin
	fig = plt.figure(figsize=(xsize, ysize))
	gs = fig.add_gridspec(ncols=nbins, nrows=nbins)

	#path to the directory containing the maps
	PATH_MAPS = f'{cf.PATH_OUT}{fd}/'
	#load the N_g and delta_g maps
	ngal_maps = hsp.HealSparseMap.read(PATH_MAPS+cf.ngal_maps)
	deltag_maps = hsp.HealSparseMap.read(PATH_MAPS+cf.deltag_maps)

	if len(cf.systs) > 0:
		#load the systematics maps and convert to full-sky realisations
		systmaps = load_maps([PATH_MAPS + s for s in cf.systs], normalise=True)
		#reshape the resultant list to have dimensions (nsyst, 1, npix)
		nsyst = len(systmaps)
		npix = len(systmaps[0])
		systmaps = np.array(systmaps).reshape([nsyst, 1, npix])
	else:
		systmaps = None

	#load the survey mask and convert to full-sky realisation
	mask, = load_maps([PATH_MAPS+cf.survey_mask])
	#identify pixels above the weight threshold
	above_thresh = mask > cf.weight_thresh
	#set all pixels below the weight threshold to 0
	mask[~above_thresh] = 0

	#full path to the output file
	outfile = f'{cf.PATH_OUT}{fd}/{cf.outfile}'
	#open the file, creating it if it doesn't exist
	with h5py.File(outfile, mode='a') as psfile:
		
		#cycle through all possible pairings of redshift bins
		for p in pairings:
			i,j = p

			#see if a group for the current pairing already exists
			p_str = str(p)
			if p_str in psfile.keys():
				#retrieve the relevant qualities for creating the plot
				cl_bias_decoupled = psfile[f'{p_str}/cl_bias_decoupled'][...]
				cl_decoupled_debiased = psfile[f'{p_str}/cl_decoupled_debiased'][...]
				err_cell = psfile[f'{p_str}/err_cell'][...]
				N_ell_decoupled = psfile[f'{p_str}/N_ell_decoupled'][...]

			else:

				#load the density maps for the current pairing and create NmtFields
				dg_i = hspToFullSky(deltag_maps[f'delta_{i}'])
				f_i = nmt.NmtField(mask, [dg_i], templates=systmaps)
				if j == i:
					f_j = nmt.NmtField(mask, [dg_i], templates=systmaps)
				else:
					dg_j = hspToFullSky(deltag_maps[f'delta_{j}'])
					f_j = nmt.NmtField(mask, [dg_j], templates=systmaps)


				cl_coupled = nmt.workspaces.compute_coupled_cell(f_i, f_j)
				#use these along with the mask to get a guess of the true C_ell
				cl_guess = cl_coupled / np.mean(mask * mask)


				#compute the mode coupling matrix
				w.compute_coupling_matrix(f_i, f_j, b)
				#compute the decoupled C_ell (w/o deprojection)
				cl_decoupled = w.decouple_cell(cl_coupled)

				#only calculate bias-related quantities if templates have been provided
				if systmaps is not None:
					#compute the deprojection bias
					cl_bias = nmt.deprojection_bias(f_i, f_j, cl_guess)
					#compute the decoupled C_ell (w/ deprojection)
					cl_decoupled_debiased = w.decouple_cell(cl_coupled, cl_bias=cl_bias)
					#decouple the bias C_ells as well
					cl_bias_decoupled = w.decouple_cell(cl_bias)
				else:
					cl_bias = cl_bias_decoupled = None
					cl_decoupled_debiased = cl_decoupled[...]


				########################
				# NOISE POWER SPECTRUM #
				########################

				#Only calculate for autocorrelations
				if i == j:
					#retrieve the Ngal map for this redshift bin
					ng_full = hspToFullSky(ngal_maps[f'ngal_{i}'])
					#calculate the mean number of galaxies per pixel (above weight threshold)
					mu_N = np.sum(ng_full[above_thresh]) / np.sum(mask[above_thresh])
					#calculate the mean value of the mask
					mu_w = np.mean(mask)
					#get pixel area in units of steradians
					Apix = hp.nside2pixarea(cf.nside_hi)

					#calculate the noise power spectrum
					N_ell_coupled = np.full(ell_max, Apix * mu_w / mu_N).reshape((1,ell_max))
					#decouple
					N_ell_decoupled = w.decouple_cell(N_ell_coupled)

				if cw is None:
					#covariance matrix calculation: create covariance workspace
					cw = nmt.NmtCovarianceWorkspace()
					#compute coupling coefficients
					cw.compute_coupling_coefficients(f_i, f_j)

				#extract (gaussian) covariance matrix
				n_ell = len(cl_decoupled[0])
				covar = nmt.gaussian_covariance(cw, 
												0, 0, 0, 0,			#spin of each field
												[cl_guess[0]],	
												[cl_guess[0]],
												[cl_guess[0]],
												[cl_guess[0]],
												w)
				#errorbars for each bandpower
				err_cell = np.diag(covar) ** 0.5


				#populate the output file with the results
				gp = psfile.require_group(p_str)
				_ = gp.create_dataset('ell_effs', data=ell_effs)
				_ = gp.create_dataset('cl_coupled', data=cl_coupled)
				_ = gp.create_dataset('cl_decoupled', data=cl_decoupled)
				_ = gp.create_dataset('cl_decoupled_debiased', data=cl_decoupled_debiased)
				_ = gp.create_dataset('cl_guess', data=cl_guess)
				_ = gp.create_dataset('cl_bias', data=cl_bias)
				_ = gp.create_dataset('cl_bias_decoupled', data=cl_bias_decoupled)
				_ = gp.create_dataset('N_ell_coupled', data=N_ell_coupled)
				_ = gp.create_dataset('N_ell_decoupled', data=N_ell_decoupled)
				_ = gp.create_dataset('covar', data=covar)
				_ = gp.create_dataset('err_cell', data=err_cell)



			####################
			# PLOTTING RESULTS #
			####################

			#add subplot to gridspec
			ax = fig.add_subplot(gs[j,i])
			#only label axes if on outer edge of figure
			if j == (nbins-1):
				ax.set_xlabel(r'$\ell$')
			if i == 0:
				ax.set_ylabel(r'$C_{\ell}$')
			#set loglog scale
			ax.set_xscale('log')
			ax.set_yscale('log')

			if cl_bias_decoupled is not None:
				#plot the deporojection bias
				bias_plot, *_ = ax.plot(b.get_effective_ells(), cl_bias_decoupled[0], c=pu.magenta)
				ax.plot(b.get_effective_ells(), -cl_bias_decoupled[0], ls='--', c=pu.magenta)
				#at this point retrieve the x limits
				#xmin, xmax = ax.get_xlim()

			#plot the debiased power spectrum, using open symbols for abs(negative) values
			mask_pve = cl_decoupled_debiased[0] > 0
			mask_nve = cl_decoupled_debiased[0] <= 0
			Y_pve = cl_decoupled_debiased[0][mask_pve]
			Y_nve = cl_decoupled_debiased[0][mask_nve]
			cell_plot = ax.errorbar(b.get_effective_ells()[mask_pve], Y_pve, yerr=err_cell[mask_pve], marker='o', c=pu.dark_blue, linestyle='none')
			ax.errorbar(b.get_effective_ells()[mask_nve], -Y_nve, yerr=err_cell[mask_nve], marker='o', markeredgecolor=pu.dark_blue, markerfacecolor='none', linestyle='none')
			

			#plot the shot noise if autocorrelation
			if i == j:
				noise_plot, *_ = ax.plot(b.get_effective_ells(), N_ell_decoupled[0], c=pu.teal)

			#reset the axis limits
			ax.set_xlim(xmin, xmax)

			#add text to the top-right corner to indicate which bins have been compared
			ax.text(0.95, 0.95, f'({i},{j})', transform=ax.transAxes, ha='right', va='top', fontsize=20.)

	#create a legend		
	handles = [
		cell_plot,
		noise_plot
		]
	labels = [
		'Signal',
		'Noise'
		]
	if cl_bias_decoupled is not None:
		handles.insert(1, bias_plot)
		labels.insert(1, 'Deprojection bias')
	

	fig.legend(handles=handles, labels=labels, loc='upper right', fontsize=28)

	plt.tight_layout()
	plt.savefig(f'{cf.PATH_PLOTS}{fd}_power_spectra.png', dpi=300)
		



