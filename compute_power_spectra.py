#####################################################################################################
# - Uses NaMaster to compute power spectra from the galaxy delta_g maps, deprojecting any 
#   systematics templates in the process.
# - TODO: ensure that script doesn't use previous output if deprojection is to occur and didn't on previous
#	run.
#####################################################################################################

import config
import healpy as hp
import healsparse as hsp
import numpy as np
from map_utils import *
import h5py
import pymaster as nmt
from matplotlib import pyplot as plt
import plot_utils as pu
import itertools
from output_utils import colour_string
import os
import glob

import faulthandler
faulthandler.enable()
### SETTINGS ###
cf = config.computePowerSpectra
plt.style.use(pu.styledict)


###################
#### FUNCTIONS ####
###################


def hspToFullSky(hsp_map):
	'''
	Converts a HealSparse map into a full-sky realisation that NaMaster can use.

	Parameters
	----------
	hsp_map: HealSparseMap
		The map being converted.

	is_systmap: bool
		If True, subtracts the mean form the value of all pixels (necessary for systematics maps).

	Returns
	-------
	fs_map: array
		Full-sky realisation containing values at every pixel (hp.UNSEEN in unoccupied pixels).
	'''
	
	#initialise an empty full-sky map (NOTE: can take a lot of memory for high nside_sparse)
	fs_map = hsp_map.generate_healpix_map(nest=False)
	vpix = fs_map != hp.UNSEEN
	fs_map[~vpix] = 0.

	return fs_map

	


def load_maps(map_paths):
	'''
	Loads any number of provided HealSparse maps and returns their pixel values in healPIX 
	RING ordering.

	TODO: add functionality for reading maps in hdf5 format?

	Parameters
	----------
	map_paths: list
		List of paths to the maps being read.

	is_systmap: bool
		If True, subtracts the mean form the value of all pixels (necessary for systematics maps).

	Returns
	-------
	maps: list
		List containing data from each map.
	'''
	maps = []

	#cycle through the list of paths
	for mp in map_paths:

		#use healSparse to read the maps
		map_data = hsp.HealSparseMap.read(mp)
		#see if multiple maps are stored in this file
		nmaps = len(map_data.dtype)
		if nmaps > 0:
			#append the data from each map to the list
			maps.extend([hspToFullSky(map_data[n]) for n in map_data.dtype.names])
		else:
			#append the data to the list
			maps.append(hspToFullSky(map_data))

	return maps





#######################################################
###############    START OF SCRIPT    #################
#######################################################

#maximum ell allowed by the resolution
ell_max = 3 * cf.nside_hi - 1
#get pixel area in units of steradians
Apix = hp.nside2pixarea(cf.nside_hi)


if cf.use_N19_bps:
	#retrieve bandpower edges from config
	bpw_edges = np.array(cf.bpw_edges).astype(int)
	#only include bandpowers < 3 * NSIDE
	bpw_edges = bpw_edges[bpw_edges <= ell_max]
else:
	if cf.log_spacing:
		bpw_edges = np.geomspace(cf.ell_min, ell_max, cf.nbpws).astype(int)
	else:
		bpw_edges = np.linspace(cf.ell_min, ell_max, cf.nbpws).astype(int)
#create pymaster NmtBin object using these bandpower objects
b = nmt.NmtBin.from_edges(bpw_edges[:-1], bpw_edges[1:])


#get the effective ells
ell_effs = b.get_effective_ells()
#use this to define the x-limits of the figures
xmin = ell_effs.min() / 1.5
xmax = ell_effs.max() * 1.2



#retrieve the number of redshift bins
nbins = len(cf.zbins) - 1
#use this to define the size of the power spectra figures
xsize = nbins * 4
ysize = nbins * 3.5


#also get all possible pairings of bins
l = list(range(nbins))
pairings = [i for i in itertools.product(l,l) if tuple(reversed(i)) >= i]


#cycle through the fields being analysed (TODO: later change to global fields)
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))

	#set up a figure for the power spectra from each redshift bin
	fig = plt.figure(figsize=(xsize, ysize))
	gs = fig.add_gridspec(ncols=nbins, nrows=nbins)

	#path to the directory containing the maps
	PATH_MAPS = f'{cf.PATH_OUT}{fd}/'

	#set up a pymaster Workspace object
	w = nmt.NmtWorkspace()
	#create a variable assignment that will later be occupied by a CovarianceWorkspace
	cw = nmt.NmtCovarianceWorkspace()

	PATH_CACHE = PATH_MAPS + 'cache/'
	#see if directory for cached workspaces exists; make it if not
	if not os.path.exists(PATH_CACHE):
		os.system(f'mkdir -p {PATH_CACHE}')
	
	#path to directory containing systematics maps
	PATH_SYST = f'{PATH_MAPS}systmaps/'
	#check for 'All' in systmaps and convert this to a list of all systematics maps
	if 'all' in map(str.lower, cf.systs):
		cf.systs = [os.path.basename(m) for m in (glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp'))]
	#if given a max number of systematics to deproject, slice the list accordingly
	if cf.Nsyst_max is not None:
		cf.systs = cf.systs[:cf.Nsyst_max]

	#file containing list of systematics maps deprojected in the previous run
	deproj_file = PATH_CACHE + cf.deproj_file
	if os.path.exists(deproj_file):
		with open(deproj_file, 'r+') as df:
			#see which (if any) systematics have been deprojected previously
			deproj_done = df.read().split('\n')
			#see if this is the same as the list specified in the config file (accounting for different ordering)
			if sorted(deproj_done) == sorted(cf.systs):
				calc = False
				print('Same systematics maps provided')
			else:
				calc = True
				print('Different systematics maps provided')
				#write the list of provided systematics to the file
				df.seek(0)
				df.truncate()
				df.write('\n'.join(cf.systs))
	else:
		if len(cf.systs) == 0:
			calc = False
			print('No systematics provided')
		else:
			calc = True
			with open(deproj_file, 'w') as df:
				df.write('\n'.join(cf.systs))
		

	#see if workspaces have already been created from a previous run
	wsp_path = PATH_CACHE + cf.wsp_file
	covwsp_path = PATH_CACHE + cf.covwsp_file
	if os.path.exists(wsp_path) and not calc:
		w.read_from(wsp_path)
	else:
		calc = True
	if os.path.exists(covwsp_path) and not calc:
		cw.read_from(covwsp_path)
	else:
		calc = True
		

	#load the delta_g maps
	deltag_maps = load_maps([PATH_MAPS+cf.deltag_maps])

	#load the survey mask and convert to full-sky realisation
	mask, = load_maps([PATH_MAPS+cf.survey_mask])
	#identify pixels above the weight threshold
	above_thresh = mask > cf.weight_thresh
	#set all pixels below the weight threshold to 0
	mask[~above_thresh] = 0

	#calculate anything to do with the mask so that it can also be cleared from memory
	sum_w_above_thresh = np.sum(mask[above_thresh])
	mu_w = np.mean(mask)
	mu_w2 = np.mean(mask * mask)


	print('Loading systematics maps...')
	if len(cf.systs) > 0:
		#load the systematics maps and convert to full-sky realisations
		systmaps = load_maps([PATH_SYST + s for s in cf.systs])
		#calculate the weighted mean of each systematics map and subtract it
		for sm in systmaps:
			mu_s = np.sum(sm[above_thresh] * mask[above_thresh]) / sum_w_above_thresh
			sm[above_thresh] -= mu_s
			print('Syst map mean: ', mu_s)
		#reshape the resultant list to have dimensions (nsyst, 1, npix)
		nsyst = len(systmaps)
		npix = len(systmaps[0])
		systmaps = np.array(systmaps).reshape([nsyst, 1, npix])
		deproj = True
		print('templates: ', np.mean(systmaps))
	else:
		systmaps = None
		deproj = False
	print('Done!')


	print('Creating NmtFields...')
	density_fields = [nmt.NmtField(mask, [d], templates=systmaps, n_iter=0) for d in deltag_maps]
	print('Done!')

	#clear some memory
	del systmaps	
	del mask
	del deltag_maps


	#load the N_g maps and calculate the mean weighted by the mask
	mu_N_all = [nmap[above_thresh].sum() / sum_w_above_thresh for nmap in load_maps([PATH_MAPS+cf.ngal_maps])]


	#full path to the output file
	outfile = f'{cf.PATH_OUT}{fd}/{cf.outfile}'
	#open the file, creating it if it doesn't exist
	with h5py.File(outfile, mode='w') as psfile:
		
		#cycle through all possible pairings of redshift bins
		for ip,p in enumerate(pairings):
			i,j = p

			#see if a group for the current pairing already exists
			p_str = str(p)
			print(colour_string(p_str, 'green'))

			f_i = density_fields[i]
			f_j = density_fields[j]
			cl_coupled = nmt.compute_coupled_cell(f_i, f_j)

			#use these along with the mask to get a guess of the true C_ell
			cl_guess = cl_coupled / mu_w2

			if ip == 0 and calc:
				#compute the mode coupling matrix (only need to compute once since same mask used for everything)
				print('Computing mode coupling matrix...')
				w.compute_coupling_matrix(f_i, f_j, b, n_iter=1)
				print('Done!')

				print('Calculating coupling coefficients...')
				#compute coupling coefficients
				cw.compute_coupling_coefficients(f_i, f_j)
				print('Done!')
			else:
				print('Using coupling matrix and coefficients from cache.')

			#only calculate bias-related quantities if templates have been provided
			if deproj:
				if calc:
					print('Calculating deprojection bias...')
					#compute the deprojection bias
					cl_bias = nmt.deprojection_bias(f_i, f_j, cl_guess)
					print('Done!')
				else:
					print('Combination of systematics matches previous run; using cached results.')
			else:
				print('No systematics maps provided; skipping deprojection bias calculation.')
				cl_bias = np.zeros_like(cl_guess)

			#multiplicative correction to delta_g of (1 / (1-Fs)) due to stars results in factor of (1 / (1 - Fs))^2 correction to Cl
			if cf.correct_for_stars:
				mult = (1 / (1 - cf.Fs_fiducial)) ** 2.
				cl_coupled *= mult
				cl_guess *= mult

			#compute the decoupled C_ell (w/o deprojection)
			cl_decoupled = w.decouple_cell(cl_coupled)
			#compute the decoupled C_ell (w/ deprojection)
			cl_decoupled_debiased = w.decouple_cell(cl_coupled, cl_bias=cl_bias)
			#decouple the bias C_ells as well
			cl_bias_decoupled = w.decouple_cell(cl_bias)


			########################
			# NOISE POWER SPECTRUM #
			########################

			#Only calculate for autocorrelations
			if i == j:
				mu_N = mu_N_all[i]
				#calculate the noise power spectrum
				N_ell_coupled = np.full(ell_max, Apix * mu_w / mu_N).reshape((1,ell_max))
				#decouple
				N_ell_decoupled = w.decouple_cell(N_ell_coupled)

			
			#######################
			# GAUSSIAN COVARIANCE #
			#######################

			#extract (gaussian) covariance matrix
			print('Calculating covariance matrix...')
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
			print('Done!')

			##################
			# SAVING RESULTS #
			##################

			#populate the output file with the results
			gp = psfile.require_group(p_str)
			_ = gp.create_dataset('ell_effs', data=ell_effs)
			_ = gp.create_dataset('cl_coupled', data=cl_coupled)
			_ = gp.create_dataset('cl_decoupled', data=cl_decoupled)
			_ = gp.create_dataset('cl_guess', data=cl_guess)
			_ = gp.create_dataset('N_ell_coupled', data=N_ell_coupled)
			_ = gp.create_dataset('N_ell_decoupled', data=N_ell_decoupled)
			_ = gp.create_dataset('covar', data=covar)
			_ = gp.create_dataset('err_cell', data=err_cell)
			_ = gp.create_dataset('cl_bias', data=cl_bias)
			_ = gp.create_dataset('cl_bias_decoupled', data=cl_bias_decoupled)
			_ = gp.create_dataset('cl_decoupled_debiased', data=cl_decoupled_debiased)



			####################
			# PLOTTING RESULTS #
			####################

			if cf.normalise:
				ylabel = r'$C_{\ell}\frac{\ell(\ell+1)}{2\pi}$'
				mfactor = ell_effs * (ell_effs + 1) / (2. * np.pi)
			else:
				ylabel = r'$C_{\ell}$'
				mfactor = np.ones_like(ell_effs)


			#add subplot to gridspec
			ax = fig.add_subplot(gs[j,i])
			#only label axes if on outer edge of figure
			if j == (nbins-1):
				ax.set_xlabel(r'$\ell$')
			if i == 0:
				ax.set_ylabel(ylabel)
			#set loglog scale
			ax.set_xscale('log')
			ax.set_yscale('log')

			if cl_bias_decoupled.any():
				#plot the deporojection bias
				bias_plot, *_ = ax.plot(ell_effs, cl_bias_decoupled[0]*mfactor, c=pu.magenta)
				ax.plot(ell_effs, -cl_bias_decoupled[0]*mfactor, ls='--', c=pu.magenta)
				#at this point retrieve the x limits
				#xmin, xmax = ax.get_xlim()

			#plot the debiased power spectrum, using open symbols for abs(negative) values
			mask_pve = cl_decoupled_debiased[0] > 0
			mask_nve = cl_decoupled_debiased[0] <= 0
			#plot the shot noise if autocorrelation
			if i == j:
				Y_pve = (cl_decoupled_debiased[0][mask_pve] - N_ell_decoupled[0][mask_pve]) * mfactor[mask_pve]
				Y_nve = (cl_decoupled_debiased[0][mask_nve] - N_ell_decoupled[0][mask_nve]) * mfactor[mask_nve]
				noise_plot, *_ = ax.plot(ell_effs, N_ell_decoupled[0]*mfactor, c=pu.teal)
			else:
				Y_pve = cl_decoupled_debiased[0][mask_pve] * mfactor[mask_pve]
				Y_nve = cl_decoupled_debiased[0][mask_nve] * mfactor[mask_nve]
			cell_plot = ax.errorbar(ell_effs[mask_pve], Y_pve, yerr=err_cell[mask_pve], marker='o', c=pu.dark_blue, linestyle='none')
			ax.errorbar(ell_effs[mask_nve], -Y_nve, yerr=err_cell[mask_nve], marker='o', markeredgecolor=pu.dark_blue, markerfacecolor='none', linestyle='none')

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
	if cl_bias_decoupled.any():
		handles.insert(1, bias_plot)
		labels.insert(1, 'Deprojection bias')
		#figure name also depends on whether deprojection has occurred
		figname = f'{cf.PATH_PLOTS}{fd}_power_spectra_{cf.nside_hi}.png'
	else:
		figname = f'{cf.PATH_PLOTS}{fd}_power_spectra_raw_{cf.nside_hi}.png'
	

	fig.legend(handles=handles, labels=labels, loc='upper right', fontsize=28)

	plt.tight_layout()
	plt.savefig(figname, dpi=300)


	######################
	# CACHING WORKSPACES #
	######################

	#write the workspaces to the cache directory
	w.write_to(f'{PATH_CACHE}{cf.wsp_file}')
	cw.write_to(f'{PATH_CACHE}{cf.covwsp_file}')
		



