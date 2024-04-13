#####################################################################################################
# - Uses NaMaster to compute power spectra from the galaxy delta_g maps, deprojecting any 
#   systematics templates in the process.
# - TODO: skip all calculations if systematics identical to previous run.
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


def load_map(map_path, apply_mask=False, is_systmap=False, mask=None):
	'''
	Loads an individual HealSparse map and returns their pixel values in healPIX 
	RING ordering. If told to, will also multiply the map by the mask and/or 
	calculate the mean of the map and subtract it from all pixels. Both of these
	operations require the mask (in full HealPIX format) as input. 


	Parameters
	----------
	map_path: str
		Path to the map being read.

	apply_mask: bool
		If True, will perform element-wise multiplication by the mask (given as separate input).

	is_systmap: bool
		If True, subtracts the mean from the value of all pixels (necessary for systematics maps).
	
	mask: MaskData
		Object containing the mask and any potentially relevant pre-computed values.

	Returns
	-------
	fs_map: np.array
		Full-sky map data (RING ordering).
	'''

	#initialise an empty full-sky map (NOTE: can take a lot of memory for high nside_sparse)
	fs_map = hsp.HealSparseMap.read(map_path).generate_healpix_map(nest=False)
	fs_map[fs_map == hp.UNSEEN] = 0.

	if is_systmap:
		if mask is not None:
			mu = np.sum(fs_map[mask.vpix] * mask.mask[mask.vpix]) / mask.sum
			fs_map[mask.vpix] -= mu
		else:
			print('Could not correct systematics map; no MaskData provided.')
	
	if apply_mask:
		if mask is not None:
			fs_map *= mask.mask
		else:
			print('Could not apply mask to map; no MaskData provided.')
		

	return fs_map


def load_tomographic_maps(map_path, apply_mask=False, mask=None):
	'''
	Loads files containing maps split into tomographic bins (e.g. delta_g) and
	returns their pixel values in healPIX RING ordering. If told to, will also 
	multiply the map by the mask, in which case a MaskData object is required
	as input. 


	Parameters
	----------
	map_path: str
		Path to the map being read.

	apply_mask: bool
		If True, will perform element-wise multiplication by the mask (given as separate input).

	mask: MaskData
		Object containing the mask and any potentially relevant pre-computed values.

	Returns
	-------
	fs_maps: list
		List containing full-sky data (RING ordering) for each tomographic map.
	'''
	#empty list in which the full-sky maps will be stored
	fs_maps = []

	#load the HealSparse file
	hsp_map = hsp.HealSparseMap.read(map_path)

	#cycle through the maps
	for d in hsp_map.dtype.names:
		#create full-sky realisation of the map
		fs_map = hsp_map[d].generate_healpix_map(nest=False)
		fs_map[fs_map == hp.UNSEEN] = 0.
		
		#multiply by the mask if told to do so
		if apply_mask:
			if mask is not None:
				fs_map *= mask.mask
			else:
				print('Could not apply mask to map; no MaskData provided.')
		
		#append the full-sky map to the list
		fs_maps.append(fs_map)
	
	return fs_maps






#######################################################
###############    START OF SCRIPT    #################
#######################################################

#maximum ell allowed by the resolution
ell_max = 3 * cf.nside_hi
#get pixel area in units of steradians
Apix = hp.nside2pixarea(cf.nside_hi)
#get the number of pixels in a full-sky map at the required resolution
npix = hp.nside2npix(cf.nside_hi)


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

#retrieve the number of redshift bins
nbins = len(cf.zbins) - 1
#also get all possible pairings of bins
l = list(range(nbins))
pairings = [i for i in itertools.product(l,l) if tuple(reversed(i)) >= i]


#cycle through the fields being analysed (TODO: later change to global fields)
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))

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
	
	#add the boolean 'lite' to the end of the list of systematics
	cf.systs.append(str(cf.lite))

	#file containing list of systematics maps deprojected in the previous run
	deproj_file = PATH_CACHE + cf.deproj_file
	if os.path.exists(deproj_file):
		with open(deproj_file, 'r+') as df:
			#see which (if any) systematics have been deprojected previously
			deproj_done = df.read().split('\n')
			#see if this is the same as the list specified in the config file (accounting for different ordering)
			if sorted(deproj_done) == sorted(cf.systs):
				print(f'Same systematics maps provided; skipping all calculations for field {fd}')
				continue
			else:
				print('Different systematics maps provided')
				#write the list of provided systematics to the file
				df.seek(0)
				df.truncate()
				df.write('\n'.join(cf.systs))
	else:
		if len(cf.systs) == 1:
			print('No systematics provided')
		else:
			with open(deproj_file, 'w') as df:
				df.write('\n'.join(cf.systs))
		
	#see if workspaces have already been created from a previous run
	wsp_path = PATH_CACHE + cf.wsp_file
	covwsp_path = PATH_CACHE + cf.covwsp_file
	if os.path.exists(wsp_path):
		w.read_from(wsp_path)
		calc = False
	else:
		calc = True
	if os.path.exists(covwsp_path):
		cw.read_from(covwsp_path)
		calc |= False
	else:
		calc |= True
		

	#load the delta_g maps
	deltag_maps = load_tomographic_maps(PATH_MAPS + cf.deltag_maps)

	#load the survey mask and convert to full-sky realisation
	mask = MaskData(PATH_MAPS + cf.survey_mask, mask_thresh=cf.weight_thresh)


	print('Loading systematics maps...')
	if len(cf.systs) > 1:
		#load the systematics maps and convert to full-sky realisations
		systmaps = [load_map(PATH_SYST + s, is_systmap=True, mask=mask) for s in cf.systs[:-1]]
		#reshape the resultant list to have dimensions (nsyst, 1, npix)
		nsyst = len(systmaps)
		systmaps = np.array(systmaps).reshape([nsyst, 1, npix])
		deproj = True
		print('templates: ', np.mean(systmaps))
	else:
		systmaps = None
		deproj = False
	print('Done!')


	print('Creating NmtFields...')
	density_fields = [nmt.NmtField(mask.mask, [d], templates=systmaps, lite=cf.lite) for d in deltag_maps]
	print('Done!')

	#delete the systematics and delta_g maps to clear some memory
	del systmaps
	del deltag_maps

	

	#retrieve the IDs of pixels above the mask threshold, as this is all that is
	#required from the mask henceforth
	above_thresh = mask.vpix
	sum_w_above_thresh = mask.sum
	mu_w = mask.mean
	mu_w2 = mask.meansq
	del mask
	
	#load the N_g maps and calculate the mean weighted by the mask
	mu_N_all = [nmap[above_thresh].sum() / sum_w_above_thresh for nmap in load_tomographic_maps(PATH_MAPS + cf.ngal_maps)]


	#full path to the output file
	outfile = f'{cf.PATH_OUT}{fd}/{cf.outfile}'
	#open the file, creating it if it doesn't exist
	with h5py.File(outfile, mode='r+') as psfile:
		
		#cycle through all possible pairings of redshift bins
		for ip,p in enumerate(pairings):
			i,j = p

			#see if a group for the current pairing already exists
			p_str = str(p)
			print(colour_string(p_str, 'green'))

			#see if group already exists for this pairing
			gp = psfile.require_group(p_str)

			f_i = density_fields[i]
			f_j = density_fields[j]
			cl_coupled = nmt.compute_coupled_cell(f_i, f_j)

			#use these along with the mask to get a guess of the true C_ell
			cl_guess = cl_coupled / mu_w2

			if ip == 0 and calc:
				#compute the mode coupling matrix (only need to compute once since same mask used for everything)
				print('Computing mode coupling matrix...')
				w.compute_coupling_matrix(f_i, f_j, b)
				print('Done!')

				print('Calculating coupling coefficients...')
				#compute coupling coefficients
				cw.compute_coupling_coefficients(f_i, f_j)
				print('Done!')
			else:
				print('Using coupling matrix and coefficients from cache.')

			#only calculate bias-related quantities if templates have been provided
			if deproj and not cf.lite:
				print('Calculating deprojection bias...')
				#compute the deprojection bias
				cl_bias = nmt.deprojection_bias(f_i, f_j, cl_guess)
				print('Done!')
			else:
				print('No systematics maps provided; skipping deprojection bias calculation.')
				cl_bias = np.zeros_like(cl_guess)
			
			#delete the group if it exists
			del psfile[p_str]

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
			gp = psfile.create_group(p_str)
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


	######################
	# CACHING WORKSPACES #
	######################

	#write the workspaces to the cache directory
	w.write_to(f'{PATH_CACHE}{cf.wsp_file}')
	cw.write_to(f'{PATH_CACHE}{cf.covwsp_file}')
		



