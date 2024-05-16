#####################################################################################################
# - Uses NaMaster to compute power spectra from the galaxy delta_g maps, deprojecting any 
#   systematics templates in the process.
#####################################################################################################

import config
import healpy as hp
import healsparse as hsp
import numpy as np
from map_utils import *
import h5py
import pymaster as nmt
from output_utils import colour_string
import os
import sys
import glob

import faulthandler
faulthandler.enable()
### SETTINGS ###
cf = config.computePowerSpectra



#############################
######### FUNCTIONS #########
#############################

def make_density_fields(deproj_file, systs, idx=None):
	if os.path.exists(deproj_file):
		with open(deproj_file, 'r+') as df:
			#see which (if any) systematics have been deprojected previously
			deproj_done = df.read().split('\n')
			#see if this is the same as the list specified in the config file (accounting for different ordering)
			if sorted(deproj_done) == sorted(systs):
				print(f'Same systematics maps provided; skipping all calculations for field {fd}')
				return None, None
			else:
				print('Different systematics maps provided')
				#write the list of provided systematics to the file
				df.seek(0)
				df.truncate()
				df.write('\n'.join(systs))
	else:
		if len(systs) == 1:
			print('No systematics provided')
		else:
			with open(deproj_file, 'w') as df:
				df.write('\n'.join(systs))
		

	#load the delta_g maps
	deltag_maps = load_tomographic_maps(PATH_MAPS + cf.deltag_maps, idx=idx)


	print('Loading systematics maps...')
	if len(systs) > 1:
		#load the systematics maps and convert to full-sky realisations
		systmaps = [load_map(PATH_SYST + s, is_systmap=True, mask=mask) for s in systs[:-1]]
		#reshape the resultant list to have dimensions (nsyst, 1, npix)
		nsyst = len(systmaps)
		systmaps = np.array(systmaps).reshape([nsyst, 1, npix])
		print('templates: ', np.mean(systmaps))
	else:
		systmaps = None
	print('Done!')


	print('Creating NmtFields...')
	density_fields = [nmt.NmtField(mask.mask, [d], templates=systmaps, lite=cf.lite) for d in deltag_maps]
	print('Done!')

	#delete the systematics and delta_g maps to clear some memory
	del systmaps
	del deltag_maps

	return density_fields, nsyst


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#retrieve the pairing being analysed from the arguments if provided
try:
	pairings = [sys.argv[1]]
	per_tomo = True
except IndexError:
	_, pairings = cf.get_bin_pairings()
	per_tomo = False

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


#cycle through the fields being analysed (TODO: later change to global fields)
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))

	#path to the directory containing the maps
	PATH_MAPS = f'{cf.PATH_OUT}{fd}/'

	#load the survey mask and convert to full-sky realisation
	mask = MaskData(PATH_MAPS + cf.survey_mask, mask_thresh=cf.weight_thresh)
	#retrieve relevant quantities from the mask data
	above_thresh = mask.vpix
	sum_w_above_thresh = mask.sum
	mu_w = mask.mean
	mu_w2 = mask.meansq

	#set up a pymaster Workspace object
	w = nmt.NmtWorkspace()
	#create a variable assignment that will later be occupied by a CovarianceWorkspace
	cw = nmt.NmtCovarianceWorkspace()

	PATH_CACHE = PATH_MAPS + 'cache/'
	#see if directory for cached workspaces exists; make it if not
	if not os.path.exists(PATH_CACHE):
		os.system(f'mkdir -p {PATH_CACHE}')
	
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
	
	#path to directory containing systematics maps
	PATH_SYST = f'{PATH_MAPS}systmaps/'
	#check for 'All' in systmaps and convert this to a list of all systematics maps
	if 'all' in map(str.lower, cf.systs):
		systs = [os.path.basename(m) for m in (glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp'))]
	
	#if given a max number of systematics to deproject, slice the list accordingly
	if cf.Nsyst_max is not None:
		systs = systs[:cf.Nsyst_max]
	
	#add the boolean 'lite' to the end of the list of systematics
	systs.append(str(cf.lite))

	#file containing list of systematics maps deprojected in the previous run
	deproj_file = PATH_CACHE + cf.deproj_file
	if not per_tomo:
		density_fields, nsyst = make_density_fields(deproj_file, systs)
		if density_fields is None:
			continue

	#full path to the output file
	outfile_main = f'{cf.PATH_OUT}{fd}/{cf.outfile}'	
	for p in pairings:
		i,j = [int(x) for x in p.strip('()').split(',')]
		outfile_now = f'{outfile_main[:-5]}_{i}_{j}.hdf5'

		print(colour_string(p, 'green'))

		if per_tomo:
			deproj_file = f'{deproj_file[:-4]}_{i}_{j}.txt'
			density_fields, nsyst = make_density_fields(deproj_file, systs, idx=[i,j])
			if density_fields is None:
				continue
			f_i, f_j = density_fields
		else:
			f_i = density_fields[i]
			f_j = density_fields[j]

		cl_coupled = nmt.compute_coupled_cell(f_i, f_j)

		#use these along with the mask to get a guess of the true C_ell
		cl_guess = cl_coupled / mu_w2

		if calc:
			#compute the mode coupling matrix (only need to compute once since same mask used for everything)
			print('Computing mode coupling matrix...')
			w.compute_coupling_matrix(f_i, f_j, b)
			print('Done!')
			#write the workspace to the cache directory
			w.write_to(f'{PATH_CACHE}{cf.wsp_file}')

			print('Calculating coupling coefficients...')
			#compute coupling coefficients
			cw.compute_coupling_coefficients(f_i, f_j)
			#write the workspace to the cache directory
			cw.write_to(f'{PATH_CACHE}{cf.covwsp_file}')
			print('Done!')

			#set calc to False for future iterations
			calc = False
		else:
			print('Using coupling matrix and coefficients from cache.')

		#only calculate bias-related quantities if templates have been provided
		if (nsyst > 0) and not cf.lite:
			print('Calculating deprojection bias...')
			#compute the deprojection bias
			cl_bias = nmt.deprojection_bias(f_i, f_j, cl_guess)
			print('Done!')
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
			#load the N_g map and calculate the mean weighted by the mask
			mu_N = load_tomographic_maps(PATH_MAPS + cf.ngal_maps, idx=i)[0][above_thresh].sum() / sum_w_above_thresh
			#calculate the noise power spectrum
			N_ell_coupled = np.full(ell_max, Apix * mu_w / mu_N).reshape((1,ell_max))
			#decouple
			N_ell_decoupled = w.decouple_cell(N_ell_coupled)
		else:
			N_ell_coupled = np.zeros((1, ell_max))
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
		
		with h5py.File(outfile_now, 'w') as psfile: 
			#populate the output file with the results
			gp = psfile.create_group(p)
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
		



