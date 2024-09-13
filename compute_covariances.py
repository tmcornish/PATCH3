#####################################################################################################
# - Uses outputs from compute_power_spectra.py to compute the covariance matrices
#   for all combinations of power spectra.
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

### SETTINGS ###
cf = config.computeCovariances


#############################
######### FUNCTIONS #########
#############################

def compute_covariance(w, cw, f_i1, f_i2, f_j1=None, f_j2=None, return_cl_coupled=False, return_cl_guess=False):
	#see if additional fields have been provided for the second power spectrum
	if f_j1 is None:
		f_j1 = f_i1
	if f_j2 is None:
		f_j2 = f_i2
	
	#compute coupled c_ells for each possible combination of i and j
	cl_coupled_i1j1 = nmt.compute_coupled_cell(f_i1, f_j1)
	cl_coupled_i1j2 = nmt.compute_coupled_cell(f_i1, f_j2)
	cl_coupled_i2j1 = nmt.compute_coupled_cell(f_i2, f_j1)
	cl_coupled_i2j2 = nmt.compute_coupled_cell(f_i2, f_j2)
	#use these along with the mask to get a guess of the true C_ell
	cl_guess_i1j1 = cl_coupled_i1j1 / mu_w2
	cl_guess_i1j2 = cl_coupled_i1j2 / mu_w2
	cl_guess_i2j1 = cl_coupled_i2j1 / mu_w2
	cl_guess_i2j2 = cl_coupled_i2j2 / mu_w2


	covar = nmt.gaussian_covariance(cw, 
									0, 0, 0, 0,			#spin of each field
									[cl_guess_i1j1[0]],	
									[cl_guess_i1j2[0]],
									[cl_guess_i2j1[0]],
									[cl_guess_i2j2[0]],
									w)
	#errorbars for each bandpower
	err_cell = np.diag(covar) ** 0.5

	to_return = [covar, err_cell]
	if return_cl_coupled:
		to_return.append([cl_coupled_i1j1, cl_coupled_i1j2, cl_coupled_i2j1, cl_coupled_i2j2])
	if return_cl_guess:
		to_return.append([cl_guess_i1j1, cl_guess_i1j2, cl_guess_i2j1, cl_guess_i2j2])


	return (*to_return,)


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#retrieve the number of tomographic bins and number of pixels per map
nbins = len(cf.zbins) - 1
npix = hp.nside2npix(cf.nside_hi)

#cycle through the fields being analysed
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))

	#path to the directory containing the maps
	PATH_MAPS = f'{cf.PATH_OUT}{fd}/'
	#systematics maps
	PATH_SYST = PATH_MAPS + 'systmaps/'
	#path to the cache for this field
	PATH_CACHE = PATH_MAPS + 'cache/'

	#load the survey mask and convert to full-sky realisation
	mask = MaskData(PATH_MAPS + cf.survey_mask)
	#retrieve relevant quantities from the mask data
	above_thresh = mask.vpix_ring
	sum_w_above_thresh = mask.sum
	mu_w = mask.mean
	mu_w2 = mask.meansq

	#load the (non-deprojected) delta_g maps
	deltag_maps = load_tomographic_maps(PATH_MAPS + cf.deltag_maps)
	#load the systematics maps
	deproj_file = PATH_CACHE + cf.deproj_file[:-4] + '_0_0.txt'
	with open(deproj_file, 'r') as depfile:
		deproj_done = depfile.read().split('\n')[:-1]
	systmaps = [load_map(PATH_SYST + s, is_systmap=True, mask=mask) for s in deproj_done]
	nsyst = len(systmaps)
	systmaps = systmaps.reshape(nsyst, 1, npix)
	#mask the systmaps
	systmaps *= mask.mask_full
	#construct density fields with and w/o deprojection
	density_fields, density_fields_nd = [], []
	for i in range(len(deltag_maps)):
		df_nd = nmt.NmtField(mask.mask_full, [deltag_maps[i]], templates=None)
		density_fields_nd.append(df_nd)
		#for the deprojected fields, load the alphas (if they exist)
		alpha_file = PATH_CACHE + cf.alphas_file[:-4] + f'bin{i}.txt'
		#create a copy of the non-deprojected deltag map
		dg_deproj = deltag_maps[i].copy()
		#multiply by the mask
		dg_deproj *= mask.mask_full
		#reshape for compatibility with other arrays
		dg_deproj = dg_deproj.reshape(1, npix)
		if os.path.exists(alpha_file):
			alphas = np.loadtxt(alpha_file)
			#deproject the systematics
			dg_deproj -= np.sum(alphas[:, None, None] * systmaps, axis=0)
			#construct the NmtField
			df = nmt.NmtField(np.ones_like(mask.mask_full), dg_deproj, templates=None)
			density_fields.append(df)

	#load the workspaces from the cache
	w = nmt.NmtWorkspace()
	w.read_from(PATH_CACHE + cf.wsp_file)
	cw = nmt.NmtCovarianceWorkspace()
	cw.read_from(PATH_CACHE + cf.covwsp_file)

	#file for containing the covariance matrices
	covar_file = f'{PATH_MAPS}{cf.covar_file}'
	#number of different fields
	nfields = len(cf.zbins) - 1
	with h5py.File(covar_file, 'w') as cvfile:
		#cycle through all possible combinations of pairs of fields
		for i1 in range(nfields):
			for i2 in range(i1, nfields):
				for j1 in range(nfields):
					for j2 in range(j1, nfields):
						#create group in the hdf5 file for this pairing
						gp = cvfile.create_group(f'{i1}{i2}-{j1}{j2}')
						#compute the covariance
						covar_now, *_ = compute_covariance(w, cw, 
										density_fields[i1],
										density_fields[i2],
										density_fields[j1],
										density_fields[j2]
										)
						if len(density_fields) > 0:
							covar_nd_now, *_ = compute_covariance(w, cw, 
										density_fields_nd[i1],
										density_fields_nd[i2],
										density_fields_nd[j1],
										density_fields_nd[j2]
										)
						else:
							covar_nd_now = covar_now
						#create datasets for the covariance matrices
						_ = gp.create_dataset('final', data=covar_now)
						_ = gp.create_dataset('no_deproj', data=covar_nd_now)