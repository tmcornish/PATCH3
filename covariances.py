#####################################################################################################
# - Uses the saved best-fit coefficients from linear deprojection to reconstruct the deprojected
# - density fields and compute the covariances for each power spectrum.
#####################################################################################################

import os
import sys
import healpy as hp
import numpy as np
import pymaster as nmt
from output_utils import colour_string
import pandas as pd
import h5py
from configuration import PipelineConfig as PC
from map_utils import MaskData, load_tomographic_maps, load_map
from cell_utils import get_bin_pairings

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='covariances')

#############################
######### FUNCTIONS #########
#############################


def make_deprojected_field(dg_map, alphas=None):
	'''
	Uses pre-saved coefficients for linear deprojection to reconstruct the systematics-deprojected
	delta_g map, and created an NmtField object from that.

	Parameters
	----------
	dg_map: numpy.ndarray
		Full-sky delta_g map, with no deprojection applied.
	
	alphas: pandas.DataFrame or None
		DataFrame containing the name of each systematic and its corresponding best-fit deprojection
		coefficient. If None, no deprojection will occur.

	Returns
	-------
	df: pymaster.NmtField
		NaMaster Field object constructed from the deprojected delta_g map.
	'''
	#mask the delta_g map and reshape it for namaster
	dg_map *= mask.mask_full
	dg_map = dg_map.reshape(1, npix)
	if alphas is not None:
		#load the systematics
		systmaps = [load_map(PATH_SYST + s, is_systmap=True, mask=mask) for s in alphas.index]
		nsyst = len(systmaps)
		systmaps = np.array(systmaps).reshape([nsyst, 1, npix])
		#apply the mask to the systematics
		systmaps *= mask.mask_full
		#deproject the systematics
		dg_map -= np.sum(alphas['alpha'].values[:, None, None] * systmaps, axis=0)

	#create an NmtField object with the deprojected map (do NOT re-mask)
	df = nmt.NmtField(np.ones_like(mask.mask_full), [dg_map], templates=None)

	return df


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

#get the bin pairings for which C_ells were computed in the previous stage
pairings, _ = get_bin_pairings(cf.nsamples)

#healpix parameters
npix = hp.nside2npix(cf.nside_hi)

#cycle through the fields being analysed
for fd in cf.fields:
	print(colour_string(fd.upper(), 'orange'))

	#path to the directory containing the outputs for this field
	PATH_FD = f'{cf.paths.out}{fd}/'
	#path to directory containing systematics maps
	PATH_SYST = PATH_FD + 'systmaps/'
	#directory of cached data from computing power spectra
	PATH_CACHE = PATH_FD + 'cache/'

	#load relevant workspaces from file
	wsp_path = PATH_CACHE + cf.cache_files.workspaces.wsp
	covwsp_path = PATH_CACHE + cf.cache_files.workspaces.covwsp
	w = nmt.NmtWorkspace()
	cw = nmt.NmtCovarianceWorkspace()
	w.read_from(wsp_path)
	cw.read_from(covwsp_path)

	#load the survey mask and convert to full-sky realisation
	mask = MaskData(PATH_FD + cf.maps.survey_mask)
	#retrieve relevant quantities from the mask data
	above_thresh = mask.vpix_ring
	sum_w_above_thresh = mask.sum
	mu_w = mask.mean
	mu_w2 = mask.meansq
	
	#load the delta_g maps
	deltag_maps = load_tomographic_maps(PATH_FD + cf.maps.deltag_maps)
	#create NmtFields without deprojection
	density_fields_nd = [nmt.NmtField(mask.mask_full, [dg], templates=None) for dg in deltag_maps]

	#load (as DataFrames) the best-fit coefficients
	alpha_dfs = []
	for i in range(cf.nsamples):
		alphas_path = PATH_FD + cf.cache_files.deproj.alphas[:-4] + f'_bin{i}.txt'
		if os.path.exists(alphas_path):
			alpha_dfs.append(pd.read_csv(alphas_path, sep='\t', index_col=0))
		else:
			alpha_dfs.append(None)
	#create NmtFields with deprojection
	density_fields = [make_deprojected_field(dg, al) for dg, al in zip(deltag_maps, alpha_dfs)]

	#file for containing the covariance matrices
	covar_file = f'{PATH_FD}{cf.cell_files.covariances}'
	with h5py.File(covar_file, 'w') as cvfile:
		#cycle through all possible combinations of pairs of fields
		for i1 in range(cf.nsamples):
			for i2 in range(i1, cf.nsamples):
				for j1 in range(cf.nsamples):
					for j2 in range(j1, cf.nsamples):
						#create group in the hdf5 file for this pairing
						gp = cvfile.create_group(f'{i1}{i2}-{j1}{j2}')
						#compute the covariance
						covar_now, *_ = compute_covariance(w, cw, 
										density_fields[i1],
										density_fields[i2],
										density_fields[j1],
										density_fields[j2]
										)
						covar_nd_now, *_ = compute_covariance(w, cw, 
									density_fields_nd[i1],
									density_fields_nd[i2],
									density_fields_nd[j1],
									density_fields_nd[j2]
									)
						#create datasets for the covariance matrices
						_ = gp.create_dataset('final', data=covar_now)
						_ = gp.create_dataset('no_deproj', data=covar_nd_now)
