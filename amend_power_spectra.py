import config
import healpy as hp
import numpy as np
import h5py
import pymaster as nmt
from map_utils import *
from output_utils import colour_string

cf = config.computePowerSpectra

#######################################################
###############    START OF SCRIPT    #################
#######################################################

#retrieve the pairing being analysed from the arguments if provided
_, pairings = cf.get_bin_pairings()

#maximum ell allowed by the resolution
ell_max = 3 * cf.nside_hi
#get pixel area in units of steradians
Apix = hp.nside2pixarea(cf.nside_hi)
#get the number of pixels in a full-sky map at the required resolution
npix = hp.nside2npix(cf.nside_hi)


#cycle through the fields being analysed (TODO: later change to global fields)
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))

	#path to the directory containing the maps
	PATH_MAPS = f'{cf.PATH_OUT}{fd}/'
	PATH_CACHE = PATH_MAPS + 'cache/'

	#load the survey mask and convert to full-sky realisation
	mask = MaskData(PATH_MAPS + cf.survey_mask, mask_thresh=cf.weight_thresh)
	#retrieve relevant quantities from the mask data
	above_thresh = mask.vpix
	sum_w_above_thresh = mask.sum
	mu_w = mask.mean
	mu_w2 = mask.meansq

	#load the delta_g maps
	deltag_maps = load_tomographic_maps(PATH_MAPS + cf.deltag_maps)

	print('Creating NmtFields...')
	density_fields = [nmt.NmtField(mask.mask, [d], templates=None, lite=cf.lite) for d in deltag_maps]
	print('Done!')
	del deltag_maps

	#read Workspaces form file
	w = nmt.NmtWorkspace()
	w.read_from(PATH_CACHE + cf.wsp_file)
	#create a variable assignment that will later be occupied by a CovarianceWorkspace
	cw = nmt.NmtCovarianceWorkspace()
	cw.read_from(PATH_CACHE + cf.covwsp_file)

	#full path to the output file
	outfile_main = f'{cf.PATH_OUT}{fd}/{cf.outfile}'
	for p in pairings:
		i,j = [int(x) for x in p.strip('()').split(',')]
		outfile_ii = f'{outfile_main[:-5]}_{i}_{i}.hdf5'
		outfile_ij = f'{outfile_main[:-5]}_{i}_{j}.hdf5'
		outfile_jj = f'{outfile_main[:-5]}_{j}_{j}.hdf5'

		f_i = density_fields[i]
		f_j = density_fields[j]
		#compute c_ells without deprojection
		cl_coupled_ii = nmt.compute_coupled_cell(f_i, f_i)
		cl_coupled_nd = nmt.compute_coupled_cell(f_i, f_j)
		cl_coupled_jj = nmt.compute_coupled_cell(f_j, f_j)
		#estimate true power spectra
		cl_guess_ii = cl_coupled_ii / mu_w2
		cl_guess_ij = cl_coupled_nd / mu_w2
		cl_guess_jj = cl_coupled_jj / mu_w2

		#decouple the c_ells (no deprojection)
		cl_decoupled_nd = w.decouple_cell(cl_coupled_nd)

		#recompute caussian covariance matrix
		print('Calculating covariance matrix (no deprojection)...')
		covar_nd = nmt.gaussian_covariance(cw,
								  		0, 0, 0, 0,
										[cl_guess_ii[0]],
										[cl_guess_ij[0]],
										[cl_guess_ij[0]],
										[cl_guess_jj[0]],
										w)
		err_cell_nd = np.diag(covar_nd) ** 0.5
		print('Done!')


		#load the coupled c_ells with deprojection
		with h5py.File(outfile_ii, 'r') as psfile:
			cl_coupled_ii = psfile[f'{i},{i}/cl_coupled'][...]
		with h5py.File(outfile_ij, 'r') as psfile:
			cl_coupled = psfile[f'{i},{j}/cl_coupled'][...]
		with h5py.File(outfile_jj, 'r') as psfile:
			cl_coupled_jj = psfile[f'{j},{j}/cl_coupled'][...]
		#use these to estimate true c_ells with deprojection
		cl_guess_ii = cl_coupled_ii / mu_w2
		cl_guess_ij = cl_coupled / mu_w2
		cl_guess_jj = cl_coupled_jj / mu_w2
		#recompute caussian covariance matrix
		print('Calculating covariance matrix (with deprojection)...')
		covar = nmt.gaussian_covariance(cw,
								  		0, 0, 0, 0,
										[cl_guess_ii[0]],
										[cl_guess_ij[0]],
										[cl_guess_ij[0]],
										[cl_guess_jj[0]],
										w)
		err_cell = np.diag(covar) ** 0.5
		print('Done!')


		with h5py.File(outfile_ij, 'a') as psfile:
			#populate the output file with the new results
			gp = psfile[p]
			_ = gp.create_dataset('cl_coupled_no_deproj', data=cl_coupled_nd)
			_ = gp.create_dataset('cl_decoupled_no_deproj', data=cl_decoupled_nd)
			_ = gp.create_dataset('covar_no_deproj', data=covar_nd)
			_ = gp.create_dataset('err_cell_no_deproj', data=err_cell_nd)
			gp['covar'][...] = covar
			gp['err_cell'][...] = err_cell
