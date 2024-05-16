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
		outfile_now = f'{outfile_main[:-5]}_{i}_{j}.hdf5'

		f_i = density_fields[i]
		f_j = density_fields[j]
		#compute c_ells without deprojection
		cl_coupled_ii = nmt.compute_coupled_cell(f_i, f_i)
		cl_coupled_ij = nmt.compute_coupled_cell(f_i, f_j)
		cl_coupled_jj = nmt.compute_coupled_cell(f_j, f_j)
		#estimate true power spectra
		cl_guess_ii = cl_coupled_ii / mu_w2
		cl_guess_ij = cl_coupled_ij / mu_w2
		cl_guess_jj = cl_coupled_jj / mu_w2

		#decouple the c_ells (note this is without deprojection)
		cl_decoupled = w.decouple_cell(cl_coupled_ij)

		#recompute caussian covariance matrix
		print('Calculating covariance matrix...')
		covar = nmt.gaussian_covariance(cw,
								  		0, 0, 0, 0,
										[cl_guess_ii[0]],
										[cl_guess_ij[0]],
										[cl_guess_ij[0]],
										[cl_guess_jj[0]],
										w)
		err_cell = np.diag(covar) ** 0.5
		print('Done!')


		with h5py.File(outfile_now, 'a') as psfile:
			#populate the output file with the new results
			gp = psfile[p]
			_ = gp.create_dataset('cl_coupled_no_deproj', data=cl_coupled_ij)
			_ = gp.create_dataset('cl_decoupled_no_deproj', data=cl_decoupled)
			gp['covar'][...] = covar
			gp['err_cell'][...] = err_cell
	