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
from mpi4py import MPI

### SETTINGS ###
cf = config.computePowerSpectra

if not cf.LOCAL and cf.NERSC:
	os.system(f'taskset -pc 0-255 {os.getpid()}')

#############################
######### FUNCTIONS #########
#############################

def load_systematics(deproj_file, systs):
	'''
	Determines whether the specified list of systematics has already been
	deprojected on a previous run, and if not then returns an array containing
	all of the systematics full-sky maps. If these systematics have been
	deprojected on a previous run, returns an array on NaNs with the same shape. 

	Parameters
	----------
	deproj_file: str
		Path to the file which -- if it exists -- contains a list of all
		systematics deprojected on a previous run.
	
	systs: list
		List of systematics to be deprojected on this run.
	
	Returns
	-------
	systmaps: numpy.ndarray
		Array containing all of the systematics full-sky maps.
	'''
	#see if deproj_file already exists
	if os.path.exists(deproj_file):
		with open(deproj_file, 'r+') as df:
			#see which (if any) systematics have been deprojected previously
			deproj_done = df.read().split('\n')
			#see if this is the same as the list specified in the config file (accounting for different ordering)
			if sorted(deproj_done) == sorted(systs):
				print(f'Same systematics maps provided; skipping all calculations for field {fd}')
				systmaps = np.full((nsyst, 1, npix), np.nan, dtype=np.float64)
				return systmaps
			else:
				print('Different systematics maps provided')
				#write the list of provided systematics to the file
				df.seek(0)
				df.truncate()
				df.write('\n'.join(systs))
	else:
		if nsyst == 0:
			print('No systematics provided')
		else:
			with open(deproj_file, 'w') as df:
				df.write('\n'.join(systs))

	print('Loading systematics maps...')
	if nsyst > 0:
		#load the systematics maps and convert to full-sky realisations
		systmaps = [load_map(PATH_SYST + s, is_systmap=True, mask=mask) for s in systs[:-1]]
		#reshape the resultant list to have dimensions (nsyst, 1, npix)
		systmaps = np.array(systmaps).reshape([nsyst, 1, npix])
		print('templates: ', np.mean(systmaps))
	else:
		systmaps = np.array(np.nan, dtype=np.float64)
	print('Done!')

	return systmaps


def compute_covariance(w, cw, f_i1, f_i2, f_j1=None, f_j2=None, return_cl_coupled=False, return_cl_guess=False):
	'''
	Computes the covariance matrix for a pair of power spectra. Can also
	optionally return the coupled C_ells used to calculate the covariances,
	as well as the guess C_ells computed from these.

	Parameters
	----------
	w: nmt.NmtWorkspace
		Workspace containing the mode coupling matrix (computed prior to this).
	
	cw: nmt.NmtCovarianceWorkspace
		Covariance workspace containing the coupling coefficients (computed 
		prior to this).
	
	f_i1, f_i2: nmt.NmtField, nmt.NmtField
		Fields involved in the calculation of the first power spectrum.
	
	f_j1, f_j2: nmt.NmtField, nmt.NmtField
		Fields involved in the calculation of the second power spectrum. If None,
		the calculated covariance will be that for the first power spectrum cross-
		correlated with itself.
	
	return_cl_coupled: bool
		If True, will return the coupled C_ells used to compute the covariance matrix.
	
	return_cl_guess: bool
		If True, will return the guess C_ells estimated from the coupled C_ells.
	
	Returns
	-------
	covar: numpy.ndarray
		Covariance matrix.
	
	err_cell: numpy.ndarray
		Square root of the diagonal of the covariance matrix (corresponds to the
		uncertainties on the power spectra data points).
	'''
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

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

#all possible pairings of tomographic bins (as tuples and as strings)
pairings, pairings_s = cf.get_bin_pairings()

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


#cycle through the fields being analysed
for fd in cf.get_global_fields():
	#path to the directory containing the maps
	PATH_MAPS = f'{cf.PATH_OUT}{fd}/'
	#path to directory of cached outputs from this script
	PATH_CACHE = PATH_MAPS + 'cache/'
	wsp_path = PATH_CACHE + cf.wsp_file
	covwsp_path = PATH_CACHE + cf.covwsp_file
	#set up NmtWorkspace and NmtCovarianceWorkspace
	w = nmt.NmtWorkspace()
	cw = nmt.NmtCovarianceWorkspace()
	
	if rank == 0:
		print(colour_string(fd.upper(), 'orange'))

		#see if directory for cached workspaces exists; make it if not
		if not os.path.exists(PATH_CACHE):
			os.system(f'mkdir -p {PATH_CACHE}')

		#load the survey mask and convert to full-sky realisation
		mask = MaskData(PATH_MAPS + cf.survey_mask)
		#temporarily create an NmtField using just the mask
		fmask = nmt.NmtField(mask.mask_full, maps=None, spin=0)

		#see if workspaces have already been created from a previous run
		if os.path.exists(wsp_path):
			w.read_from(wsp_path)
		else:
			print('Computing mode coupling matrix...')
			w.compute_coupling_matrix(fmask, fmask, b)
			print('Done!')
			#write the workspace to the cache directory
			w.write_to(wsp_path)
		#either load or compute coupling coefficients
		if os.path.exists(covwsp_path):
			cw.read_from(covwsp_path)
		else:
			print('Calculating coupling coefficients...')
			cw.compute_coupling_coefficients(fmask, fmask)
			print('Done!')
			#write the workspace to the cache directory
			cw.write_to(covwsp_path)
		
		#delete the NmtField to conserve memory
		del fmask

		#path to directory containing systematics maps
		PATH_SYST = f'{PATH_MAPS}systmaps/'
		systs = []
		#check for 'All' in systmaps and convert this to a list of all systematics maps
		if 'all' in map(str.lower, cf.systs):
			systs = [os.path.basename(m) for m in (glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp'))]
		#if given a max number of systematics to deproject, slice the list accordingly
		if cf.Nsyst_max is not None:
			systs = systs[:cf.Nsyst_max]
		#get the number of systematics
		nsyst = len(systs)
		#add the boolean 'lite' to the end of the list of systematics
		systs.append(str(cf.lite))
		#file containing list of systematics maps deprojected in the previous run
		deproj_file = PATH_CACHE + cf.deproj_file
		#load the systematics maps
		systmaps = load_systematics(deproj_file, systs)
		
	else:
		mask = None
		nsyst = None
		while True:
			try:
				w.read_from(wsp_path)
				cw.read_from(covwsp_path)
				break
			except (FileNotFoundError, RuntimeError):
				continue

	#broadcast the mask 
	mask = comm.bcast(mask, root=0)
	#retrieve relevant quantities from the mask data
	above_thresh = mask.vpix_ring
	sum_w_above_thresh = mask.sum
	mu_w = mask.mean
	mu_w2 = mask.meansq	

	#broadcast the number of systematics
	nsyst = comm.bcast(nsyst, root=0)
	#broadcast the array of systematics maps
	if rank != 0:
		systmaps = np.empty((nsyst, 1, npix), dtype=np.float64)
	comm.Bcast(systmaps, root=0)
	#if systmaps is all NaN, continue to next field
	if np.isnan(systmaps).all():
		continue

	#create density fields, with and without deprojection, splitting the work between ranks if possible
	if size >= cf.nbins:
		if rank in range(cf.nbins):
			dg_map, = load_tomographic_maps(PATH_MAPS + cf.deltag_maps, idx=rank)
			print(f'Creating NmtField for bin {rank} (without deprojection)...')
			density_fields_nd = nmt.NmtField(mask.mask_full, [dg_map], templates=None, lite=cf.lite)
			print(f'Creating NmtField for bin {rank} (with deprojection)...')
			density_fields = nmt.NmtField(mask.mask_full, [dg_map], templates=systmaps, lite=cf.lite)
		density_fields_nd = comm.gather(density_fields_nd, root=0)
		density_fields = comm.gather(density_fields, root=0)
	else:
		if rank == 0:
			dg_map = load_tomographic_maps(PATH_MAPS + cf.deltag_maps)
			print(f'Creating NmtFields for all bins (without deprojection)...')
			density_fields_nd = [nmt.NmtField(mask.mask_full, [dg], templates=None, lite=cf.lite)
									for dg in dg_map]
			print(f'Creating NmtField for all bins (with deprojection)...')
			density_fields = [nmt.NmtField(mask.mask_full, [dg], templates=systmaps, lite=cf.lite)
					 				for dg in dg_map]
	
	#broadcast the list of NmtFields
	if rank != 0:
		density_fields_nd = None
		density_fields = None
	density_fields_nd = comm.bcast(density_fields_nd, root=0)
	density_fields = comm.bcast(density_fields, root=0)

	#compute the C_ells, splitting the work across ranks if possible
	if size >= len(pairings):
		#current bin pairing
		i, j = pairings[rank]
		#current density fields (no deprojection)
		f_i_nd = density_fields_nd[i]
		f_j_nd = density_fields_nd[j]
		#current density fields (with deprojection)
		f_i = density_fields[i]
		f_j = density_fields[j]

		print(f'Calculating covar(C^{i}{j}, C^{i}{j}) (without deprojection)...')
		covar_nd, err_cell_nd, (_, cl_coupled_ij_nd, _, _), (_, cl_guess_ij_nd, _, _) = compute_covariance(w, cw, f_i_nd, f_j_nd,
																					return_cl_coupled=True,
																					return_cl_guess=True)
		print('Done!')

		print(f'Calculating covar(C^{i}{j}, C^{i}{j}) (with deprojection)...')
		if nsyst > 0:
			covar, err_cell, (_, cl_coupled_ij, _, _), (_, cl_guess_ij, _, _) = compute_covariance(w, cw, f_i, f_j,
																					return_cl_coupled=True,
																					return_cl_guess=True)
		else:
			covar = covar_nd
			err_cell = err_cell_nd
			cl_coupled_ij = cl_coupled_ij_nd
			cl_guess_ij = cl_guess_ij_nd
		print('Done!')
	else:
		if rank == 0:
			for p in pairings:
				i, j = p
			#current density fields (no deprojection)
			f_i_nd = density_fields_nd[i]
			f_j_nd = density_fields_nd[j]
			#current density fields (with deprojection)
			f_i = density_fields[i]
			f_j = density_fields[j]

			print(f'Calculating covar(C^{i}{j}, C^{i}{j}) (without deprojection)...')
			covar_nd, err_cell_nd, (_, cl_coupled_ij_nd, _, _), (_, cl_guess_ij_nd, _, _) = compute_covariance(w, cw, f_i_nd, f_j_nd,
																						return_cl_coupled=True,
																						return_cl_guess=True)
			print('Done!')

			print(f'Calculating covar(C^{i}{j}, C^{i}{j}) (with deprojection)...')
			if nsyst > 0:
				covar, err_cell, (_, cl_coupled_ij, _, _), (_, cl_guess_ij, _, _) = compute_covariance(w, cw, f_i, f_j,
																						return_cl_coupled=True,
																						return_cl_guess=True)
			else:
				covar = covar_nd
				err_cell = err_cell_nd
				cl_coupled_ij = cl_coupled_ij_nd
				cl_guess_ij = cl_guess_ij_nd
			print('Done!')
