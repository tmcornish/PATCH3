#####################################################################################################
# - Uses NaMaster to compute power spectra from the galaxy delta_g maps, deprojecting any 
#   systematics templates in the process.
#####################################################################################################

from configuration import PipelineConfig as PC
import healpy as hp
import healsparse as hsp
import numpy as np
from map_utils import load_map, load_tomographic_maps, MaskData
import cell_utils as cu
import h5py
import pymaster as nmt
from output_utils import colour_string
import os
import sys
import glob
from mpi4py import MPI
import sacc
import pandas as pd

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='computePowerSpectra')

if cf.platform == 'nersc':
	os.system(f'taskset -pc 0-255 {os.getpid()}')

#############################
######### FUNCTIONS #########
#############################

def load_systematics(deproj_file, systs):
	'''
	Determines whether the specified list of systematics has already been
	deprojected on a previous run, and if not then returns an array containing
	all of the systematics full-sky maps. If these systematics have been
	deprojected on a previous run, returns an array of NaNs with the same shape. 

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
	nsyst = len(systs) - 1
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

def split_list(a, n):
	'''
	Splits a list into n approximately equal parts. For example, a list
	of length 10 split into 3 parts would be returned as chunks of size
	(4, 3, 3).

	Parameters
	----------
	a: array-like
		List (or array) to be split.
	
	n: int
		Number of chunks into which the list is to be split.
	
	Returns
	-------
	a_split: list
		List containing the separate chunks of the input list.
	'''
	#compute the quotient and remainder of the list length divided by the number of chunks
	k, m = divmod(len(a), n)
	#split into n approximately equal-sized chunks
	a_split = [a[i*k+min(i,m):(i+1)*k+min(i+1,m)] for i in range(n)]
	return a_split


#######################################################
###############    START OF SCRIPT    #################
#######################################################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

#all possible pairings of tomographic bins (as tuples and as strings)
pairings, _, label_pairs = cu.get_bin_pairings(cf.nsamples, labels=list(cf.samples))
npairs = len(pairings)
#split the lists of pairings and labels depending on the comm size
pairings_split = split_list(pairings, size)
label_pairs_split = split_list(label_pairs, size)

#maximum ell allowed by the resolution
ell_max = 3 * cf.nside_hi
#get pixel area in units of steradians
Apix = hp.nside2pixarea(cf.nside_hi)
#get the number of pixels in a full-sky map at the required resolution
npix = hp.nside2npix(cf.nside_hi)

#set the bandpower edges
bpw_edges = cu.get_bpw_edges(ell_max, ell_min=cf.ell_min, nbpws=cf.nbpws, spacing=cf.bpw_spacing)
#create pymaster NmtBin object using these bandpower objects
b = nmt.NmtBin.from_edges(bpw_edges[:-1], bpw_edges[1:])
#get the effective ells
ell_effs = b.get_effective_ells()


#cycle through the fields being analysed
for fd in cf.fields:
	#path to the directory containing the maps
	PATH_FD = f'{cf.paths.out}{fd}/'
	#path to directory of cached outputs from this script
	PATH_CACHE = PATH_FD + 'cache/'
	wsp_path = PATH_CACHE + cf.cache_files.workspaces.wsp
	covwsp_path = PATH_CACHE + cf.cache_files.workspaces.covwsp
	#set up NmtWorkspace and NmtCovarianceWorkspace
	w = nmt.NmtWorkspace()
	cw = nmt.NmtCovarianceWorkspace()
	
	#load the survey mask and convert to full-sky realisation
	mask = MaskData(PATH_FD + cf.maps.survey_mask)
	#retrieve relevant quantities from the mask data
	above_thresh = mask.vpix_ring
	sum_w_above_thresh = mask.sum
	mu_w = mask.mean
	mu_w2 = mask.meansq	

	if rank == 0:
		print(colour_string(fd.upper(), 'orange'))

		#see if directory for cached workspaces exists; make it if not
		if not os.path.exists(PATH_CACHE):
			os.system(f'mkdir -p {PATH_CACHE}')

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
		PATH_SYST = f'{PATH_FD}systmaps/'
		systs = []
		#check for 'All' in systmaps and convert this to a list of all systematics maps
		if 'all' in map(str.lower, cf.systs):
			systs = [os.path.basename(m) for m in glob.glob(f'{PATH_SYST}*_nside{cf.nside_hi}*.hsp')]
		#if given a max number of systematics to deproject, slice the list accordingly
		if cf.Nsyst_max is not None:
			systs = systs[:cf.Nsyst_max]
		#add the boolean 'lite' to the end of the list of systematics
		systs.append(str(cf.lite))
		#file containing list of systematics maps deprojected in the previous run
		deproj_file = PATH_CACHE + cf.cache_files.deproj.deprojected
		#load the systematics maps
		systmaps = load_systematics(deproj_file, systs)
		
	else:
		nsyst = None
		systs = None
		while True:
			try:
				w.read_from(wsp_path)
				cw.read_from(covwsp_path)
				break
			except (FileNotFoundError, RuntimeError):
				continue

	#broadcast the list of systematics
	syst = comm.bcast(systs, root=0)
	#get the number of systematics
	nsyst = len(syst) - 1
	#see if no systematics have been provided
	if nsyst == 0:
		systmaps = None
	else:
		#broadcast the array of systematics maps
		if rank != 0:
			systmaps = np.empty((nsyst, 1, npix), dtype=np.float64)
		comm.Bcast(systmaps, root=0)
		#if systmaps is all NaN, continue to next field
		if np.isnan(systmaps).all():
			continue

	#number of iterations required per node to compute all bin pairs
	niter = int(np.ceil(npairs / size))
	#select a set of pairings according to the current rank
	pairings_now = pairings_split[rank]
	label_pairs_now = label_pairs_split[rank]
	#set up arrays for the c_ells calculated on the current node
	cls_now = np.full((niter, 1, cf.nbpws), np.nan)
	cls_nd_now = np.full((niter, 1, cf.nbpws), np.nan)
	cls_noise_now = np.full((niter, 1, cf.nbpws), np.nan)
	cls_bias_now = np.full((niter, 1, cf.nbpws), np.nan)

	#set up buffers for gathering results
	cl_buff = None			#buffer for the main C_ells
	cl_nd_buff = None		#buffer for the C_ells without deprojection
	cl_noise_buff = None 	#buffer for the noise power spectra
	cl_bias_buff = None		#buffer for the deprojection bias
	if rank == 0:
		cl_buff = np.empty((size, niter, 1, cf.nbpws), dtype=np.float64)			
		cl_nd_buff = np.empty((size, niter, 1, cf.nbpws), dtype=np.float64)					
		cl_noise_buff = np.empty((size, niter, 1, cf.nbpws), dtype=np.float64)			 	
		cl_bias_buff = np.empty((size, niter, 1, cf.nbpws), dtype=np.float64)

	#set up dictionaries for storing computed fields
	fdone_nd = {lp: None for lp in cf.samples}
	fdone = {lp: None for lp in cf.samples}
	#cycle through the pairings assigned to this node
	for idx in range(len(pairings_now)):
		#current bin pairing
		i, j = pairings_now[idx]
		label_i, label_j = label_pairs_now[idx]
		#check if field i has already been computed
		if fdone[label_i] is None:
			#delta_g map
			dg_i = load_tomographic_maps(PATH_FD + cf.maps.deltag_maps, idx=i)[0]

			print(f'Creating NmtField for {label_i} (without deprojection)...')
			##########################################################################
			f_i_nd = nmt.NmtField(mask.mask_full, [dg_i], templates=None, lite=cf.lite)
			
			if nsyst == 0:
				f_i = f_i_nd
			else:
				print(f'Creating NmtField for {label_i} (with deprojection)...')
				##########################################################################
				f_i = nmt.NmtField(mask.mask_full, [dg_i], templates=systmaps, lite=cf.lite)
			
			#store these in the dictionary iff they are required in later calculations on this node
			if any([label_i in p for p in label_pairs_now[idx:]]):
				fdone_nd[label_i] = f_i_nd
				fdone[label_i] = f_i
		else:
			#retrieve the fields from the dictionary rather than reconstruct them
			f_i_nd = fdone_nd[label_i]
			f_i = fdone[label_i]
		
		#do the same for field j
		if fdone[label_j] is None:
			#delta_g map
			dg_j = load_tomographic_maps(PATH_FD + cf.maps.deltag_maps, idx=j)[0]
			
			print(f'Creating NmtField for {label_j} (without deprojection)...')
			##########################################################################
			f_j_nd = nmt.NmtField(mask.mask_full, [dg_j], templates=None, lite=cf.lite)
			
			if nsyst == 0:
				f_j = f_j_nd
			else:
				print(f'Creating NmtField for {label_j} (with deprojection)...')
				##########################################################################
				f_j = nmt.NmtField(mask.mask_full, [dg_j], templates=systmaps, lite=cf.lite)

			#store these in the dictionary iff they are required in later calculations on this node
			if any([label_j in p for p in label_pairs_now[idx:]]):
				fdone_nd[label_j] = f_j_nd
				fdone[label_j] = f_j
		else:
			#retrieve the fields from the dictionary rather than reconstruct them
			f_j_nd = fdone_nd[label_j]
			f_j = fdone[label_j]

		print(f'Calculating coupled C_ells for pairing {label_i},{label_j}...')
		###########################################################
		#without deprojection
		cl_coupled_nd = nmt.compute_coupled_cell(f_i_nd, f_j_nd)
		cl_guess_nd = cl_coupled_nd / mu_w2

		#with deprojection
		if nsyst > 0:
			cl_coupled = nmt.compute_coupled_cell(f_i, f_j)
			cl_guess = cl_coupled / mu_w2
		else:
			cl_coupled = cl_coupled_nd
			cl_guess = cl_guess_nd

		if i == j:
			if nsyst > 0:
				print(f'Saving deprojection coefficients for {label_i}...')
				##############################################################
				alphas = f_i.alphas
				with open(PATH_CACHE + cf.cache_files.deproj.alphas[:-4] + f'_{label_i}.txt', 'w') as alphas_file:
					alphas_file.write('Sytematic\talpha\n')
					for k in range(nsyst):
						alphas_file.write(f'{syst[k]}\t{alphas[k]}\n')
			
			print(f'Calculating shot noise for {label_i}...')
			###############################################
			#load the N_g map and calculate the mean weighted by the mask
			mu_N = load_tomographic_maps(PATH_FD + cf.maps.ngal_maps, idx=i)[0][above_thresh].sum() / sum_w_above_thresh
			#calculate the noise power spectrum
			cl_noise_coupled = np.full(ell_max, Apix * mu_w / mu_N).reshape((1,ell_max))
			#decouple
			cl_noise_decoupled = w.decouple_cell(cl_noise_coupled)
		else:
			cl_noise_coupled = np.zeros((1, ell_max))
			cl_noise_decoupled = np.zeros((1, cf.nbpws))
		
		print(f'Calculating deprojection bias for pairing {label_i},{label_j}...')
		##############################################################
		if nsyst > 0 and not cf.lite:
			cl_bias = nmt.deprojection_bias(f_i, f_j, cl_guess)
		else:
			print('No systematics maps provided; skipping deprojection bias calculation.')
			cl_bias = np.zeros_like(cl_guess)
		
		#multiplicative correction to delta_g of (1 / (1-Fs)) due to stars results in factor of (1 / (1 - Fs))^2 correction to Cl
		if cf.correct_for_stars:
			mult = (1 / (1 - cf.Fs_fiducial)) ** 2.
			cl_coupled *= mult
			cl_guess *= mult
		
		print(f'Calculating decoupled C_ells for pairing {label_i},{label_j}...')
		###########################################################
		#compute the decoupled C_ell (w/o deprojection)
		cl_decoupled_nd = w.decouple_cell(cl_coupled_nd)
		#compute the decoupled, debiased C_ell (w/ deprojection)
		cl_decoupled = w.decouple_cell(cl_coupled, cl_bias=cl_bias)
		#decouple the bias C_ells as well
		cl_bias_decoupled = w.decouple_cell(cl_bias)

		#check if fields requied for future calculations
		if not any([label_i in p for p in label_pairs_now[idx:]]):
			f_i_nd[label_i] = None
			f_i[label_i] = None
		if not any([label_j in p for p in label_pairs_now[idx:]]):
			f_j_nd[label_j] = None
			f_j[label_j] = None

		#update the relevant arrays with the computed C_ells
		cls_now[idx,:,:] = cl_decoupled
		cls_nd_now[idx,:,:] = cl_decoupled_nd
		cls_noise_now[idx,:,:] = cl_noise_decoupled
		cls_bias_now[idx,:,:] = cl_bias_decoupled
	
	#delete systmaps to free up some memory
	del systmaps			
	#delete the NamasterFields to save some memory
	del f_i, f_j, f_i_nd, f_j_nd, fdone, fdone_nd


	#gather (decoupled) results
	comm.Gather(cls_now, cl_buff, root=0)
	comm.Gather(cls_nd_now, cl_nd_buff, root=0)
	comm.Gather(cls_noise_now, cl_noise_buff, root=0)
	comm.Gather(cls_bias_now, cl_bias_buff, root=0)

	if rank == 0:
		#remove all-zero rows from the Gathered arrays 
		cl_buff = cl_buff[~np.all(np.isnan(cl_buff), axis=3)].reshape(npairs, 1, cf.nbpws)
		cl_nd_buff = cl_nd_buff[~np.all(np.isnan(cl_nd_buff), axis=3)].reshape(npairs, 1, cf.nbpws)
		cl_noise_buff = cl_noise_buff[~np.all(np.isnan(cl_noise_buff), axis=3)].reshape(npairs, 1, cf.nbpws)
		cl_bias_buff = cl_bias_buff[~np.all(np.isnan(cl_bias_buff), axis=3)].reshape(npairs, 1, cf.nbpws)

		print('Calculating covariances (without deprojection)...')
		##########################################################

		#reload maps and deproject using saved information
		dg_maps = load_tomographic_maps(PATH_FD + cf.maps.deltag_maps)
		density_fields_nd = [nmt.NmtField(mask.mask_full, [dg], templates=None) for dg in dg_maps]
		#set up an array for the covariance matrices
		covar_all_nd = np.zeros((npairs, cf.nbpws, npairs, cf.nbpws))

		#cycle through the possible combinations of pairs of fields
		id_i = 0
		for i1 in range(cf.nsamples):
			for i2 in range(i1, cf.nsamples):
				id_j = 0
				for j1 in range(cf.nsamples):
					for j2 in range(j1, cf.nsamples):
						covar_all_nd[id_i, :, id_j, :] = cu.compute_covariance(w, cw,
													   density_fields_nd[i1],
													   density_fields_nd[i2],
													   density_fields_nd[j1],
													   density_fields_nd[j2]
													   )[0]
						id_j += 1
				id_i += 1

		#reshape the covariance matrix
		covar_all_nd = covar_all_nd.reshape((npairs * cf.nbpws, npairs * cf.nbpws))

		if nsyst > 0:
			print('Calculating covariances (with deprojection)...')
			##########################################################

			#retrieve the deprojection coefficients and make deprojected fields
			alphas_dfs = [pd.read_csv(PATH_CACHE + cf.cache_files.deproj.alphas[:-4] + f'_{k}.txt', sep='\t', index_col=0)
							for k in cf.samples]
			density_fields = [make_deprojected_field(dg, al) for dg, al in zip(dg_maps, alphas_dfs)]
			#set up an array for the covariance matrices
			covar_all = np.zeros((npairs, cf.nbpws, npairs, cf.nbpws))

			#cycle through the possible combinations of pairs of fields
			id_i = 0
			for i1 in range(cf.nsamples):
				for i2 in range(i1, cf.nsamples):
					id_j = 0
					for j1 in range(cf.nsamples):
						for j2 in range(j1, cf.nsamples):
							covar_all[id_i, :, id_j, :] = cu.compute_covariance(w, cw,
														density_fields[i1],
														density_fields[i2],
														density_fields[j1],
														density_fields[j2],
														f_sky=mu_w2
														)[0]
							id_j += 1
					id_i += 1

			#reshape the covariance matrix
			covar_all = covar_all.reshape((npairs * cf.nbpws, npairs * cf.nbpws))
		else:
			covar_all = covar_all_nd
		
	
		print('Constructing SACCs...')
		##############################
		#get the bandpower window functions
		wins = w.get_bandpower_windows()[0, :, 0, :].T
		wins = sacc.BandpowerWindow(np.arange(ell_max), wins)

		#set up the various SACC files
		s_main = sacc.Sacc()		# main results (i.e. w/ deprojection)
		s_nodeproj = sacc.Sacc()	# results w/o deprojection
		s_noise = sacc.Sacc()		# noise power spectra
		s_bias = sacc.Sacc()		# deprojection bias

		#get the n(z) distributions
		if cf.use_dir:
			nofz_info = cf.nofz_files.nz_dir
		else:
			nofz_info = f'{PATH_FD}{cf.nofz_files.nz_mc}'
		with h5py.File(nofz_info, 'r') as hf:
			#get the redshifts at which n(z) distributions are defined
			z = hf['z'][:]
			#add tracers to the Sacc object (one for each redshift bin)
			for samp in cf.samples:
				#get the n(z) distribution for this bin
				nz = hf[f'nz_{samp}'][:]
				s_main.add_tracer('NZ',	#n(z)-type tracer
							samp,	#tracer name
							quantity='galaxy_density', #quantity
							spin=0,
							z=z,
							nz=nz
							)
				s_nodeproj.add_tracer('NZ',	#n(z)-type tracer
							samp,	#tracer name
							quantity='galaxy_density', #quantity
							spin=0,
							z=z,
							nz=nz
							)
				s_noise.add_tracer('NZ',	#n(z)-type tracer
							samp,	#tracer name
							quantity='galaxy_density', #quantity
							spin=0,
							z=z,
							nz=nz
							)
				s_bias.add_tracer('NZ',	#n(z)-type tracer
							samp,	#tracer name
							quantity='galaxy_density', #quantity
							spin=0,
							z=z,
							nz=nz
							)

		#cycle through the bin pairings
		for ip, ((i,j), (label_i, label_j)) in enumerate(zip(pairings, label_pairs)):
			#add the relevant c_ell info to the Sacc
			s_main.add_ell_cl('cl_00',
						label_i, label_j,
						ell_effs,
						cl_buff[ip][0]-cl_noise_buff[ip][0],
						window=wins
						)
			s_nodeproj.add_ell_cl('cl_00',
						label_i, label_j,
						ell_effs,
						cl_nd_buff[ip][0]-cl_noise_buff[ip][0],
						window=wins
						)
			s_bias.add_ell_cl('cl_00',
							label_i, label_j,
							ell_effs,
							cl_bias_buff[ip][0],
							window=wins
							)
			s_noise.add_ell_cl('cl_00',
						label_i, label_j,
						ell_effs,
						cl_noise_buff[ip][0],
						window=wins)

		#add the covariance matrix to the Sacc
		s_main.add_covariance(covar_all)
		s_nodeproj.add_covariance(covar_all_nd)

		#save the SACC files
		s_main.save_fits(f'{PATH_FD}{cf.sacc_files.main}', overwrite=True)
		s_nodeproj.save_fits(f'{PATH_FD}{cf.sacc_files.nodeproj}', overwrite=True)
		s_noise.save_fits(f'{PATH_FD}{cf.sacc_files.noise}', overwrite=True)
		s_bias.save_fits(f'{PATH_FD}{cf.sacc_files.bias}', overwrite=True)
