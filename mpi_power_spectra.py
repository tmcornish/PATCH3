#####################################################################################################
# - Uses NaMaster to compute power spectra from the galaxy delta_g maps, deprojecting any 
#   systematics templates in the process.
# TODO: 
# - include possibility that number of nodes < number of bin pairings 
#####################################################################################################

import config
import healpy as hp
import healsparse as hsp
import numpy as np
from map_utils import load_map, load_tomographic_maps, MaskData
from cell_utils import get_bpw_edges, compute_covariance
import h5py
import pymaster as nmt
from output_utils import colour_string
import os
import sys
import glob
from mpi4py import MPI
import sacc

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



#######################################################
###############    START OF SCRIPT    #################
#######################################################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

#all possible pairings of tomographic bins (as tuples and as strings)
pairings, pairings_s = cf.get_bin_pairings()
npairs = len(pairings)

#maximum ell allowed by the resolution
ell_max = 3 * cf.nside_hi - 1
#get pixel area in units of steradians
Apix = hp.nside2pixarea(cf.nside_hi)
#get the number of pixels in a full-sky map at the required resolution
npix = hp.nside2npix(cf.nside_hi)

if cf.use_N19_bps:
	#retrieve bandpower edges from config
	bpw_edges = np.array(cf.bpw_edges).astype(int)
	#only include bandpowers < 3 * NSIDE (NOTE: NaMaster treats the upper edges as exclusive)
	bpw_edges = bpw_edges[bpw_edges <= ell_max+1]
else:
	bpw_edges = get_bpw_edges(cf.nside_hi, cf.nbpws, ell_min=cf.ell_min, log_spacing=cf.log_spacing)
#create pymaster NmtBin object using these bandpower objects
b = nmt.NmtBin.from_edges(bpw_edges[:-1], bpw_edges[1:])
#get the effective ells
ell_effs = b.get_effective_ells()


#cycle through the fields being analysed
for fd in cf.get_global_fields():
	#path to the directory containing the maps
	PATH_FD = f'{cf.PATH_OUT}{fd}/'
	#path to directory of cached outputs from this script
	PATH_CACHE = PATH_FD + 'cache/'
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
		mask = MaskData(PATH_FD + cf.survey_mask)
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
			dg_map, = load_tomographic_maps(PATH_FD + cf.deltag_maps, idx=rank)
			print(f'Creating NmtField for bin {rank} (without deprojection)...')
			####################################################################
			density_fields_nd = nmt.NmtField(mask.mask_full, [dg_map], templates=None, lite=cf.lite)
			
			print(f'Creating NmtField for bin {rank} (with deprojection)...')
			################################################################
			density_fields = nmt.NmtField(mask.mask_full, [dg_map], templates=systmaps, lite=cf.lite)

		density_fields_nd = comm.gather(density_fields_nd, root=0)
		density_fields = comm.gather(density_fields, root=0)

	else:
		if rank == 0:
			dg_map = load_tomographic_maps(PATH_FD + cf.deltag_maps)
			print(f'Creating NmtFields for all bins (without deprojection)...')
			###################################################################
			density_fields_nd = [nmt.NmtField(mask.mask_full, [dg], templates=None, lite=cf.lite)
									for dg in dg_map]
			
			print(f'Creating NmtField for all bins (with deprojection)...')
			###############################################################
			density_fields = [nmt.NmtField(mask.mask_full, [dg], templates=systmaps, lite=cf.lite)
					 				for dg in dg_map]
	

	#broadcast the list of NmtFields
	if rank != 0:
		density_fields_nd = None
		density_fields = None
	density_fields_nd = comm.bcast(density_fields_nd, root=0)
	density_fields = comm.bcast(density_fields, root=0)


	#set up buffers for gathering results
	cl_buff = None				#buffer for the main C_ells
	cl_nd_buff = None		#buffer for the C_ells without deprojection
	cl_noise_buff = None 	#buffer for the noise power spectra
	cl_bias_buff = None		#buffer for the deprojection bias
	if rank == 0:
		cl_buff = np.empty((npairs, 1, cf.nbpws), dtype=np.float64)			
		cl_nd_buff = np.empty((npairs, 1, cf.nbpws), dtype=np.float64)					
		cl_noise_buff = np.empty((npairs, 1, cf.nbpws), dtype=np.float64)			 	
		cl_bias_buff = np.empty((npairs, 1, cf.nbpws), dtype=np.float64)			

	#current bin pairing
	i, j = pairings[rank]
	#current density fields (no deprojection)
	f_i_nd = density_fields_nd[i]
	f_j_nd = density_fields_nd[j]
	#current density fields (with deprojection)
	f_i = density_fields[i]
	f_j = density_fields[j]

	print(f'Calculating coupled C_ells for pairing {i},{j}...')
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
	#multiplicative correction to delta_g of (1 / (1-Fs)) due to stars results in factor of (1 / (1 - Fs))^2 correction to Cl
	if cf.correct_for_stars:
		mult = (1 / (1 - cf.Fs_fiducial)) ** 2.
		cl_coupled *= mult
		cl_guess *= mult


	print(f'Calculating deprojection bias for pairing {i},{j}...')
	##############################################################
	if nsyst > 0 and not cf.lite:
		cl_bias = nmt.deprojection_bias(f_i, f_j, cl_guess)
	else:
		print('No systematics maps provided; skipping deprojection bias calculation.')
		cl_bias = np.zeros_like(cl_guess)

	if i == j:
		print(f'Saving deprojection coefficients for bin {i}...')
		##############################################################
		alphas = f_i.alphas
		with open(PATH_CACHE + cf.alphas_file[:-4] + f'_bin{i}.txt', 'w') as alphas_file:
			alphas_file.write('Sytematic\talpha\n')
			for k in range(nsyst):
				alphas_file.write(f'{systs[k]}\t{alphas[k]}\n')

	print(f'Calculating decoupled C_ells for pairing {i},{j}...')
	###########################################################
	#compute the decoupled C_ell (w/o deprojection)
	cl_decoupled_nd = w.decouple_cell(cl_coupled_nd)
	#compute the decoupled, debiased C_ell (w/ deprojection)
	cl_decoupled = w.decouple_cell(cl_coupled, cl_bias=cl_bias)
	#decouple the bias C_ells as well
	cl_bias_decoupled = w.decouple_cell(cl_bias)

	if i == j:
		print(f'Calculating shot noise for bin {i}...')
		###############################################
		#load the N_g map and calculate the mean weighted by the mask
		mu_N = load_tomographic_maps(PATH_FD + cf.ngal_maps, idx=i)[0][above_thresh].sum() / sum_w_above_thresh
		#calculate the noise power spectrum
		cl_noise_coupled = np.full(ell_max+1, Apix * mu_w / mu_N).reshape((1,ell_max+1))
		#decouple
		cl_noise_decoupled = w.decouple_cell(cl_noise_coupled)
	else:
		cl_noise_coupled = np.zeros((1, ell_max+1))
		cl_noise_decoupled = np.zeros((1, cf.nbpws))
	
	#gather (decoupled) results
	comm.Gather(cl_decoupled, cl_buff, root=0)
	comm.Gather(cl_decoupled_nd, cl_nd_buff, root=0)
	comm.Gather(cl_noise_decoupled, cl_noise_buff, root=0)
	comm.Gather(cl_bias_decoupled, cl_bias_buff, root=0)

	if rank == 0:
		print('Calculating covariances...')
		###################################
		#set up an array for the covariance matrices
		covar_all = np.zeros((npairs, cf.nbpws, npairs, cf.nbpws))
		covar_all_nd = np.zeros((npairs, cf.nbpws, npairs, cf.nbpws))

		#cycle through the possible combinations of pairs of fields
		id_i = 0
		for i1 in range(cf.nbins):
			for i2 in range(i1, cf.nbins):
				id_j = 0
				for j1 in range(cf.nbins):
					for j2 in range(j1, cf.nbins):
						covar_all[id_i, :, id_j, :] = compute_covariance(w, cw,
													   density_fields[i1],
													   density_fields[i2],
													   density_fields[j1],
													   density_fields[j2]
													   )[0]
						covar_all_nd[id_i, :, id_j, :] = compute_covariance(w, cw,
													   density_fields_nd[i1],
													   density_fields_nd[i2],
													   density_fields_nd[j1],
													   density_fields_nd[j2]
													   )[0]
						id_j += 1
				id_i += 1

		#reshape the covariance matrix
		covar_all = covar_all.reshape((npairs * cf.nbpws, npairs * cf.nbpws))
		covar_all_nodeproj = covar_all_nodeproj.reshape((npairs * cf.nbpws, npairs * cf.nbpws))
	
		print('Constructing SACCs...')
		##############################
		#get the bandpower window functions
		wins = w.get_bandpower_windows()[0, :, 0, :].T
		wins = sacc.BandpowerWindow(np.arange(1, ell_max+1), wins)

		#set up the various SACC files
		s_main = sacc.Sacc()		# main results (i.e. w/ deprojection)
		s_nodeproj = sacc.Sacc()	# results w/o deprojection
		s_noise = sacc.Sacc()		# noise power spectra
		s_bias = sacc.Sacc()		# deprojection bias

		#get the n(z) distributions
		if cf.use_dir:
			nofz_info = cf.nz_dir_file
		else:
			nofz_info = f'{PATH_FD}{cf.nz_mc_file}'
		with h5py.File(nofz_info, 'r') as hf:
			#get the redshifts at which n(z) distributions are defined
			z = hf['z'][:]
			#add tracers to the Sacc object (one for each redshift bin)
			for i in range(len(cf.zbins)-1):
				#get the n(z) distribution for this bin
				nz = hf[f'nz_{i}'][:]
				s_main.add_tracer('NZ',	#n(z)-type tracer
							f'cl{i}',	#tracer name
							quantity='galaxy_density', #quantity
							spin=0,
							z=z,
							nz=nz
							)
				s_nodeproj.add_tracer('NZ',	#n(z)-type tracer
							f'cl{i}',	#tracer name
							quantity='galaxy_density', #quantity
							spin=0,
							z=z,
							nz=nz
							)
				s_noise.add_tracer('NZ',	#n(z)-type tracer
							f'cl{i}',	#tracer name
							quantity='galaxy_density', #quantity
							spin=0,
							z=z,
							nz=nz
							)
				s_bias.add_tracer('NZ',	#n(z)-type tracer
							f'cl{i}',	#tracer name
							quantity='galaxy_density', #quantity
							spin=0,
							z=z,
							nz=nz
							)

		#cycle through the bin pairings
		for ip, (i,j) in enumerate(pairings):
			#add the relevant c_ell info to the Sacc
			s_main.add_ell_cl('cl_00',
						f'cl{i}', f'cl{j}',
						ell_effs,
						cl_buff[ip],
						window=wins
						)
			s_nodeproj.add_ell_cl('cl_00',
						f'cl{i}', f'cl{j}',
						ell_effs,
						cl_nd_buff[ip],
						window=wins
						)
			s_bias.add_ell_cl('cl_00',
							f'cl{i}', f'cl{j}',
							ell_effs,
							cl_bias_buff[ip],
							window=wins
							)
			s_noise.add_ell_cl('cl_00',
						f'cl{i}', f'cl{i}',
						ell_effs,
						cl_noise_buff[ip],
						window=wins)

		#add the covariance matrix to the Sacc
		s_main.add_covariance(covar_all)
		s_nodeproj.add_covariance(covar_all_nodeproj)

		#save the SACC files
		s_main.save_fits(f'{PATH_FD}{cf.outsacc}', overwrite=True)
		s_nodeproj.save_fits(f'{PATH_FD}{cf.outsacc[:-5]}_nodeproj.fits', overwrite=True)
		s_noise.save_fits(f'{PATH_FD}{cf.outsacc[:-5]}_noise.fits', overwrite=True)
		s_bias.save_fits(f'{PATH_FD}{cf.outsacc[:-5]}_deprojbias.fits', overwrite=True)
