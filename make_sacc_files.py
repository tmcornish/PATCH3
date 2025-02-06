#####################################################################################################
# - Consolidates all information about the C_ells into SACC files.
#####################################################################################################

import sys
from configuration import PipelineConfig as PC
import h5py
import sacc
import pymaster as nmt
import numpy as np
from output_utils import colour_string
from cell_utils import get_bin_pairings, get_bpw_edges

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='makeSaccFiles')

#######################################################
###############    START OF SCRIPT    #################
#######################################################

#get the possible bin pairings
nbins = cf.nsamples
pairings, _ = get_bin_pairings(nbins)
ncross = len(pairings)
#bandpower edges and maximum ell at which the c_ells are defined
bpw_edges = get_bpw_edges(3*cf.nside_hi, ell_min=cf.ell_min, nbpws=cf.nbpws, spacing=cf.bpw_spacing)
ell_max = bpw_edges[-1]
#get the effective ells
b = nmt.NmtBin.from_edges(bpw_edges[:-1], bpw_edges[1:])
ell_effs = b.get_effective_ells()
n_ells = len(ell_effs)


#cycle through the fields being analysed
for fd in cf.fields:
	print(colour_string(fd.upper(), 'orange'))

	#relevant directories
	PATH_INFO = f'{cf.paths.out}{fd}/'
	PATH_CACHE = f'{PATH_INFO}cache/'

	#load the relevant workspaces from the cache
	w = nmt.NmtWorkspace()
	cw = nmt.NmtCovarianceWorkspace()
	w.read_from(f'{PATH_CACHE}{cf.cache_files.workspaces.wsp}')
	cw.read_from(f'{PATH_CACHE}{cf.cache_files.workspaces.covwsp}')
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
		nofz_info = cf.nofz_files.nz_dir
	else:
		nofz_info = f'{PATH_INFO}{cf.nofz_files.nz_mc}'
	with h5py.File(nofz_info, 'r') as hf:
		#get the redshifts at which n(z) distributions are defined
		z = hf['z'][:]
		#add tracers to the Sacc object (one for each redshift bin)
		for i, samp in enumerate(cf.samples):
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
	for i,j in pairings:
		#retrieve the saved c_ell info
		cell_info = f'{cf.paths.out}{fd}/{cf.cell_files.main[:-5]}_{i}_{j}.hdf5'
		with h5py.File(cell_info, 'r') as hf:
			gp = hf[f'{i},{j}']
			#retrieve the relevant c_ell info
			cell_nodeproj = gp['cl_decoupled_no_deproj'][0]
			cell_debiased = gp['cl_decoupled_debiased'][0]
			nell = gp['N_ell_decoupled'][0]
			bias = gp['cl_bias_decoupled'][0]
			cell_final = cell_debiased - nell
			cell_final_nodeproj = cell_nodeproj - nell
		#add the relevant c_ell info to the Sacc
		s_main.add_ell_cl('cl_00',
					f'cl{i}', f'cl{j}',
					ell_effs,
					cell_final,
					window=wins
					)
		s_nodeproj.add_ell_cl('cl_00',
						f'cl{i}', f'cl{j}',
						ell_effs,
						cell_final_nodeproj,
						window=wins
						)
		s_bias.add_ell_cl('cl_00',
						f'cl{i}', f'cl{j}',
						ell_effs,
						bias,
						window=wins
						)
		if i == j:
			s_noise.add_ell_cl('cl_00',
						f'cl{i}', f'cl{i}',
						ell_effs,
						nell,
						window=wins)
		
	#set up an array for the covariance matrices
	covar_all = np.zeros((ncross, n_ells, ncross, n_ells))
	covar_all_nodeproj = np.zeros((ncross, n_ells, ncross, n_ells))
	#open the file containing all the covariance matrices
	with h5py.File(f'{PATH_INFO}{cf.cell_files.covariances}', 'r') as hf:
		#cycle through the possible combinations of pairs of fields
		id_i = 0
		for i1 in range(nbins):
			for i2 in range(i1, nbins):
				id_j = 0
				for j1 in range(nbins):
					for j2 in range(j1, nbins):
						covar_all[id_i, :, id_j, :] = hf[f'{i1}{i2}-{j1}{j2}/final'][:]
						covar_all_nodeproj[id_i, :, id_j, :] = hf[f'{i1}{i2}-{j1}{j2}/no_deproj'][:]
						id_j += 1
				id_i += 1

	#reshape the covariance matrix
	covar_all = covar_all.reshape((ncross * n_ells, ncross * n_ells))
	covar_all_nodeproj = covar_all_nodeproj.reshape((ncross * n_ells, ncross * n_ells))
	#add the covariance matrix to the Sacc
	s_main.add_covariance(covar_all)
	s_nodeproj.add_covariance(covar_all_nodeproj)

	s_main.save_fits(f'{PATH_INFO}{cf.sacc_files.main}', overwrite=True)
	s_nodeproj.save_fits(f'{PATH_INFO}{cf.sacc_files.nodeproj}', overwrite=True)
	s_noise.save_fits(f'{PATH_INFO}{cf.sacc_files.noise}', overwrite=True)
	s_bias.save_fits(f'{PATH_INFO}{cf.sacc_files.bias}', overwrite=True)