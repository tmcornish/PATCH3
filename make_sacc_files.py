#####################################################################################################
# - Consolidates all information about the C_ells into SACC files.
#####################################################################################################

import config
import h5py
import sacc
import pymaster as nmt
import numpy as np
from output_utils import colour_string

### SETTINGS ###
cf = config.makeSaccFiles

#######################################################
###############    START OF SCRIPT    #################
#######################################################

#get the possible bin pairings
pairings, _ = cf.get_bin_pairings()
#bandpower edges and maximum ell at which the c_ells are defined
bpw_edges = cf.get_bpw_edges()
ell_max = bpw_edges[-1]
#get the effective ells
b = nmt.NmtBin.from_edges(bpw_edges[:-1], bpw_edges[1:])
ell_effs = b.get_effective_ells()


#cycle through the fields being analysed
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))

	#relevant directories
	PATH_INFO = f'{cf.PATH_OUT}{fd}/'
	PATH_CACHE = f'{PATH_INFO}cache/'

	#load the relevant workspaces from the cache
	w = nmt.NmtWorkspace()
	cw = nmt.NmtCovarianceWorkspace()
	w.read_from(f'{PATH_CACHE}{cf.wsp_file}')
	cw.read_from(f'{PATH_CACHE}{cf.covwsp_file}')
	#get the bandpower window functions
	wins = w.get_bandpower_windows()[0, :, 0, :].T
	wins = sacc.BandpowerWindow(np.arange(1, ell_max+1), wins)

	#set up the main SACC file
	s_main = sacc.Sacc()
	#get the n(z) distributions
	nofz_info = f'{PATH_INFO}{cf.nofz_file}'
	with h5py.File(nofz_info, 'r') as hf:
		#add tracers to the Sacc object (one for each redshift bin)
		for i in range(len(cf.zbins)-1):
			s_main.add_tracer('NZ',	#n(z)-type tracer
						f'gc_{i}',	#tracer name
						quantity='galaxy_density', #quantity
						spin=0,
						z=hf['z'][:],
						nz=hf[f'bin{i}'][:]
						)
	
	#cycle through the bin pairings
	for i,j in pairings:
		#retrieve the saved c_ell info
		cell_info = f'{cf.PATH_OUT}{fd}/{cf.outfile[:-5]}_{i}_{j}.hdf5'
		with h5py.File(cell_info, 'r') as hf:
			gp = hf[f'{i},{j}']
			#retrieve the relevant c_ell info
			cell_debiased = gp['cl_decoupled_debiased'][0]
			nell = gp['N_ell_decoupled'][0]
			cell_final = cell_debiased - nell
		#add the relevant c_ell info to the Sacc
		s_main.add_ell_cl('galaxy_density_cl',
					f'gc_{i}', f'gc_{j}',
					ell_effs,
					cell_final,
					window=wins
					)
	s_main.save_fits(f'{PATH_INFO}{cf.outsacc}', overwrite=True)