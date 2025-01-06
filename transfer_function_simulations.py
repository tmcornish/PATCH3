#####################################################################################################
# - Uses simulated maps to compare deprojection bias vs transfer function as a
#   means of obtaining unbiased angular power spectra. This is a standalone script and
#   not intended as part of the pipeline, but may be worked into it eventually.
# - Produces up to 3 suites of simulations:
#	- DB suite: input C_ell is a best-fit HOD to the HSC data
# 	- TF suite 1: additional suite using the same input C_ell as suite 1
#	- TF suite 2: (optional) input C_ell = 1 / (ell + 10)
#####################################################################################################

import healpy as hp
import numpy as np
from scipy.optimize import minimize
import h5py
import pyccl as ccl
import pymaster as nmt
import matplotlib.pyplot as plt
import os
import config
import map_utils as mu
import output_utils as outils
import warnings
import glob
import pandas as pd


### SETTINGS & GLOBALS ###
cf  = config.fitHods

pairings = [(0,0)]
npairs = len(pairings)
#fiducial comsology
cosmo = ccl.Cosmology(**cf.cosmo_fiducial)

#range of wavenumbers and scale factors over which theory power spectra will be computed
lk_arr = np.log(np.geomspace(1E-4, 100, 256))
k_arr = np.e ** lk_arr
z_arr = np.linspace(0, 6, 16)[::-1]
a_arr = 1. / (1. + z_arr)

#create a halo model
m200def = ccl.halos.MassDef200m							#halo mass definition
hmf = ccl.halos.MassFuncTinker08(mass_def=m200def)		#halo mass function
hbf = ccl.halos.HaloBiasTinker10(mass_def=m200def)		#halo bias function
conc = ccl.halos.ConcentrationDuffy08(mass_def=m200def)	#halo concentration function
prf = ccl.halos.HaloProfileNFW(mass_def=m200def, concentration=conc)	#halo NFW profile
#feed this information to a HaloMassCalculator
hmc = ccl.halos.HMCalculator(
	mass_function=hmf,
	halo_bias=hbf,
	mass_def=m200def
	)
#2pt correlator
prof2pt = ccl.halos.Profile2ptHOD()

#correction using halofit power spectrum
pk_mm_halofit = ccl.power.nonlin_matter_power(cosmo, k_arr, a_arr)
pk_mm_halomodel = ccl.halos.halomod_power_spectrum(cosmo, hmc, k_arr, a_arr, prf)
R_kz = pk_mm_halofit / pk_mm_halomodel
#restrict the correction to the k-range where it is most relevant
R_kz[:, (k_arr < 0.02)+(k_arr >= 8)] = 1
#create Pk2D object containing the correction values
pk_R = ccl.pk2d.Pk2D(a_arr=a_arr, lk_arr=lk_arr, pk_arr=R_kz, is_logp=False)

#directory containing the measured HSC C_ells and maps
PATH_FD = cf.PATH_OUT + 'combined/'
#directory containing the systematics maps
PATH_SYST = PATH_FD + 'systmaps/'
#directory for the synthetic maps
PATH_SIMS = PATH_FD + 'sims/'
if not os.path.exists(PATH_SIMS):
	os.system(f'mkdir {PATH_SIMS}')

#number of maps to simulate (per suite)
nsim = 100
#map/C_ell parameters
ell_max = 3 * cf.nside_hi - 1
ells_theory = np.arange(0, ell_max+1)
npix = hp.nside2npix(cf.nside_hi)
bpw_edges = np.linspace(cf.ell_min, ell_max+1, cf.nbpws).astype(int)
b = nmt.NmtBin.from_edges(bpw_edges[:-1], bpw_edges[1:])
ell_effs = b.get_effective_ells()

#number of digits to use for string-formatted numbers when naming files
ndigit = int(np.floor(np.log10(nsim))) + 1

#whether to use a different C_ell as input for the transfer function simulations
diff_cl_tf = True
#C_ell to use if the above is true
cl_in_alt = 1 / (ells_theory + 10)

#############################
######### FUNCTIONS #########
#############################

def get_zeff(i):
	'''
	Computes the effective redshift for the specified tomographic bin.
	'''
	#compute comoving distance for each tomographic bin included in pairings
	z = tracers[i].z
	nz = tracers[i].nz
	#effective redshift of bin
	zeff = np.sum(z * nz) / np.sum(nz)
	return zeff

def scale_cuts():
	'''
	Computes a suitable scale cut to apply prior to fitting, depending on the effective
	redshifts of the tomographic bins.

	Parameters
	----------
	pairings: list[tuple]
		List of bin pairings involved in the fitting.
	
	fd: str
		Name of the field currently being analysed.
	
	Returns
	-------
	cuts: dict
		Dictionary of scale cuts to use for each bin pairing.
	'''
	#set up dictionary for the scale cuts
	ell_max_dict = {}
	#compute comoving distance for each tomographic bin included in pairings
	for i in np.unique(pairings):
		z = tracers[i].z
		nz = tracers[i].nz
		#effective redshift of bin
		zeff = np.sum(z * nz) / np.sum(nz)
		#comoving distance (Mpc)
		chi = cosmo.comoving_radial_distance(1. / (1. + zeff))
		#ell corresponding to kmax (Mpc^{-1}) at zeff
		ell_max = cf.kmax * chi
		ell_max_dict[i] = ell_max
	#now cycle through each bin pairing
	cuts = {}
	for p in pairings:
		i, j = p
		#take the minimum lmax for each bin pairing
		cuts[p] = min(ell_max_dict[i], ell_max_dict[j])
	return cuts


def apply_scale_cuts(ells, cells, cov, return_cuts=False, compute_cuts=True, hard_lmax=2000):
	'''
	Applies scale cuts to the data involved in the fit.

	Parameters
	----------
	cuts: dict
		Dictionary of scale cuts for each bin pairing.
	
	pairings: list[tuple]
		All bin pairings involved in the fit.
	
	ells: list[numpy.ndarray]
		List of arrays containing the effective multipoles for each 
		bin pairing.
	
	cells: list[numpy.ndarray]
		List of arrays containing the angular power spectra for each 
		bin pairing.
	
	cov: numpy.ndarray
		Covariance matrix for the angular power spectra.
	
	Returns
	-------

	'''
	#decide on scale cuts to use
	if compute_cuts:
		cuts = scale_cuts()
	else:
		cuts = {p : hard_lmax for p in pairings}

	#copy the inputs
	ells_cut = ells.copy()
	cells_cut = cells.copy()
	cov_cut = cov.copy()

	#set up a list for containing the masks for each bin pairing
	masks = []
	#cycle through the bin pairings
	for ip, p in enumerate(pairings):
		#get the maximum multipole allowed for the fit
		lmax = cuts[p]
		#mask the higher multipoles
		lmask = ells_cut[ip] <= lmax
		ells_cut[ip] = ells_cut[ip][lmask]
		cells_cut[ip] = cells_cut[ip][lmask]
		#append the mask to the list
		masks.append(lmask)
	#to mask the covariance matrix, combine and flatten all masks, then
	#take the outer product with itself
	masks = np.array(masks).flatten()
	nkeep = int(masks.sum())
	covmask = np.outer(masks, masks)
	cov_cut = cov[covmask].reshape((nkeep, nkeep))
	
	if return_cuts:
		return ells_cut, cells_cut, cov_cut, cuts
	else:
		return ells_cut, cells_cut, cov_cut


#######################################################
###############    START OF SCRIPT    #################
#######################################################

print('Creating input C_ell for simulations...')
################################################################################

#retrieve relevant data for the specified bin pairing
s = outils.select_from_sacc(PATH_FD + cf.outsacc,
							[(f'cl{i}', f'cl{j}') for i,j in pairings],
							'cl_00')
ell_cl = [s.get_ell_cl('cl_00', i, j) for i,j in s.get_tracer_combinations()]
ells_all = [ell_cl[i][0] for i in range(npairs)]
cells_all = [ell_cl[i][1] for i in range(npairs)]
cov_all = s.covariance.covmat

#construct NumberCountsTracer objects from the saved n(z) info
tracers = [s.tracers[i] for i in s.tracers.keys()]
NCT = [
	ccl.NumberCountsTracer(
		cosmo,
		has_rsd=False,
		dndz=(t.z, t.nz),
		bias=(t.z, np.ones_like(t.z))
	) for t in tracers
]

#apply scale cuts
ells, cells, cov = apply_scale_cuts(ells_all, cells_all, cov_all, compute_cuts=False, hard_lmax=2*cf.nside_hi)
err_cell = np.sqrt(np.diag(cov))
icov = np.linalg.inv(cov)

a_pivot = 1 / (1 + 0.65)
suppress_1h = lambda a: 0.01

def model_cell(theta, ells):
	#free parameters in the model
	logM0, logM1, mu0p, mu1p, alpha = theta
	smooth_transition = lambda a: alpha
	#halo profile
	prof = ccl.halos.HaloProfileHOD(
		mass_def=m200def,
		concentration=conc,
		log10Mmin_0=logM0,
		log10M0_0=logM0,
		log10M1_0=logM1,
		log10M0_p=mu0p,
		log10Mmin_p=mu0p,
		log10M1_p=mu1p,
		a_pivot=a_pivot
		)
	#halo-model power spectrum for galaxies
	pk = ccl.halos.halomod_Pk2D(
		cosmo,
		hmc,
		prof,
		prof_2pt=prof2pt,
		prof2=prof,
		a_arr=a_arr,
		lk_arr=lk_arr,
		smooth_transition=smooth_transition,
		suppress_1h=suppress_1h
	)
	#compute theory C_ells
	theory_cells = [
		ccl.angular_cl(
			cosmo,
			NCT[i],
			NCT[j],
			ells,
			p_of_k_a=pk
		)
		for i,j in pairings
		]
	
	return theory_cells


#best-fit parameters from HOD fitting (performed elsewhere)
hod_params = (
	11.15910346,	#mu_min
	12.97450235,	#mu_1
	-4.45496618,	#mu_{min,p}
	-4.91124055,	#mu_{1,p}
	0.38082766		#alpha_smooth
)

#compute theory C_ells
cl_in = model_cell(hod_params, ells=np.array(ells_theory, dtype=int))[0]

print('Loading HSC mask...')
############################
maskdata = mu.MaskData(PATH_FD + cf.survey_mask)
mask = maskdata.mask_full
mu_w2 = maskdata.meansq

print('Setting up workspaces...')
#################################
df = nmt.NmtField(mask, None, spin=0)

print('Workspaces: calculating mode coupling matrix...')
########################################################
w = nmt.NmtWorkspace()
wsp_file = PATH_SIMS + 'sims_workspace.fits'
if os.path.exists(wsp_file):
	w.read_from(wsp_file)
else:
	w.compute_coupling_matrix(df, df, b)
	w.write_to(wsp_file)

print('Workspaces: calculating coupling coefficients...')
########################################################
cw = nmt.NmtCovarianceWorkspace()
covwsp_file = PATH_SIMS + 'sims_covworkspace.fits'
if os.path.exists(covwsp_file):
	cw.read_from(covwsp_file)
else:
	cw.compute_coupling_coefficients(df, df)
	cw.write_to(covwsp_file)

#convolve the input C_ell(s) with the bandpower window functions
cl_in_decoupled = w.decouple_cell(w.couple_cell(cl_in.reshape(1,ell_max+1)))
if diff_cl_tf:
	cl_in_alt_decoupled = w.decouple_cell(w.couple_cell(cl_in_alt.reshape(1,ell_max+1)))

print('Loading systematics maps...')
####################################
if 'all' in map(str.lower, cf.systs):
	systs = glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp')
else:
	systs = [PATH_SYST + s for s in cf.systs]

#shuffle the systematics
np.random.seed(0)
np.random.shuffle(systs)

if cf.Nsyst_max:
	systs = systs[:cf.Nsyst_max]
nsyst = len(systs)

systmaps = np.array(
				[mu.load_map(s, is_systmap=True, mask=maskdata, apply_mask=True) for s in systs]
				).reshape([nsyst, 1, npix])


print(f'Loaded {nsyst} systematics.')

#retrieve example deprojection coefficients from relevant file
alphas_df = pd.read_csv(
						PATH_FD + 'cache/deprojection_alphas_1024_bin0.txt',
						sep='\t',
						index_col=0
						)
alphas_in = alphas_df.loc[[os.path.basename(s) for s in systs]]['alpha'].values
#multiply by 5 to exaggerate effect of each systematic
alphas_in *= 5
print('Input deprojection coefficients:')
print(alphas_in)


#list of required datasets to search for from existing simulaitons
out_required = [
	'map_in',
	'cl_meas',
	'cl_meas_nd',
	'cov',
	'cov_nd',
	'err_cell',
	'err_cell_nd',
	'cl_clean',
	'cl_clean_coupled',
	'cl_meas_coupled',
	'cl_meas_coupled_nd',
	'cl_guess',
	'cl_guess_nd',
	'cl_bias',
	'alphas_meas'	
]

#total number of simulations
if diff_cl_tf:
	nsuites = 3
else:
	nsuites = 2
nsim_tot = nsuites * nsim

#whether to compute deprojection bias
compute_db = True

print(f'Generating {nsim_tot}+1 synthetic maps...')
#############################################
for i in range(nsim_tot+1):
	#string form of iteration index
	i_str = str(i % nsim).zfill(ndigit)
	#define distinctive string for output filename
	if i == 0:
		id_str = f'test_sim'
		cl_in_now = cl_in
		cl_in_decoupled_now = cl_in_decoupled
	elif 0 < i <= nsim:
		id_str = f'sim{i_str}_DBsuite'
	else:
		#switch off deprojection bias calculation
		if i == nsim + 1:
			compute_db = False
			out_required.remove('cl_bias')
		if nsim < i <= (2*nsim):
			id_str = f'sim{i_str}_TFsuite1'
		else:
			if i == (nsim * 2) + 1:
				cl_in_now = cl_in_alt
			cl_in_decoupled_now = cl_in_alt_decoupled
			id_str = f'sim{i_str}_TFsuite2'
		
	#filename for outputs from this simulation
	outfile = f'{PATH_SIMS}{id_str}_nside{cf.nside_hi}.hdf5'

	#set up dictionary for storing outputs
	out_dict = {
		'ell_effs' : ell_effs,
		'cl_in' : cl_in_now,
		'cl_in_decoupled' : cl_in_decoupled_now,
		'systs' : systs,
		'alphas_in' : alphas_in
	}
	#if outputs exist for this iteration, skip
	if os.path.exists(outfile):
		print(f'Output file found for {id_str}; checking for required outputs...')
		with h5py.File(outfile, 'r') as hf:
			for k in out_required:
				if k in hf.keys():
					out_dict[k] = hf[k][:]
				else:
					out_dict[k] = None
	else:
		out_dict = {**out_dict, **{k : None for k in out_required}}

	if out_dict['map_in'] is None:
		print(f'{id_str}: synthesising, masking and contaminating map...')
		#synthesise the map from the input C_ells
		np.random.seed(i)
		out_dict['map_in'] = hp.synfast(cl_in_now, nside=cf.nside_hi) 

	#apply the survey mask
	masked_pix = (mask == 0.)
	map_masked = out_dict['map_in'] * mask

	#contaminate with systematics
	map_cont = map_masked + np.sum(
								alphas_in[:,None,None] * systmaps,
								axis=0
								)
	if out_dict['cl_clean_coupled'] is None:
		print(f'{id_str}: Computing C_ells prior to contamination...')
		#calculate (deprojected) angular power spectra
		df = nmt.NmtField(mask, [out_dict['map_in']], templates=None, masked_on_input=False)
		out_dict['cl_clean_coupled'] = nmt.compute_coupled_cell(df, df)
	if out_dict['cl_clean'] is None:
		out_dict['cl_clean'] = w.decouple_cell(out_dict['cl_clean_coupled'])

	if out_dict['cl_meas_coupled_nd'] is None:
		print(f'{id_str}: Computing C_ells without deprojection...')
		#calculate (deprojected) angular power spectra
		df = nmt.NmtField(mask, [map_cont], templates=None, masked_on_input=True)
		out_dict['cl_meas_coupled_nd'] = nmt.compute_coupled_cell(df, df)
	if out_dict['cl_meas_nd'] is None:
		out_dict['cl_meas_nd'] = w.decouple_cell(out_dict['cl_meas_coupled_nd'])
	
	if out_dict['cl_meas_coupled'] is None:
		print(f'{id_str}: Computing deprojected C_ells...')
		#calculate (deprojected) angular power spectra
		df = nmt.NmtField(mask, [map_cont], templates=systmaps, masked_on_input=True)
		out_dict['alphas_meas'] = df.alphas
		out_dict['cl_meas_coupled'] = nmt.compute_coupled_cell(df, df)
	
	if out_dict['cl_meas'] is None:
		out_dict['cl_meas'] = w.decouple_cell(out_dict['cl_meas_coupled'])

	if out_dict['cl_guess'] is None:
		#covariances and errors
		out_dict['cl_guess'] = out_dict['cl_meas_coupled'] / mu_w2
	
	if out_dict['cl_guess_nd'] is None:
		#covariances and errors
		out_dict['cl_guess_nd'] = out_dict['cl_meas_coupled_nd'] / mu_w2

	if out_dict['cov_nd'] is None:
		print(f'{id_str}: computing covariances (no deprojection)...')
		out_dict['cov_nd'] = nmt.gaussian_covariance(cw,
									0, 0, 0, 0,
									[out_dict['cl_guess_nd'][0]],
									[out_dict['cl_guess_nd'][0]],
									[out_dict['cl_guess_nd'][0]],
									[out_dict['cl_guess_nd'][0]],
									w)
	if out_dict['err_cell_nd'] is None:
		out_dict['err_cell_nd'] = np.sqrt(np.diag(out_dict['cov_nd']))

	if out_dict['cov'] is None:
		print(f'{id_str}: computing covariances...')
		out_dict['cov'] = nmt.gaussian_covariance(cw,
									0, 0, 0, 0,
									[out_dict['cl_guess'][0]],
									[out_dict['cl_guess'][0]],
									[out_dict['cl_guess'][0]],
									[out_dict['cl_guess'][0]],
									w)
	if out_dict['err_cell'] is None:
		out_dict['err_cell'] = np.sqrt(np.diag(out_dict['cov']))

	if compute_db:
		if out_dict['cl_bias'] is None:
			print(f'{id_str}: computing deprojection bias...')
			#deprojection bias, using measured C_ell as estimate for true
			cl_bias_coupled = nmt.deprojection_bias(df, df, out_dict['cl_guess'])
			out_dict['cl_bias'] = w.decouple_cell(cl_bias_coupled)


	print(f'{id_str}: Saving outputs...')
	with h5py.File(outfile, 'w') as hf:
		for k in out_dict.keys():
			hf.create_dataset(k, data=out_dict[k])
