#####################################################################################################
# - Uses simulated maps to compare deprojection bias vs transfer function as a
#   means of obtaining unbiased angular power spectra. This is a standalone script and
#   not intended as part of the pipeline, but may be worked into it eventually.
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

#number of maps to simulate
nsim = 100
#map/C_ell parameters
ell_max = 3 * cf.nside_hi - 1
ells_theory = range(0, ell_max+1)
npix = hp.nside2npix(cf.nside_hi)
bpw_edges = np.linspace(cf.ell_min, ell_max+1, cf.nbpws).astype(int)
b = nmt.NmtBin.from_edges(bpw_edges[:-1], bpw_edges[1:])
ell_effs = b.get_effective_ells()

#number of digits to use for string-formatted numbers when naming files
ndigit = int(np.floor(np.log10(nsim))) + 1

#############################
######### FUNCTIONS #########
#############################

def get_zeff(i):
	'''
	Computes the effective redshift for the specified tomographic bin.
	'''
	#file containing DIR-based n(z) distributions
	if cf.use_dir:
		nofz_file = cf.nz_dir_file
	else:
		nofz_file = PATH_FD + cf.nz_mc_file
	#open the n(z) file
	with h5py.File(nofz_file, 'r') as hf:
		#compute comoving distance for each tomographic bin included in pairings
		z = hf['z'][:]
		nz = hf[f'nz_{i}'][:]
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
	#file containing DIR-based n(z) distributions
	if cf.use_dir:
		nofz_file = cf.nz_dir_file
	else:
		nofz_file = PATH_FD + cf.nz_mc_file
	#open the n(z) file
	with h5py.File(nofz_file, 'r') as hf:
		#compute comoving distance for each tomographic bin included in pairings
		for i in np.unique(pairings):
			z = hf['z'][:]
			nz = hf[f'nz_{i}'][:]
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

def apply_scale_cuts(ells, cells, cov):
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
	if cf.compute_scale_cuts:
		cuts = scale_cuts()
	else:
		cuts = {p : cf.hard_lmax for p in pairings}

	#set up a list for containing the masks for each bin pairing
	masks = []
	#cycle through the bin pairings
	for ip, p in enumerate(pairings):
		#get the maximum multipole allowed for the fit
		lmax = cuts[p]
		#mask the higher multipoles
		lmask = ells[ip] <= lmax
		ells[ip] = ells[ip][lmask]
		cells[ip] = cells[ip][lmask]
		#append the mask to the list
		masks.append(lmask)
	#to mask the covariance matrix, combine and flatten all masks, then
	#take the outer product with itself
	masks = np.array(masks).flatten()
	nkeep = int(masks.sum())
	covmask = np.outer(masks, masks)
	cov = cov[covmask].reshape((nkeep, nkeep))
	
	return ells, cells, cov


#######################################################
###############    START OF SCRIPT    #################
#######################################################

print('Using observed HSC C_ells data to derive input C_ell for simulations...')
################################################################################

#retrieve relevant data for the specified bin pairing
s = outils.select_from_sacc(PATH_FD + cf.outsacc,
							[(f'cl{i}', f'cl{j}') for i,j in pairings],
							'cl_00')
ell_cl = [s.get_ell_cl('cl_00', i, j) for i,j in s.get_tracer_combinations()]
ells_all = [ell_cl[i][0] for i in range(npairs)]
cells_all = [ell_cl[i][1] for i in range(npairs)]
cov_all = s.covariance.covmat
#apply scale cuts
ells, cells, cov = apply_scale_cuts(ells_all, cells_all, cov_all)
err_cell = np.sqrt(np.diag(cov))
icov = np.linalg.inv(cov)

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

def log_prior(theta):
	'''
	Defines the priors on the free parameters in the HOD model
	(in this case, logM0 and logM1).

	Parameters
	----------
	theta: tuple
		Tuple containing values for the free parameters logM0 and logM1.
	
	Returns
	-------
	logp: float
		Either 0 or negative infinity, depending on whether the parameter
		values are within the bound sof the flat priors.
	'''
	logM0, logM1 = theta
	if (0. < logM0 < 15.) and (0. < logM1 < 17.):
		logp = 0.
	else:
		logp = -np.inf
	return logp

def log_likelihood(theta):
	'''
	Defines the priors on the free parameters in the HOD model
	(in this case, logM0 and logM1).

	Parameters
	----------
	theta: tuple
		Tuple containing values for the free parameters logM0 and logM1.
	
	Returns
	-------
	logL: float
		Logarithm of the likelihood.
	'''
	#free parameters in the model
	logM0, logM1 = theta
	#halo profile
	prof = ccl.halos.HaloProfileHOD(
		mass_def=m200def,
		concentration=conc,
		log10Mmin_0=logM0,
		log10M0_0=logM0,
		log10M1_0=logM1
		)
	#halo-model power spectrum for galaxies
	pk = ccl.halos.halomod_Pk2D(
		cosmo,
		hmc,
		prof,
		prof_2pt=prof2pt,
		prof2=prof,
		a_arr=a_arr,
		lk_arr=lk_arr
	)
	#correct using halofit
	pk = pk * pk_R
	#compute theory C_ells
	theory_cells = [
		ccl.angular_cl(
			cosmo,
			NCT[i],
			NCT[j],
			ells[ip],
			p_of_k_a=pk
		)
		for ip, (i,j) in enumerate(pairings)
		]
		
	#residuals
	diff = np.concatenate(cells) - np.concatenate(theory_cells)
	
	logL = -np.dot(diff, np.dot(icov, diff)) / 2.
	
	return logL

def log_probability(theta):
	'''
	Calculate the logged probability, factoring in the priors
	and the likelihood.

	Parameters
	----------
	theta: tuple
		Tuple containing values for the free parameters logM0 and logM1.
	
	Returns
	-------
	log_prob: float
		Logarithm of the total probability.
	'''
	#compute the logged prior
	lp = log_prior(theta)
	if not np.isfinite(lp):
		log_prob = -np.inf
	else:
		log_prob = lp + log_likelihood(theta)
	return log_prob

def nll(*args):
	'''
	Converts log_probability into a function that can be minimised
	to find the best-fit parameters.
	'''
	return -log_probability(*args)

#fit HOD
initial = [12, 12]
ndim = len(initial)
with warnings.catch_warnings():
	warnings.simplefilter('ignore')
	soln = minimize(nll, initial)
logM0, logM1 = soln.x

print('Results of HOD fit:')
print(soln)

#use the best-fit to generate C_ell^in (just use the first bin pairing since
#maps can only be generated for one C_ell at a time)

#halo profile
prof = ccl.halos.HaloProfileHOD(
	mass_def=m200def,
	concentration=conc,
	log10Mmin_0=logM0,
	log10M0_0=logM0,
	log10M1_0=logM1
	)
#halo-model power spectrum for galaxies
pk = ccl.halos.halomod_Pk2D(
	cosmo,
	hmc,
	prof,
	prof_2pt=prof2pt,
	prof2=prof,
	a_arr=a_arr,
	lk_arr=lk_arr
) * pk_R
#compute theory C_ells
cl_in = ccl.angular_cl(
		cosmo,
		NCT[pairings[0][0]],
		NCT[pairings[0][1]],
		np.array(ells_theory, dtype=int),
		p_of_k_a=pk
	)


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
				[mu.load_map(s, is_systmap=True, mask=maskdata, apply_mask=False) for s in systs]
				).reshape([nsyst, 1, npix])


print(f'Loaded {nsyst} systematics.')

#generate random deprojection coefficients for each systematic
np.random.seed(nsim+1)
alphas_in = np.random.randn(nsyst) * 1.5 - 1 
print('Input deprojection coefficients:')
print(alphas_in)

#lists for storing results from each simulation
cl_best_db = []		#cl_best from deprojection bias
cl_ratio = []		#ratio of measured cl to input cl

#list of required datasets to search for from existing simulaitons
out_required = [
	'map_in',
	'cl_meas',
	'cov',
	'err_cell',
	'cl_meas_coupled',
	'cl_guess',
	'cl_bias',
	'alphas_meas'	
]

print(f'Generating {nsim} synthetic maps...')
#############################################
for i in range(1,nsim+1):
	#string form of iteration index
	i_str = str(i).zfill(ndigit)
	#filename for outputs from this simulation
	outfile = f'{PATH_SIMS}sim{i_str}_nside{cf.nside_hi}.hdf5'

	#set up dictionary for storing outputs
	out_dict = {
		'ell_effs' : ell_effs,
		'cl_in' : cl_in,
		'systs' : systs,
		'alphas_in' : alphas_in
	}
	#if outputs exist for this iteration, skip
	if os.path.exists(outfile):
		print(f'Output file found for iteration {i}; checking for required outputs...')
		with h5py.File(outfile, 'r') as hf:
			for k in out_required:
				if k in hf.keys():
					out_dict[k] = hf[k][:]
				else:
					out_dict[k] = None

	if out_dict['map_in'] is None:
		print(f'Synthesising, masking and contaminating map {i}...')
		#synthesise the map from the input C_ells
		np.random.seed(i)
		out_dict['map_in'] = hp.synfast(cl_in, nside=cf.nside_hi) 

	#apply the survey mask
	masked_pix = (mask == 0.)
	map_masked = out_dict['map_in'] * mask

	#contaminate with systematics
	map_cont = map_masked + np.sum(
								alphas_in[:,None,None] * mask[None, None, :] * systmaps,
								axis=0
								)

	if out_dict['cl_meas_coupled'] is None:
		print(f'Sim {i}: Computing deprojected C_ells...')
		#calculate (deprojected) angular power spectra
		df = nmt.NmtField(mask, [map_cont], templates=systmaps)
		out_dict['alphas_meas'] = df.alphas
		out_dict['cl_meas_coupled'] = nmt.compute_coupled_cell(df, df)
	
	if out_dict['cl_meas'] is None:
		out_dict['cl_meas'] = w.decouple_cell(out_dict['cl_meas_coupled'])

	if out_dict['cl_guess'] is None:
		#covariances and errors
		cl_guess = out_dict['cl_meas_coupled'] / mu_w2

	if out_dict['cov'] is None:
		print(f'Sim {i}: computing covariances...')
		out_dict['cov'] = nmt.gaussian_covariance(cw,
									0, 0, 0, 0,
									[cl_guess[0]],
									[cl_guess[0]],
									[cl_guess[0]],
									[cl_guess[0]],
									w)
	if out_dict['err_cell'] is None:
		out_dict['err_cell'] = np.sqrt(np.diag(out_dict['cov']))

	if out_dict['cl_bias'] is None:
		print(f'Sim {i}: computing deprojection bias...')
		#deprojection bias, using measured C_ell as estimate for true
		cl_bias_coupled = nmt.deprojection_bias(df, df, out_dict['cl_guess'])
		out_dict['cl_bias'] = w.decouple_cell(cl_bias_coupled)

	#best estimate of unbiased C_ell
	cl_best = out_dict['cl_meas'] - out_dict['cl_bias']
	cl_best_db.append(cl_best)

	print(f'Sim {i}: calculating ratio of measured to input C_ell...')
	r_cl = w.decouple_cell(out_dict['cl_meas_coupled'] / cl_in)
	cl_ratio.append(r_cl)

	print(f'Sim {i}: Saving outputs...')
	with h5py.File(outfile, 'w') as hf:
		for k in out_dict.keys():
			hf.create_dataset(k, data=out_dict[k])
