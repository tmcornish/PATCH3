#####################################################################################################
# - Fits a halo occupation distribution model to the measured angular power spectra.
#####################################################################################################

import sys
from configuration import PipelineConfig as PC
from mpi4py import MPI
import numpy as np
import pyccl as ccl
import multiprocessing as mp
from scipy.optimize import minimize
from output_utils import colour_string
from cell_utils import select_from_sacc, get_bin_pairings, apply_scale_cuts
from cobaya.run import run
from cobaya.log import LoggedError

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='fitHods')
#fiducial comsology
cosmo = ccl.Cosmology(**cf.cosmo_fiducial)

#determine the bin pairings
pairings, _, label_pairs = get_bin_pairings(cf.nsamples, cf.auto_only, labels=list(cf.samples))
ncombos = len(pairings)

#range of wavenumbers and scale factors over which theory power spectra will be computed
lk_arr = np.log(np.geomspace(1E-4, 100, 256))
a_arr = 1. / (1. + np.linspace(0, 6, 16)[::-1])

#create a halo model
m200def = ccl.halos.MassDef200m							#halo mass definition
hmf = ccl.halos.MassFuncTinker08(mass_def=m200def)		#halo mass function
hbf = ccl.halos.HaloBiasTinker10(mass_def=m200def)		#halo bias function
conc = ccl.halos.ConcentrationDuffy08(mass_def=m200def)	#halo concentration function
#feed this information to a HaloMassCalculator
hmc = ccl.halos.HMCalculator(
	mass_function=hmf,
	halo_bias=hbf,
	mass_def=m200def
	)
#2pt correlator
prof2pt = ccl.halos.Profile2ptHOD()

#for redshift-dependent parameters
a_pivot = 1 / (1 + cf.z_pivot)

#############################
######### FUNCTIONS #########
#############################

def scale_cuts():
	'''
	Computes a suitable scale cut to apply prior to fitting, depending on the effective
	redshifts of the tomographic bins.
	
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
		#use this value only if it is less than 2*NSIDE
		ell_max_dict[i] = min(ell_max, 2*cf.nside_hi)
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
	ells: list[numpy.ndarray]
		List of arrays containing the effective multipoles for each 
		bin pairing.
	
	cells: list[numpy.ndarray]
		List of arrays containing the angular power spectra for each 
		bin pairing.
	
	cov: numpy.ndarray
		Covariance matrix for the angular power spectra.
	
	return_cuts: bool
		Returns the applied cuts if True.
	
	compute_cuts: bool
		If True, will compute the scale cuts using a value of k_max specified
		in the config file. If False, will use the value specified by the
		hard_lmax argument.
	
	hard_lmax: int
		Maximum multipole to use in compute_cuts is False.
		
	Returns
	-------
	ells_cut: list[numpy.ndarray]
		List of effective multipoles after applying the scale cuts.
	
	cells_cut: list[numpy.ndarray]
		List of C_ells after applying the scale cuts.
	
	cov_cut: numpy.ndarray
		Covariance matrix after applying the scale cuts.
	
	cuts: dict (optional)
		Dictionary containing the ell_max values used. 
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

	

def log_prior(theta):
	'''
	Defines the priors on the free parameters in the HOD model
	(in this case, mu_min and mu_1).

	Parameters
	----------
	theta: tuple
		Tuple containing values for the free parameters.
	
	Returns
	-------
	logp: float
		Either 0 or negative infinity, depending on whether the parameter
		values are within the bound sof the flat priors.
	'''
	mu_min, mu_1, mu_minp, mu_1p, alpha_smooth = theta
	if (0. < mu_min < 15.) and (0. < mu_1 < 17.) and (-10. < mu_minp < 10.) and (-12. < mu_1p < 15.) and (0.01 < alpha_smooth < 4):
		logp = 0.
	else:
		logp = -np.inf
	return logp


def suppress_1h(a):
	'''
	Returns the wavenumber below which the 1-halo term should be suppressed in
	the halo model power spectrum. Needs to be a function of the scale factor
	a in order to be compliant with CCL, but currently returns a constant.

	Parameters
	----------
	a: float
		Cosmic scale factor.
	
	Returns
	-------
	k_max: float
		Wavenumber (in Mpc^{-1}) below which the 1-halo term will be suppressed.
	'''
	k_max = 0.01
	return k_max


def log_likelihood(theta):
	'''
	Defines the priors on the free parameters in the HOD model.

	Parameters
	----------
	theta: tuple
		Tuple containing values for the free parameters.
	
	Returns
	-------
	logL: float
		Logarithm of the likelihood.
	'''
	#free parameters in the model
	mu_min, mu_1, mu_minp, mu_1p, alpha_smooth = theta
	smooth_transition = lambda a: alpha_smooth
	#halo profile
	prof = ccl.halos.HaloProfileHOD(
		mass_def=m200def,
		concentration=conc,
		log10Mmin_0=mu_min,
		log10M0_0=mu_min,
		log10M1_0=mu_1,
		log10M0_p=mu_minp,
		log10Mmin_p=mu_minp,
		log10M1_p=mu_1p,
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
		Tuple containing values for the free parameters.
	
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

#######################################################
###############    START OF SCRIPT    #################
#######################################################

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#cycle through the fields being analysed
for fd in cf.fields:
	
	PATH_FD = cf.paths.out + fd + '/'

	#load the Sacc file containing the power spectrum info
	s = select_from_sacc(PATH_FD + cf.sacc_files.main, label_pairs, 'cl_00')
	#get ells and cells (no scale cuts applied at this stage)
	ell_cl = [s.get_ell_cl('cl_00', i, j) for i, j in s.get_tracer_combinations()]
	ells_all = [ell_cl[i][0] for i in range(len(ell_cl))]
	cells_all = [ell_cl[i][1] for i in range(len(ell_cl))]
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
	ells, cells, cov, ell_cuts = apply_scale_cuts(ells_all, 
												cells_all, 
												cov_all, 
												return_cuts=True, 
												compute_cuts=cf.compute_scale_cuts 
												)
	#invert the covariance matrix
	icov = np.linalg.inv(cov)

	if rank == 0:
		print(colour_string(fd.upper(), 'orange'))
		#retrieve cobaya config options from config file
		info = cf.cobaya_info
		#add the likelihood
		info['likelihood'] = {'hod': log_likelihood}

		#perform initial fit to get initial positions of each walker
		initial = [p['ref'] for p in info['params']]
		ndim = len(initial)
		if cf.compute_initial:
			print('Estimating intial best fit...')
			initial = minimize(nll, initial).x
		#otherwise search for values in the config file
		else:
			print('Using initial guesses from config file...')
			
		print(f'Initial best-fit values:\n')
		for i,p in enumerate(info['params']):
			info['params'][p]['ref'] = initial.x[i]
			print(f'{p} = {initial.x[i]:.3f}')
		
		#ensure chains are saved to the correct place
		info['output'] = f'{cf.paths.out}{fd}/{info["output"]}{cf.suffix}'
	else:
		info = None
	info = comm.bcast(info, root=0)

	#run the sampler
	success = False
	try:
		updated_info, sampler = run(info)
		success = True
	except LoggedError as err:
		pass

	#did it work? (e.g. did not get stuck)
	success = all(comm.allgather(success))

	if not success and rank == 0:
		print('Sampling failed!')
	print('Sampling successful!')	
