#####################################################################################################
# - Fits a halo occupation distribution model to the measured angular power spectra.
#####################################################################################################

import config
import emcee
import sacc
import numpy as np
import itertools
import pyccl as ccl
import h5py
import multiprocessing as mp
from scipy.optimize import minimize
from output_utils import colour_string

### SETTINGS ###
cf = config.fitHods
#fiducial comsology
cosmo = ccl.Cosmology(**cf.cosmo_fiducial)

#range of wavenumbers and scale factors over which theory power spectra will be computed
lk_arr = np.log(np.geomspace(1E-4, 100, 256))
a_arr = 1. / (1. + np.linspace(0, 6, 100)[::-1])

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


#############################
######### FUNCTIONS #########
#############################

def get_data(s):
	'''
	Retrieves the relevant covariance informaiton, which depends on the
	user settings in config.py.

	Parameters
	----------
	s: sacc.sacc.Sacc
		Sacc object containing the angular power spectrum information.
	
	Returns
	-------
	ells: np.ndarray
		Array of multipoles at which the power spectra are defined (i.e. effective
		multipoles of each bandpower).
	
	cells: np.ndarray
		1D array containing the concatenation of all desired power spectra.
	
	icov: np.ndarray
		Inverted matrix of the relevant convariances.
	'''

	#get the full covariance matrix from the Sacc
	cov_full = s.covariance.covmat

	if cf.auto_only:
		#get data for the auto-correlations only
		ell_cl = np.array([s.get_ell_cl('cl_00', f'cl{i}', f'cl{i}') for i in range(cf.nbins)])
		#number of bandpowers (get from first ell-C_ell pair)
		nbpws = len(ell_cl[0,0])
		#length of the covariance matrix along each side
		L_cov = nbpws * cf.nbins
		cov = np.zeros((L_cov, L_cov))
		#map indices in the full covariance matrix to indices in the covariance matrix for the
		#auto power spectra
		i_auto = np.array(range(cf.nbins))
		i_full = np.array([np.sum(cf.nbins - i_auto[:j]) for j in i_auto])
		ij_auto = list(itertools.product(i_auto, i_auto))
		ij_full = list(itertools.product(i_full, i_full))
		#use this to retrieve the covariance matrix information for the auto-power spectra only
		for k in range(len(ij_auto)):
			i_old, j_old = ij_full[k]
			i_new, j_new = ij_auto[k]

			cov[i_new*nbpws:(i_new+1)*nbpws, j_new*nbpws:(j_new+1)*nbpws] = \
				cov_full[i_old*nbpws:(i_old+1)*nbpws, j_old*nbpws:(j_old+1)*nbpws]
		#possible bin pairings
		pairings = [(i,i) for i in range(cf.nbins)]
	
	else:
		#get the possible bin pairings from the config
		pairings, _ = cf.get_bin_pairings()
		#get data for all power spectra
		ell_cl = np.array([s.get_ell_cl('cl_00', i, j) for i, j in s.get_tracer_combinations()])
		cov = cov_full
	
	#same multipoles are used for all power spectra; retrieve these only
	ells = ell_cl[0,0]
	#flatten the C_ells into one array
	cells = ell_cl[:,1,:].flatten()
	#invert the covariance matrix
	icov = np.linalg.inv(cov)

	#return the ells, C_ells and inverted covariance matrix, along with the bin pairings
	return ells, cells, icov, pairings


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
	#compute theory C_ells
	theory_cells = np.array(
		[
		ccl.angular_cl(
			cosmo,
			NCT[i],
			NCT[j],
			ells,
			p_of_k_a=pk
		)
		for i,j in pairings
		]
		).flatten()
	#residuals
	diff = cells - theory_cells
	
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

#######################################################
###############    START OF SCRIPT    #################
#######################################################

if __name__ == '__main__':
	#cycle through the fields being analysed
	for fd in cf.get_global_fields():
		print(colour_string(fd.upper(), 'orange'))
		PATH_FD = cf.PATH_OUT + fd + '/'

		#load the Sacc file containing the power spectrum info
		s = sacc.Sacc.load_fits(PATH_FD + cf.outsacc)

		#get the relevant data
		ells, cells, icov, pairings = get_data(s)

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

		
		print('Estimating intial best fit...')
		######################################
		initial = [10, 10]
		ndim = len(initial)
		#soln = minimize(nll, initial)
		class soln():
			x = np.array([12.5, 14.2])

		print('Setting up the MCMC...')
		###############################
		#initial positions for the walkers
		p0 = [soln.x + np.array([cf.dlogM0 * np.random.rand(),
								cf.dlogM1 * np.random.rand()])
								for _ in range(cf.nwalkers)]
		p0 = np.array(p0)
		p0 -= np.array([cf.dlogM0/2., cf.dlogM1/2.])

		#set up an HDF5 backend
		backend = emcee.backends.HDFBackend(PATH_FD + cf.backend_file)
		backend.reset(cf.nwalkers, ndim)
		
		################################
		with mp.get_context('fork').Pool(20) as pool:
			#initialise the sampler
			sampler = emcee.EnsembleSampler(40, ndim, log_probability, pool=pool, backend=backend)
			#run the sampler
			print('Running sampler...')
			old_tau = np.inf
			for sample in sampler.sample(p0, iterations=cf.niter_max, progress=True):
				#check convergence time every 20 steps
				if sampler.iteration % 20:
					continue
				tau = sampler.get_autocorr_time(tol=0)
				#"convergence" reached if: 
				# - autocorrelation time < 1/100th the current iteration
				# - fractional change in autocorrelation time < 0.01
				converged = np.all(tau * 100 < sampler.iteration) \
							& np.all(np.abs(old_tau - tau) / tau < 0.01)
				if converged:
					break
				old_tau = tau

			#raise alert if convergence not reached
			if sampler.iteration == cf.niter_max - 1:
				print('WARNING: chains have not converged. Best-fit values may be inaccurate.')

			#print best-fit values
			theta0 = sampler.flatchain[np.argmax(sampler.flatlnprobability)]
			print(f'''Best-fit values:
				logM0: {theta0[0]:.3f}
				logM1: {theta0[1]:.3f}'''
				)
