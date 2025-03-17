###############################################################################
# - Fits a halo occupation distribution model to the measured angular power
#   spectra.
###############################################################################

import sys
from configuration import PipelineConfig as PC
from mpi4py import MPI
import numpy as np
import pyccl as ccl
from scipy.optimize import minimize
from output_utils import colour_string
from cell_utils import select_from_sacc, get_bin_pairings
from cobaya.run import run
from cobaya.log import LoggedError
import cell_utils as cu

# SETTINGS #
config_file = sys.argv[1]
cf = PC(config_file, stage='fitHods')
# Fiducial comsology
cosmo = ccl.Cosmology(**cf.cosmo_fiducial)

# Determine the bin pairings
pairings, _, label_pairs = get_bin_pairings(cf.nsamples,
                                            cf.auto_only,
                                            labels=list(cf.samples))
ncombos = len(pairings)

# Range of wavenumbers and scale factors over which theory power spectra will
# be computed
lk_arr = np.log(np.geomspace(1E-4, 100, 256))
a_arr = 1. / (1. + np.linspace(0, 6, 16)[::-1])

# Create a halo model
m200def = ccl.halos.MassDef200m							 # Halo mass definition
hmf = ccl.halos.MassFuncTinker08(mass_def=m200def)		 # Halo mass function
hbf = ccl.halos.HaloBiasTinker10(mass_def=m200def)		 # Halo bias function
conc = ccl.halos.ConcentrationDuffy08(mass_def=m200def)	 # Halo conc. function
# Feed this information to a HaloMassCalculator
hmc = ccl.halos.HMCalculator(
    mass_function=hmf,
    halo_bias=hbf,
    mass_def=m200def
    )
# 2pt correlator
prof2pt = ccl.halos.Profile2ptHOD()

# For redshift-dependent parameters
a_pivot = 1 / (1 + cf.z_pivot)

#############################
#         FUNCTIONS         #
#############################


def get_lmax():
    '''
    Computes a suitable scale cut to apply prior to fitting, depending on the
    effective redshifts of the tomographic bins.

    Returns
    -------
    cuts: dict
        Dictionary of lmax to use for each bin pairing.
    '''
    # Set up dictionary for the scale cuts
    lmax_dict = {}
    # Compute comoving distance for each tomographic bin included in pairings
    for i in np.unique(pairings):
        z = tracers[i].z
        nz = tracers[i].nz
        # Effective redshift of bin
        zeff = np.sum(z * nz) / np.sum(nz)
        # Comoving distance (Mpc)
        chi = cosmo.comoving_radial_distance(1. / (1. + zeff))
        # Ell corresponding to kmax (Mpc^{-1}) at zeff
        lmax = cf.kmax * chi
        # Use this value only if it is less than 2*NSIDE
        lmax_dict[i] = min(lmax, 2*cf.nside_hi)
    # Now cycle through each bin pairing
    cuts = []
    for p in pairings:
        i, j = p
        # Take the minimum lmax for each bin pairing
        cuts.append(min(lmax_dict[i], lmax_dict[j]))
    return cuts


def log_prior(theta):
    '''
    Defines the priors on the free parameters in the HOD model
    (in this case, mu_min and mu_1). NOTE: assumes uniform priors.

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
    logp = 0.
    for i, p in enumerate(cf.cobaya_info.params):
        if (theta[i] <= cf.cobaya_info.params[p]['prior']['min'])\
         or (theta[i] >= cf.cobaya_info.params[p]['prior']['max']):
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
        Wavenumber (in Mpc^{-1}) below which the 1-halo term will be
        suppressed.
    '''
    k_max = 0.01
    return k_max


def log_likelihood(mu_min, mu_1, mu_min_p, mu_1_p, alpha_smooth):
    '''
    Defines the priors on the free parameters in the HOD model.

    Parameters
    ----------
    mu_min: float
        log_10 of the constant component of M_min, where M_min is the
        characteristic mass above which halos can host central galaxies.

    mu_1: float
        log_10 of the constant component of M_1, where M_1 is a normalisation
        factor for the satellite occupation distribution.

    mu_minp: float
        log_10 of the time-evolving component of M_min.

    mu_1p: float
        log_10 of the time-evolving component of M_1.

    alpha_smooth: float
        Smoothing between the 1-halo and 2-halo regimes of the matter power
        spectrum.

    Returns
    -------
    logL: float
        Logarithm of the likelihood.
    '''
    def smooth_transition(a):
        return alpha_smooth

    # Halo profile
    prof = ccl.halos.HaloProfileHOD(
        mass_def=m200def,
        concentration=conc,
        log10Mmin_0=mu_min,
        log10M0_0=mu_min,
        log10M1_0=mu_1,
        log10M0_p=mu_min_p,
        log10Mmin_p=mu_min_p,
        log10M1_p=mu_1_p,
        a_pivot=a_pivot
        )
    # Halo-model power spectrum for galaxies
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
    # Compute theory C_ells
    theory_cells = [
        ccl.angular_cl(
            cosmo,
            NCT[i],
            NCT[j],
            ells[ip],
            p_of_k_a=pk
        )
        for ip, (i, j) in enumerate(pairings)
        ]

    # Residuals
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
    # Compute the logged prior
    lp = log_prior(theta)
    if not np.isfinite(lp):
        log_prob = -np.inf
    else:
        log_prob = lp + log_likelihood(*theta)
    return log_prob


def nll(*args):
    '''
    Converts log_probability into a function that can be minimised
    to find the best-fit parameters.
    '''
    return -log_probability(*args)

#######################################################
#                  START OF SCRIPT                    #
#######################################################


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Cycle through the fields being analysed
for fd in cf.fields:

    PATH_FD = cf.paths.out + fd + '/'

    # Load the Sacc file containing the power spectrum info
    s = select_from_sacc(PATH_FD + cf.sacc_files.main, label_pairs, 'cl_00')
    # Get ells and cells (no scale cuts applied at this stage)
    ell_cl = [s.get_ell_cl('cl_00', i, j)
              for i, j in s.get_tracer_combinations()]
    ells_all = [ell_cl[i][0] for i in range(len(ell_cl))]
    cells_all = [ell_cl[i][1] for i in range(len(ell_cl))]
    cov_all = s.covariance.covmat

    # Construct NumberCountsTracer objects from the saved n(z) info
    tracers = [s.tracers[i] for i in s.tracers.keys()]
    NCT = [
        ccl.NumberCountsTracer(
            cosmo,
            has_rsd=False,
            dndz=(t.z, t.nz),
            bias=(t.z, np.ones_like(t.z))
        ) for t in tracers
    ]

    # Compute scale cuts from k_max?
    if cf.compute_lmax:
        lmax_cuts = get_lmax()
    else:
        lmax_cuts = [cf.lmax for _ in range(ncombos)]
    # Apply scale cuts
    ells, cells, cov = cu.apply_scale_cuts(
                                        ells_all,
                                        cells_all,
                                        cov_all,
                                        cf.lmin,
                                        lmax_cuts
                                        )
    # Invert the covariance matrix
    icov = np.linalg.inv(cov)
    print('Success!')
    exit()
    if rank == 0:
        print(colour_string(fd.upper(), 'orange'))
        # Retrieve cobaya config options from config file
        info = cf.cobaya_info
        # Add the likelihood
        info['likelihood'] = {'hod': log_likelihood}

        # Perform initial fit to get initial positions of each walker
        initial = [info['params'][p]['ref'] for p in info['params']]
        ndim = len(initial)
        if cf.compute_initial:
            print('Estimating intial best fit...')
            initial = minimize(nll, initial).x
        # Otherwise search for values in the config file
        else:
            print('Using initial guesses from config file...')

        print('Initial best-fit values:\n')
        for i, p in enumerate(info['params']):
            info['params'][p]['ref'] = initial[i]
            print(f'{p} = {initial[i]:.3f}')

        # Ensure chains are saved to the correct place
        info['output'] = f'{cf.paths.out}{fd}/{info["output"]}{cf.suffix}'
    else:
        info = None
    info = comm.bcast(info, root=0)

    # Run the sampler
    success = False
    try:
        updated_info, sampler = run(info)
        success = True
    except LoggedError:
        pass

    # Did it work? (e.g. did not get stuck)
    success = all(comm.allgather(success))

    if not success and rank == 0:
        print('Sampling failed!')
    print('Sampling successful!')
