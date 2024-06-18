#####################################################################################################
# - Estimates the n(z) distributions and then computes theoretical predictions 
#   for the power spectra.
#####################################################################################################

import pyccl as ccl
import h5py
import numpy as np
import config

### SETTINGS ###
cf = config.theoryPredictions


###################
#### FUNCTIONS ####
###################

def get_nofz(fd, group='', zlims=None):
	'''
	Computes estimates of the n(z) distributions in each redshift bin for a given
	field and returns them in a dictionary.

	Parameters
	----------
	fd: str
		Name of the field being analysed.
	
	group: str
		Group within which the relevant data are expected to reside.
	
	Returns
	-------
	nofz: dict
		Dictionary containing the redshift bin centres and the n(z) values in each
		of those bins for each tomographic bin.
	'''
	#define lists to contain the best z estimates and the random MC draws
	z_best, z_mc = [], []

	#if provided field is 'combined', need to load data from all fields
	if fd == 'combined':
		for f in ['hectomap', 'equatora', 'equatorb']:
			with h5py.File(cf.PATH_OUT + f + '/' + cf.cat_main, 'r') as hf:
				z_best.append(hf[f'{group}/{cf.zcol}'][:])
				z_mc.append(hf[f'{group}/{cf.z_mc_col}'][:])
	else:
		with h5py.File(cf.PATH_OUT + fd + '/' + cf.cat_main, 'r') as hf:
			z_best.append(hf[f'{group}/{cf.zcol}'][:])
			z_mc.append(hf[f'{group}/{cf.z_mc_col}'][:])
	
	#concatenate the lists of arrays
	z_best = np.concatenate(z_best)
	z_mc = np.concatenate(z_mc)

	#determine the bin edges and centres to use for the n(z) histograms
	if zlims is None:
		bins = np.linspace(z_best.min(), z_best.max(), (cf.nbins_nofz + 1))
	else:
		bins = np.linspace(*zlims, (cf.nbins_nofz + 1))
	bin_centres = (bins[1:] + bins[:-1]) / 2
	
	#set up a dictionary for containing the n(z)s for the different bins
	nofz = {'z' : bin_centres}

	#generate the histograms and store them in the dictionary
	for i in range(len(cf.zbins) - 1):
		zmask = (z_best >= cf.zbins[i]) * (z_best < cf.zbins[i+1])
		nofz[f'bin{i}'] = np.histogram(z_mc[zmask], bins=bins, density=True)[0]
	
	return nofz


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#retrieve the bin pairings
pairings, pairings_s = cf.get_bin_pairings()

#define the fiducial cosmology
cosmo = ccl.Cosmology(**cf.cosmo_fiducial)

#range of wavenumbers and scale factors over which theory power spectra will be computed
lk_arr = np.log(np.geomspace(1E-4, 100, 256))
a_arr = 1. / (1. + np.linspace(0, 6, 100)[::-1])

#cycle through the specified fields
for fd in cf.get_global_fields():
	#retrieve estimates of the n(z) distributions in each bin
	nofz = get_nofz(fd, group='photometry')
	
	#create a dictionary containing the tracers in each redshift bin
	tracers = {
		f'bin{i}' : ccl.NumberCountsTracer(
											cosmo, 
											has_rsd=False, 
											dndz=(nofz['z'], nofz[f'bin{i}']), 
											bias=(nofz['z'], np.ones_like(nofz['z']))
											)
		for i in range(len(cf.zbins) - 1)
	}

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
	#create halo profile
	prof = ccl.halos.HaloProfileHOD(mass_def=m200def, concentration=conc)
	#2pt correlator
	prof2pt = ccl.halos.Profile2ptHOD()

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

	#compute the theory c_ells
	lmax = 3. * cf.nside_hi - 1
	ells = np.unique(np.geomspace(1, lmax, 128).astype(int)).astype(float)
	#open the output file, compute and save the theory cells
	outfile = f'{cf.PATH_OUT}{fd}/{cf.theory_out}' 
	with h5py.File(outfile, 'w') as psfile:
		_ = psfile.create_dataset('ells', data=ells)
		for i, j in pairings:
			c_ells = ccl.angular_cl(
				cosmo,
				tracers[f'bin{i}'],
				tracers[f'bin{j}'],
				ells,
				p_of_k_a=pk
			)
			_ = psfile.create_dataset(f'bin{i}-bin{j}', data=c_ells)
	

