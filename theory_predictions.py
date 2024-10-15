#####################################################################################################
# - Estimates the n(z) distributions and then computes theoretical predictions 
#   for the power spectra.
#####################################################################################################

import pyccl as ccl
import h5py
import numpy as np
from output_utils import colour_string
import config

### SETTINGS ###
cf = config.theoryPredictions


###################
#### FUNCTIONS ####
###################

def compute_nofz(fd):
	'''
	Computes estimates of the n(z) distributions in each redshift bin for a given
	field and returns them in a dictionary.

	Parameters
	----------
	fd: str
		Name of the field being analysed.
	
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
		for f in ['hectomap', 'spring', 'autumn']:
			with h5py.File(cf.PATH_OUT + f + '/' + cf.cat_main, 'r') as hf:
				gr = hf['photometry']
				z_best.append(gr[f'{cf.zcol}'][:])
				z_mc.append(gr[f'{cf.z_mc_col}'][:])
	else:
		with h5py.File(cf.PATH_OUT + fd + '/' + cf.cat_main, 'r') as hf:
			gr = hf['photometry']
			z_best.append(gr[f'{cf.zcol}'][:])
			z_mc.append(gr[f'{cf.z_mc_col}'][:])
	
	#concatenate the lists of arrays
	z_best = np.concatenate(z_best)
	z_mc = np.concatenate(z_mc)

	#determine the bin edges and centres to use for the n(z) histograms
	bins = np.arange(0., z_mc.max()+cf.dz, cf.dz)
	bin_centres = (bins[1:] + bins[:-1]) / 2
	
	#set up a dictionary for containing the n(z)s for the different bins
	nofz = {'z' : bin_centres}

	#generate the histograms and store them in the dictionary
	for i in range(cf.nbins):
		zmask = (z_best >= cf.zbins[i]) * (z_best < cf.zbins[i+1])
		nofz[f'nz_{i}'] = np.histogram(z_mc[zmask], bins=bins, density=True)[0]
	
	return nofz


def get_tracers(nofz, cosmo):
	'''
	Creates and returns a dictionary of NumberCountsTracer objects for
	each tomographic bin.

	Parameters
	----------
	nofz: dict
		Dictionary containing n(z) distributions for each bin and the
		redshifts at which they are defined.
	
	cosmo: ccl.Cosmology
		Fiducial cosmology object.
	
	Returns
	-------
	tracers: dict
		Dictionary containing tracer information for each bin.
	'''

	#create a dictionary containing the tracers in each redshift bin
	tracers = {
		f'bin{i}' : ccl.NumberCountsTracer(
											cosmo, 
											has_rsd=False, 
											dndz=(nofz['z'], nofz[f'nz_{i}']), 
											bias=(nofz['z'], np.ones_like(nofz['z']))
											)
		for i in range(cf.nbins)
	}

	return tracers


def get_theory_cells(tracers, pairings, ells, cosmo, pk):
	'''
	Computes theoretical predictions for the angular power spectra for the 
	specified pairs of tomographic bins.

	Parameters
	----------
	tracers: dict
		Dictionary containing NumberCountsTracer objects for each tomographic
		bin.
	
	pairings: list
		List of tuples, with each tuple containing the indices of the bins
		being paired.
	
	ells: array-like
		Multipoles at which the theory C_ells will be defined.
	
	cosmo: ccl.Cosmology
		Fiducial cosmology.
	
	pk: ccl.halos.halomod_Pk2D
		Galaxy halo-model power spectrum.
	
	Returns
	-------
	theory_cells: dict
		Dictionary containing the theoretical angular power spectra for each
		bin pair, as well as the multipoles at which they are defined.

	'''
	theory_cells = {'ells' : ells}
	for i,j in pairings:
		c_ell = ccl.angular_cl(
							cosmo,
							tracers[f'bin{i}'],
							tracers[f'bin{j}'],
							ells,
							p_of_k_a=pk
						)
		theory_cells[f'bin{i}-bin{j}'] = c_ell
	
	return theory_cells


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

#theory c_ell settings
lmax = 3. * cf.nside_hi - 1
ells = np.unique(np.geomspace(1, lmax, 128).astype(int)).astype(float)

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

#if using DIR outputs, load them now
if cf.use_dir:
	print('Retrieving DIR n(z) distributions...')
	with h5py.File(cf.nz_dir_file, 'r') as hf:
		nofz = {k : hf[k][:] for k in hf.keys()}
	#create a dictionary containing the tracers in each redshift bin
	tracers = get_tracers(nofz, cosmo)
	print('Computing theory C_ells...')
	#compute theory C_ells
	theory_cells = get_theory_cells(tracers, pairings, ells, cosmo, pk)
	
#cycle through the specified fields
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))
	if not cf.use_dir:
		print('Computing n(z) distributions from Monte-Carlo draws...')
		#compute estimates of the n(z) distributions in each bin
		nofz = compute_nofz(fd)
		#save the n(z) info to a file
		outfile = f'{cf.PATH_OUT}{fd}/{cf.nz_mc_file}'
		with h5py.File(outfile, 'w') as hf:
			for k in nofz.keys():
				hf.create_dataset(k, data=nofz[k])
		#create a dictionary containing the tracers in each redshift bin
		tracers = get_tracers(nofz, cosmo)
		print('Computing theory C_ells...')
		#compute theory C_ells
		theory_cells = get_theory_cells(tracers, pairings, ells, cosmo, pk)

	print('Saving theory C_ells...')
	#open the output file, compute and save the theory cells
	outfile = f'{cf.PATH_OUT}{fd}/{cf.theory_out}' 
	with h5py.File(outfile, 'w') as psfile:
		for k in theory_cells.keys():
			psfile.create_dataset(k, data=theory_cells[k])
	

