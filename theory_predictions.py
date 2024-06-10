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

def get_nofz(fd, nbins=50, group='', zlims=None):
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
		bins = np.linspace(z_best.min(), z_best.max(), (nbins + 1))
	else:
		bins = np.linspace(*zlims, (nbins + 1))
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