#####################################################################################################
# - Reads output from compute_power_spectra and plots the results.
#####################################################################################################

import config
import numpy as np
from map_utils import *
import h5py
from matplotlib import pyplot as plt
import plot_utils as pu
from output_utils import colour_string
import itertools

### SETTINGS ###
cf = config.plotPowerSpectra
plt.style.use(pu.styledict)


#######################################################
###############    START OF SCRIPT    #################
#######################################################

#retrieve the number of redshift bins
nbins = len(cf.zbins) - 1
#use this to define the size of the power spectra figures
xsize = nbins * 4
ysize = nbins * 3.5

#also get all possible pairings of bins
l = list(range(nbins))
pairings = [i for i in itertools.product(l,l) if tuple(reversed(i)) >= i]


#cycle through the fields being analysed 
for fd in cf.get_global_fields():
	print(colour_string(fd.upper(), 'orange'))
	