#####################################################################################################
# - Fits a halo occupation distribution model to the measured angular power spectra.
#####################################################################################################

import config
import emcee
import sacc
import numpy as np
import itertools
import sys
import pyccl as ccl

### SETTINGS ###
cf = config.fitHods
#fiducial comsology
cosmo = ccl.Cosmology(**cf.cosmo_fiducial)

#######################################################
###############    START OF SCRIPT    #################
#######################################################


