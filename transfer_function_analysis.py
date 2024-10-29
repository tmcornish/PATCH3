#####################################################################################################
# - Uses simulated maps to compare deprojection bias vs transfer function as a
#   means of obtaining unbiased angular power spectra. This is a standalone script and
#   not intended as part of the pipeline, but may be worked into it eventually.
#####################################################################################################

import numpy as np
import h5py
import matplotlib.pyplot as plt
import pymaster as nmt
import config
import plot_utils as pu

cf = config.fitHods
plt.style.use(pu.styledict)

#######################################################
###############    START OF SCRIPT    #################
#######################################################

nsim = 100
PATH_FD = cf.PATH_OUT + 'combined/'
PATH_SIMS = PATH_FD + 'sims/'

#load workspace
w = nmt.NmtWorkspace()
w.read_from(PATH_SIMS + 'sims_workspace.fits')

#set up lists for the results of each method
cl_best_db = []
err_cl_best_db = []
cl_ratio = []
#cycle through the simulations
for i in range(len(nsim)):
	simfile = f'{PATH_SIMS}sim{i}_nside{cf.nside_hi}.hdf5'
	with h5py.File(simfile, 'r') as hf:
		#get the input C_ell (same in all output files)
		if i == 0:
			cl_in = hf['cl_in'][:]
		#calculate the best estimate of C_ell using deprojection bias
		cl_best = hf['cl_meas'][:] - hf['cl_bias'][:]
		cl_best_db.append(cl_best)
		#append the uncertainties on C_ell to the list
		err_cl_best_db.append(hf['err_cell'][:])
		#calculate the ratio of the measured to input C_ell
		r_cl = w.decouple_cell(hf['cl_meas_coupled'] / cl_in)
		cl_ratio.append(r_cl)

