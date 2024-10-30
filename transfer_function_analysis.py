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
ndigit = int(np.floor(np.log10(nsim))) + 1

PATH_FD = cf.PATH_OUT + 'combined/'
PATH_SIMS = PATH_FD + 'sims/'

ell_max = 3 * cf.nside_hi - 1
ells_theory = range(0, ell_max+1)

#load workspace
w = nmt.NmtWorkspace()
w.read_from(PATH_SIMS + 'sims_workspace.fits')

#set up lists for the results of each method
cl_best_db = []
err_cl_best_db = []
cl_ratio = []
#cycle through the simulations
for i in range(1,nsim+1):
	i_str = str(i).zfill(ndigit)
	id_str = f'sim{i_str}'
	simfile = f'{PATH_SIMS}{id_str}_nside{cf.nside_hi}.hdf5'
	with h5py.File(simfile, 'r') as hf:
		#get the input C_ell (same in all output files)
		if i == 1:
			ells = hf['ell_effs'][:]
			cl_in = hf['cl_in'][:]
		#calculate the best estimate of C_ell using deprojection bias
		cl_best = hf['cl_meas'][:] - hf['cl_bias'][:]
		cl_best_db.append(cl_best)
		#append the uncertainties on C_ell to the list
		err_cl_best_db.append(hf['err_cell'][:])
		#calculate the ratio of the measured to input C_ell
		r_cl = w.decouple_cell(hf['cl_meas_coupled'] / cl_in)
		cl_ratio.append(r_cl)

#convert lists to arrays
cl_best_db = np.array(cl_best_db)
err_cl_best_db = np.array(err_cl_best_db)
cl_ratio = np.array(cl_ratio)
print(cl_best_db.shape)

#compute the relevant summary statistics
cl_best_db_av = np.mean(cl_best_db, axis=0)
err_cl_best_db_av = (1 / nsim) * np.sqrt(np.sum(err_cl_best_db ** 2, axis=0))
T_ell = np.mean(cl_ratio, axis=0)

#set up a figure for plotting the results
f, ax = plt.subplots()
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$C_\ell$')

xmin = min(ells) / 1.5
xmax = max(ells) * 1.2

ax.set_xscale('log')
ax.set_yscale('log')

#plot input C_ell
ax.plot(ells_theory, cl_in, c=pu.red, label=r'$C_\ell^{\rm in}$')

#plot average best C_ell obtained via deprojection bias
mask_pve = cl_best_db_av[0] > 0
mask_nve = cl_best_db_av[0] <= 0
ax.errorbar(ells[mask_pve], cl_best_db_av[0][mask_pve], yerr=err_cl_best_db_av[mask_pve], marker='o', c=pu.dark_blue, linestyle='none', label=r'$\langle\bar{C_\ell}^{\rm deproj. bias}\rangle$')
ax.errorbar(ells[mask_nve], -cl_best_db_av[0][mask_nve], yerr=err_cl_best_db_av[mask_nve], marker='o', markeredgecolor=pu.dark_blue, markerfacecolor='none', linestyle='none')

plt.show()