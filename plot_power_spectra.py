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

    #set up a figure for the power spectra from each redshift bin
    fig = plt.figure(figsize=(xsize, ysize))
    gs = fig.add_gridspec(ncols=nbins, nrows=nbins)


    #full path to the output file
    outfile_main = f'{cf.PATH_OUT}{fd}/{cf.outfile}'

    #cycle through all possible pairings of redshift bins
    for ip,p in enumerate(pairings):
        outfile = f'{outfile_main[:-5]}_{i}_{j}.hdf5'
        #open the file, creating it if it doesn't exist
        with h5py.File(outfile, mode='r') as psfile:
            i,j = p
            p_str = str(p)
            print(colour_string(p_str, 'green'))
            
            #retrieve the relevant information
            gp = psfile[p_str]
            ell_effs = gp['ell_effs'][:]
            Cell_pre = gp['cl_decoupled'][:]
            Cell_post = gp['cl_decoupled_debiased'][:]
            bias = gp['cl_bias_decoupled'][:]
            err_cell = gp['err_cell'][:]
            Nell = gp['N_ell_decoupled'][:]

            #define the x-limits using the effective ells
            xmin = ell_effs.min() / 1.5
            xmax = ell_effs.max() * 1.2

            if cf.normalise:
                ylabel = r'$C_{\ell}\frac{\ell(\ell+1)}{2\pi}$'
                mfactor = ell_effs * (ell_effs + 1) / (2. * np.pi)
            else:
                ylabel = r'$C_{\ell}$'
                mfactor = np.ones_like(ell_effs)

            #add subplot to gridspec
            ax = fig.add_subplot(gs[j,i])
            #only label axes if on outer edge of figure
            if j == (nbins-1):
                ax.set_xlabel(r'$\ell$')
            if i == 0:
                ax.set_ylabel(ylabel)
            #set loglog scale
            ax.set_xscale('log')
            ax.set_yscale('log')


            if bias.any():
                #plot the deporojection bias
                bias_plot, *_ = ax.plot(ell_effs, bias[0]*mfactor, c=pu.magenta)
                ax.plot(ell_effs, -bias[0]*mfactor, ls='--', c=pu.magenta)
            
            #plot the debiased power spectrum, using open symbols for abs(negative) values
            mask_pve = Cell_post[0] > 0
            mask_nve = Cell_post[0] <= 0
            #plot the shot noise if autocorrelation
            if i == j:
                Y_pve = (Cell_post[0][mask_pve] - Nell[0][mask_pve]) * mfactor[mask_pve]
                Y_nve = (Cell_post[0][mask_nve] - Nell[0][mask_nve]) * mfactor[mask_nve]
                noise_plot, *_ = ax.plot(ell_effs, Nell[0]*mfactor, c=pu.teal)
            else:
                Y_pve = Cell_post[0][mask_pve] * mfactor[mask_pve]
                Y_nve = Cell_post[0][mask_nve] * mfactor[mask_nve]
            cell_plot = ax.errorbar(ell_effs[mask_pve], Y_pve, yerr=err_cell[mask_pve], marker='o', c=pu.dark_blue, linestyle='none')
            ax.errorbar(ell_effs[mask_nve], -Y_nve, yerr=err_cell[mask_nve], marker='o', markeredgecolor=pu.dark_blue, markerfacecolor='none', linestyle='none')

            if cf.show_pre_deproj:
                #plot the power spectrum pre-deprojection, using open symbols for abs(negative) values
                mask_pve = Cell_pre[0] > 0
                mask_nve = Cell_pre[0] <= 0
                #subtract the shot noise if autocorrelation
                if i == j:
                    Y_pve = (Cell_pre[0][mask_pve] - Nell[0][mask_pve]) * mfactor[mask_pve]
                    Y_nve = (Cell_pre[0][mask_nve] - Nell[0][mask_nve]) * mfactor[mask_nve]
                else:
                    Y_pve = Cell_pre[0][mask_pve] * mfactor[mask_pve]
                    Y_nve = Cell_pre[0][mask_nve] * mfactor[mask_nve]
                cell_plot_pre = ax.errorbar(ell_effs[mask_pve]*1.05, Y_pve, yerr=err_cell[mask_pve], marker='o', c=pu.dark_blue, linestyle='none', alpha=0.4)
                ax.errorbar(ell_effs[mask_nve]*1.05, -Y_nve, yerr=err_cell[mask_nve], marker='o', markeredgecolor=pu.dark_blue, markerfacecolor='none', linestyle='none', alpha=0.4)


            #reset the axis limits
            ax.set_xlim(xmin, xmax)

            #add text to the top-right corner to indicate which bins have been compared
            ax.text(0.95, 0.95, f'({i},{j})', transform=ax.transAxes, ha='right', va='top', fontsize=20.)

    #create a legend		
    handles = [
        cell_plot,
        noise_plot
        ]
    labels = [
        'Signal',
        'Noise'
        ]
    
    if bias.any():
        handles.insert(1, bias_plot)
        labels.insert(1, 'Deprojection bias')
        #figure name also depends on whether deprojection has occurred
        figname = f'{cf.PATH_PLOTS}{fd}_power_spectra_{cf.nside_hi}.png'
    else:
        figname = f'{cf.PATH_PLOTS}{fd}_power_spectra_raw_{cf.nside_hi}.png'

    if cf.show_pre_deproj:
        handles.insert(1, cell_plot_pre)
        labels.insert(1, 'Biased signal')

    fig.legend(handles=handles, labels=labels, loc='upper right', fontsize=28)

    plt.tight_layout()
    plt.savefig(figname, dpi=300)