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
pairings, pairings_s = cf.get_bin_pairings()

#ylabel for figure(s) depends on normalisation of the C_ells
if cf.normalise:
    ylabel = r'$C_{\ell}\frac{\ell(\ell+1)}{2\pi}$'
else:
    ylabel = r'$C_{\ell}$'

#if told to make a figure showing all results, set it up
if cf.make_combined:
    fig_comb = plt.figure(figsize=(xsize, ysize))
    gs_comb = fig_comb.add_gridspec(ncols=nbins, nrows=nbins)
    axes_comb = [fig_comb.add_subplot(gs_comb[j,i]) for i,j in pairings]

    #set up lists to which data, handles and labels will be appended
    ell_effs_all = []
    Cell_all = []
    Cell_debiased_all = []
    Nell_all = []
    err_cell_all = []
    bias_all = []
    handles_all = []
    labels_all = []

    #colours and markers to cycle through
    colours_cell = ['#0000c9', '#94000d', '#048101', pu.green]
    markers_cell = ['o', 's', '^', 'v']
    colours_nell = ['#0e78d1', '#cd2b70', '#39b523', pu.ruby]
    colours_bias = ['#34dff0', '#ec52af', '#61f33e', pu.coral]

    

#cycle through the fields being analysed 
for idx, fd in enumerate(cf.get_global_fields()):
    print(colour_string(fd.upper(), 'orange'))

    #set up a figure for the power spectra from each redshift bin
    fig = plt.figure(figsize=(xsize, ysize))
    gs = fig.add_gridspec(ncols=nbins, nrows=nbins)


    #full path to the output file
    outfile_main = f'{cf.PATH_OUT}{fd}/{cf.outfile}'

    #cycle through all possible pairings of redshift bins
    for ip, (p, p_str) in enumerate(zip(pairings, pairings_s)):
        i,j = p
        outfile = f'{outfile_main[:-5]}_{i}_{j}.hdf5'
        #open the file, creating it if it doesn't exist
        with h5py.File(outfile, mode='r') as psfile:
            print(colour_string(p_str, 'green'))
            
            #retrieve the relevant information
            gp = psfile[p_str]
            ell_effs = gp['ell_effs'][:]
            Cell = gp['cl_decoupled'][:]
            Cell_debiased = gp['cl_decoupled_debiased'][:]
            Cell_no_deproj = gp['cl_decoupled_no_deproj'][:]
            bias = gp['cl_bias_decoupled'][:]
            err_cell = gp['err_cell'][:]
            err_cell_no_deproj = gp['err_cell_no_deproj'][:]
            Nell = gp['N_ell_decoupled'][:]

            #define the x-limits using the effective ells
            xmin = ell_effs.min() / 1.5
            xmax = ell_effs.max() * 1.2

            if cf.normalise:
                mfactor = ell_effs * (ell_effs + 1) / (2. * np.pi)
            else:
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
            mask_pve = Cell_debiased[0] > 0
            mask_nve = Cell_debiased[0] <= 0
            #plot the shot noise if autocorrelation
            if i == j:
                Y_pve = (Cell_debiased[0][mask_pve] - Nell[0][mask_pve]) * mfactor[mask_pve]
                Y_nve = (Cell_debiased[0][mask_nve] - Nell[0][mask_nve]) * mfactor[mask_nve]
                noise_plot, *_ = ax.plot(ell_effs, Nell[0]*mfactor, c=pu.teal)
            else:
                Y_pve = Cell_debiased[0][mask_pve] * mfactor[mask_pve]
                Y_nve = Cell_debiased[0][mask_nve] * mfactor[mask_nve]
            cell_plot = ax.errorbar(ell_effs[mask_pve], Y_pve, yerr=err_cell[mask_pve], marker='o', c=pu.dark_blue, linestyle='none')
            ax.errorbar(ell_effs[mask_nve], -Y_nve, yerr=err_cell[mask_nve], marker='o', markeredgecolor=pu.dark_blue, markerfacecolor='none', linestyle='none')

            if cf.show_pre_debias:
                #plot the power spectrum pre-deprojection, using open symbols for abs(negative) values
                mask_pve_pre = Cell[0] > 0
                mask_nve_pre = Cell[0] <= 0
                #subtract the shot noise if autocorrelation
                if i == j:
                    Y_pve_pre = (Cell[0][mask_pve_pre] - Nell[0][mask_pve_pre]) * mfactor[mask_pve_pre]
                    Y_nve_pre = (Cell[0][mask_nve_pre] - Nell[0][mask_nve_pre]) * mfactor[mask_nve_pre]
                else:
                    Y_pve_pre = Cell[0][mask_pve_pre] * mfactor[mask_pve_pre]
                    Y_nve_pre = Cell[0][mask_nve_pre] * mfactor[mask_nve_pre]
                cell_plot_pre = ax.errorbar(ell_effs[mask_pve_pre]*1.05, Y_pve_pre, yerr=err_cell[mask_pve_pre], marker='o', c=pu.dark_blue, ecolor=pu.dark_blue, linestyle='none', alpha=0.4)
                ax.errorbar(ell_effs[mask_nve_pre]*1.05, -Y_nve_pre, yerr=err_cell[mask_nve_pre], marker='o', markeredgecolor=pu.dark_blue, ecolor=pu.dark_blue, markerfacecolor='none', linestyle='none', alpha=0.4)

            if cf.show_no_deproj:
                #plot the power spectrum pre-deprojection, using open symbols for abs(negative) values
                mask_pve_nd = Cell_no_deproj[0] > 0
                mask_nve_nd = Cell_no_deproj[0] <= 0
                #subtract the shot noise if autocorrelation
                if i == j:
                    Y_pve_nd = (Cell_no_deproj[0][mask_pve_nd] - Nell[0][mask_pve_nd]) * mfactor[mask_pve_nd]
                    Y_nve_nd = (Cell_no_deproj[0][mask_nve_nd] - Nell[0][mask_nve_nd]) * mfactor[mask_nve_nd]
                else:
                    Y_pve_nd = Cell_no_deproj[0][mask_pve_nd] * mfactor[mask_pve_nd]
                    Y_nve_nd = Cell_no_deproj[0][mask_nve_nd] * mfactor[mask_nve_nd]
                cell_plot_nd = ax.errorbar(ell_effs[mask_pve_nd]*1.05, Y_pve_nd, yerr=err_cell_no_deproj[mask_pve_nd], marker='o', c=pu.dark_blue, ecolor=pu.dark_blue, linestyle='none', alpha=0.4)
                ax.errorbar(ell_effs[mask_nve_nd]*1.05, -Y_nve_nd, yerr=err_cell[mask_nve_nd], marker='o', markeredgecolor=pu.dark_blue, ecolor=pu.dark_blue, markerfacecolor='none', linestyle='none', alpha=0.4)

            #reset the axis limits
            ax.set_xlim(xmin, xmax)

            #add text to the top-right corner to indicate which bins have been compared
            ax.text(0.95, 0.95, f'({i},{j})', transform=ax.transAxes, ha='right', va='top', fontsize=20.)

            if cf.make_combined:
                #get subplot axes
                ax_comb = axes_comb[ip]
                if idx == 0:
                    #only label axes if on outer edge of figure
                    if j == (nbins-1):
                        ax_comb.set_xlabel(r'$\ell$')
                    if i == 0:
                        ax_comb.set_ylabel(ylabel)
                    #set loglog scale
                    ax_comb.set_xscale('log')
                    ax_comb.set_yscale('log')

                    ax_comb.set_xlim(xmin, xmax)
                
                #plot c_ells
                handles_all.append(ax_comb.errorbar(ell_effs[mask_pve], Y_pve, yerr=err_cell[mask_pve], marker=markers_cell[idx], c=colours_cell[idx], linestyle='none')[0])
                ax_comb.errorbar(ell_effs[mask_nve], -Y_nve, yerr=err_cell[mask_nve], marker=markers_cell[idx], markeredgecolor=colours_cell[idx], markerfacecolor='none', linestyle='none')
                labels_all.append(f'Signal ({fd})')

                #plot n_ells
                if i == j:
                    handles_all.append(ax_comb.plot(ell_effs, Nell[0]*mfactor, c=colours_nell[idx])[0])
                    labels_all.append(f'Noise ({fd})')
                
                #plot deprojection bias
                if bias.any():
                    handles_all.append(ax_comb.plot(ell_effs, bias[0]*mfactor, c=colours_bias[idx])[0])
                    ax_comb.plot(ell_effs, -bias[0]*mfactor, ls='--', c=colours_bias[idx])
                    labels_all.append(f'Deproj. bias ({fd})')
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
        figname = f'{cf.PATH_PLOTS}power_spectra_{fd}_{cf.nside_hi}.png'
    else:
        figname = f'{cf.PATH_PLOTS}power_spectra_raw_{fd}_{cf.nside_hi}.png'

    if cf.show_pre_debias:
        handles.insert(1, cell_plot_pre)
        labels.insert(1, 'Biased signal')
    
    if cf.show_no_deproj:
        handles.insert(1, cell_plot_nd)
        labels.insert(1, 'Signal pre-deproj.')

    fig.legend(handles=handles, labels=labels, loc='upper right', fontsize=28)

    fig.tight_layout()
    fig.savefig(figname, dpi=300)


if cf.make_combined:
    by_label = dict(zip(labels_all, handles_all))
    fig_comb.legend(by_label.values(), by_label.keys(), loc='upper right', ncols=idx+1, fontsize=19)

    figname = f'{cf.PATH_PLOTS}power_spectra_all_{cf.nside_hi}.png'
    fig_comb.tight_layout()
    fig_comb.savefig(figname, dpi=300)