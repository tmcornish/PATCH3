#####################################################################################################
# - 
#####################################################################################################

import config
import healpy as hp
import healsparse as hsp
import numpy as np
from map_utils import *
from output_utils import colour_string
import os
import glob

### SETTINGS ###
cf = config.pcaSystematics

if cf.plot_eigen:
    import plot_utils as pu
    from matplotlib import pyplot as plt
    plt.style.use(pu.styledict)


###################
#### FUNCTIONS ####
###################

def load_systmap_hsp(map_path, mask):
    '''
    Loads a HealSparse map, and uses the mask to retrieve only the valid pixels,
    and also compute and subtract the mean. 
    
    Parameters
    ----------
    map_path: str
        Path to the HealSparse map being loaded.
    
    mask: MaskData
        MaskData object containing all information about the mask.
    
    Returns
    -------
    systmap: np.array
        Systematics map data in all valid pixels, 
    '''
    #load the systematics map
    systmap = hsp.HealSparseMap.read(map_path)
    #get the valid pixels from the mask (NEST ordering)
    vpix = mask.vpix_nest
    #keep only these pixels and convert any hp.UNSEEN or NaNs to zero
    systmap = systmap[vpix].astype(np.float64)
    systmap[systmap == hp.UNSEEN] = 0.
    systmap[np.isnan(systmap)] = 0.
    #calculate the weighted mean
    mu = np.sum(systmap * mask.mask[mask.vpix]) / mask.sum
    #subtract this from the data
    systmap -= mu
    return systmap



    

#######################################################
###############    START OF SCRIPT    #################
#######################################################


for fd in cf.get_global_fields():
    print(colour_string(fd.upper(), 'orange'))

    #relevant directories
    PATH_MAPS = cf.PATH_OUT + f'{fd}/'
    PATH_SYST = PATH_MAPS + 'systmaps/'
    PATH_PCA = PATH_SYST + 'pca/'

    #load the mask data
    MD = MaskData(PATH_MAPS + f'survey_mask_{cf.nside_hi}.hsp', mask_thresh=cf.weight_thresh)
    mask = MD.mask[MD.vpix]
    #get the number of (valid) pixels
    Npix = len(MD.vpix_nest)

    #get the filenames of all systematics maps at the correct NSIDE
    systs = [os.path.basename(m) for m in (glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp'))]
    Nsyst = len(systs)

    #set up an array for the data
    data_orig = np.array(np.zeros((Npix, Nsyst)))
    #load each systematic and store in the array
    for i, s in enumerate(systs):
        sm = load_systmap_hsp(PATH_SYST + s, mask=MD)
        data_orig[:,i] = sm
    
    #standardise the data
    data_orig /= data_orig.std(axis=0)

    #compute the weighted covariance matrix
    cov = (1 / MD.sum) * data_orig.T * mask[None, :]
    cov = np.matmul(cov, data_orig)

    #compute the eigenvectors and eigenvalues
    eigvals, eigvecs = np.linalg.eig(cov)

    #sort the eigenvalues and eigenvectors by order of decreasing eigenvalue
    idx_sort = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx_sort]
    eigvecs = eigvecs[:,idx_sort]

    if cf.plot_eigen:
        f, ax = plt.subplots()
        ax.set_xlabel('Principal component')
        ax.set_ylabel('Fraction of total variance')

        #along the x axis just plot integer label for each principal component
        X = np.arange(Nsyst).astype(int) + 1
        #plot the fraction of the total variance contained within each PC
        Y = eigvals.real / eigvals.real.sum()

        ax.plot(X, Y)

        #determine the number of PCs required to contain the desired fraction of the total variance
        PC_max = np.argwhere(np.cumsum(Y) > cf.var_thresh).min()

        #plot the cutoff as a vertical line
        ax.axvline(PC_max + 0.5, c=pu.red, linestyle='--')
        ax.text(PC_max+1.7, Y.max()*0.95, f'{int(100*cf.var_thresh)}\% variance', rotation=0, va='center', ha='left', color=pu.red)

        #set y-axis limits
        #ax.set_ylim(0., Y.max()*1.02)

        f.tight_layout()
        f.savefig(cf.PATH_PLOTS + f'systematics_pca_variance_{fd}.png', dpi=300)
    
    
    

