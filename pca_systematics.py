##############################################################################
# - Uses principal component analysis to create a new set of orthogonal
#   templates for systematics mitigation.
##############################################################################

import os
import sys
from configuration import PipelineConfig as PC
import healpy as hp
import healsparse as hsp
import numpy as np
from map_utils import MaskData
from output_utils import colour_string
import glob

# SETTINGS #
config_file = sys.argv[1]
cf = PC(config_file, stage='pcaSystematics')

if cf.plot_eigen:
    import plot_utils as pu
    from matplotlib import pyplot as plt
    plt.style.use(pu.styledict)


###################
#    FUNCTIONS    #
###################

def load_systmap_hsp(map_path, mask):
    '''
    Loads a HealSparse map, and uses the mask to retrieve only the valid
    pixels, and also compute and subtract the mean.

    Parameters
    ----------
    map_path: str
        Path to the HealSparse map being loaded.

    mask: MaskData
        MaskData object containing all information about the mask.

    Returns
    -------
    systmap: np.array
        Systematics map data in all valid pixels.
    '''
    # Load the systematics map
    systmap = hsp.HealSparseMap.read(map_path)
    # Get the valid pixels from the mask (NEST ordering)
    vpix = mask.vpix_nest
    # Keep only these pixels and convert any hp.UNSEEN or NaNs to zero
    systmap = systmap[vpix].astype(np.float64)
    systmap[systmap == hp.UNSEEN] = 0.
    systmap[np.isnan(systmap)] = 0.
    # Calculate the weighted mean
    mu = np.sum(systmap * mask.mask[mask.vpix]) / mask.sum
    # Subtract this from the data
    systmap -= mu
    return systmap


#######################################################
#                  START OF SCRIPT                    #
#######################################################

for fd in cf.fields:
    print(colour_string(fd.upper(), 'orange'))

    # Relevant directories
    PATH_MAPS = cf.paths.out + f'{fd}/'
    PATH_SYST = PATH_MAPS + 'systmaps/'
    PATH_PCA = PATH_SYST + 'pca/'

    # Load the mask data
    MD = MaskData(PATH_MAPS + f'survey_mask_{cf.nside_hi}.hsp',
                  mask_thresh=cf.weight_thresh)
    mask = MD.mask[MD.vpix]
    # Get the number of (valid) pixels
    Npix = len(MD.vpix_nest)

    # Get the filenames of all systematics maps at the correct NSIDE
    systs = [os.path.basename(m)
             for m in (glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp')
                       + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp'))]
    Nsyst = len(systs)

    # Set up an array for the data
    data_orig = np.array(np.zeros((Npix, Nsyst)))
    # Load each systematic and store in the array
    for i, s in enumerate(systs):
        sm = load_systmap_hsp(PATH_SYST + s, mask=MD)
        data_orig[:, i] = sm

    # Standardise the data
    data_orig /= data_orig.std(axis=0)

    # Compute the weighted covariance matrix
    cov = (1 / MD.sum) * data_orig.T * mask[None, :]
    cov = np.matmul(cov, data_orig)

    # Compute the eigenvectors and eigenvalues
    eigvals, eigvecs = np.linalg.eig(cov)

    # Sort the eigenvalues and eigenvectors by order of decreasing eigenvalue
    idx_sort = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx_sort]
    eigvecs = eigvecs[:, idx_sort]

    if cf.plot_eigen:
        f, ax = plt.subplots()
        ax.set_xlabel('Principal component')
        ax.set_ylabel('Fraction of total variance')

        # Along the x axis just plot integer label for each principal component
        X = np.arange(Nsyst).astype(int) + 1
        # Plot the fraction of the total variance contained within each PC
        Y = eigvals.real / eigvals.real.sum()

        ax.plot(X, Y)

        # Determine the number of PCs required to contain the desired fraction
        # of the total variance
        PC_max = np.argwhere(np.cumsum(Y) > cf.var_thresh).min()

        # Plot the cutoff as a vertical line
        ax.axvline(PC_max + 0.5, c=pu.red, linestyle='--')
        ax.text(PC_max+1.7, Y.max()*0.95,
                f'{int(100*cf.var_thresh)}% variance',
                rotation=0, va='center', ha='left', color=pu.red)

        f.tight_layout()
        f.savefig(cf.PATH_PLOTS + f'systematics_pca_variance_{fd}.png',
                  dpi=300)
