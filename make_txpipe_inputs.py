###############################################################################
# Converts relevant catalogues and maps into the correct formats for TXPipe.
# Currently only set up to convert the systematics maps, as these are the only
# ones not created in TXPipe.
###############################################################################

import os
import sys
from configuration import PipelineConfig as PC
import h5py
import healsparse as hsp
import numpy as np
import map_utils as mu
import glob


# SETTINGS #
config_file = sys.argv[1]
cf = PC(config_file, stage='makeTXPipeInputs')


###################
#    FUNCTIONS    #
###################

def make_tomography_cat(cat, fname):
    '''
    Makes an HDF file containing the information required by TXPipe's
    TXLensMaps stage. NOTE: This information includes data described as
    'lens_weight', which is not relevant for this study. It is therefore
    assumed here that this quantity is 1 for all sources.

    Parameters
    ----------
    cat: h5py Dataset or Group
        Catalogue of galaxies, previously sorted into tomographic bins.

    fname: str
        Filename to be given to the output file. If the file already exists,
        will try to append to existing data in the file.
    '''

    # Create an array to contain the counts in each bin
    counts = np.zeros(cf.nsamples, dtype='i8')

    # Create a 1D array where galaxies can be flagged with integers if
    # belonging to a bin
    nsources = len(cat[cat.keys()[0]][:])
    zbin_flags = np.full(nsources, -1).astype(int)
    # Cycle through the redshift bins
    for i, samp in enumerate(cf.samples):
        zmask = cat[samp][:]
        zbin_flags[zmask] = i
        counts[i] = zmask.sum()

    # Create an array for containing the total counts in all bins
    counts_2d = np.array([counts.sum()])

    # Create a column for 'lens_weight' and simply assign a value of 1 for
    # every source
    lens_weight = np.ones(nsources, dtype='f8')

    # Compile the relevant data into a list and assign them names
    data = [zbin_flags, counts, counts_2d, lens_weight]
    names = ['bin', 'counts', 'counts_2d', 'lens_weight']

    # Write to the output file
    with h5py.File(fname, 'w') as hf:
        g = hf.create_group('tomography')
        g.attrs['nbin'] = cf.nsamples
        for d, n in zip(data, names):
            g.create_dataset(n, data=d, shape=(len(d),), dtype=d.dtype)


#######################################################
#                  START OF SCRIPT                    #
#######################################################

# Cycle throught the fields
for fd in cf.fields:
    # Path to directory containing the outputs for this field
    OUT = cf.paths.out + fd + '/'
    # Systematics maps directory
    PATH_SYST = OUT + 'systmaps/'
    # Path in which the reformatted systematics maps will be placed
    PATH_TX = PATH_SYST + 'for_txpipe/'
    if not os.path.exists(PATH_TX):
        os.system(f'mkdir -p {PATH_TX}')

    # Bright object mask, to be written to HDF5
    bo_mask = hsp.HealSparseMap.read(f'{OUT}{cf.maps.bo_mask}')
    mu.healsparseToHDF(bo_mask,
                       f'{PATH_TX}{cf.maps.bo_mask[:-4]}.hdf5',
                       group='maps/mask')

    # Survey mask, also to be written to HDF5
    survey_mask = hsp.HealSparseMap.read(f'{OUT}{cf.maps.survey_mask}')
    # Calculate area above the mask threshold and fractional sky coverage
    A_unmasked, f_sky = mu.maskAreaSkyCoverage(survey_mask,
                                               thresh=cf.weight_thresh)
    mask_meta = {'area': A_unmasked, 'f_sky': f_sky}
    mu.healsparseToHDF(survey_mask,
                       f'{OUT}/{cf.maps.survey_mask[:-4]}.hdf5',
                       group='maps/mask',
                       metadata=mask_meta)

    # Systematics maps
    maps = sorted(glob.glob(f'{PATH_SYST}*_nside{cf.nside_hi}*.hsp'))\
        + [f'{OUT}']
    # Cycle through the maps
    for m in maps:
        # Load the HealSparse map
        hsp_map = hsp.HealSparseMap.read(m)
        # Retrieve the basename of the file
        bn = os.path.basename(m)
        # See if the input map contains a recarray
        if len(hsp_map.dtype) > 0:
            for n in hsp_map.dtype.names:
                fname = f'{PATH_TX}{bn[:-4]}_{n}.fits'
                mu.healsparseToFITS(hsp_map[n], fname, nest=False)
        else:
            fname = f'{PATH_TX}{bn[:-4]}.fits'
            mu.healsparseToFITS(hsp_map, fname, nest=False)

    # Load the cleaned galaxy catalogue
    galcat = f'{OUT}{cf.cats.main}'
    with h5py.File(galcat, 'r') as hf:
        # Retrieve the tomographic bin labels
        zbin_flags = galcat['photometry/zbin'][:]
    # Save relevant data to a tomographic catalogue
    tomo_out = PATH_TX + cf.cats.tomography
    make_tomography_cat(zbin_flags, tomo_out)
