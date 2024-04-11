#####################################################################################################
# Converts relevant catalogues and maps into the correct formats for TXPipe.
# Currently only set up to convert the systematics maps, as these are the only ones not created
# in TXPipe.
#####################################################################################################

import os
import config
import healpy as hp
import healsparse as hsp
import numpy as np
from output_utils import colour_string
from map_utils import *
import h5py
import glob


### SETTINGS ###
cf = config.makeTXPipeInputs



#######################################################
###############    START OF SCRIPT    #################
#######################################################


#cycle throught the fields
for fd in cf.get_global_fields():
    #path to directory containing the outputs for this field
    OUT = cf.PATH_OUT + fd + '/'
    #systematics maps directory
    PATH_SYST = OUT + 'systmaps/'
    #path in which the reformatted systematics maps will be placed
    PATH_TX = PATH_SYST + 'for_txpipe/'
    if not os.path.exists(PATH_TX):
        os.system(f'mkdir -p {PATH_TX}')

    #systematics maps
    maps = sorted(glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp'))
    
    #cycle through the maps
    for m in maps:
        #load the HealSparse map
        hsp_map = hsp.HealSparseMap.read(m)
        #retrieve the basename of the file
        bn = os.path.basename(m)
        #see if the input map contains a recarray
        if len(hsp_map.dtype) > 0:
            for n in hsp_map.dtype.names:
                fname = f'{PATH_TX}{bn[:-4]}_{n}.fits'
                healsparseToFITS(hsp_map[n], fname, nest=False)
        else:
            fname = f'{PATH_TX}{bn[:-4]}.fits'
            healsparseToFITS(hsp_map, fname, nest=False)


    