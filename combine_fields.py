#####################################################################################################
# - Combines maps from all individual fields specified in the config file into single maps (i.e. 
#   one map per quantity across all fields simultaneously).
#####################################################################################################

import os
import sys
from configuration import PipelineConfig as PC
import healsparse as hsp
import glob
import map_utils as mu
from functools import reduce
import numpy as np

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='combineFields')

#############################
######### FUNCTIONS #########
#############################

def combine_maps(map_name, fields):
    '''
    Given the base name for a map, will load that map for each of the
    specified fields (if it exists) and return the combination of those
    maps.

    Parameters
    ----------
    map_name: str
        The base name of the maps to be loaded.
    
    fields: list
        A list of the fields for which the maps are to be loaded.
    
    Returns
    -------
    union: HealSparseMap
        Combination of the loaded HealSparse maps.
    '''

    #load the maps
    maps = [hsp.HealSparseMap.read(f'{cf.paths.out}{fd}/{map_name}') 
            for fd in fields
            if os.path.exists(f'{cf.paths.out}{fd}/{map_name}')
            ]
    #if map only exists for one field, return that as is
    if len(maps) == 1:
        return maps[0]
    #if no maps are found, return None
    elif len(maps) == 0:
        return None
    
    #check if these maps are single maps or recarrays of multiple maps
    ndtype = len(maps[0].dtype)
    if ndtype > 0:
        #cycle through the individual maps and combine from each field
        unions_indiv = [hsp.operations.sum_union([m[m.dtype.names[i]] for m in maps]) for i in range(ndtype)]
        #initialise a recarray for the output
        union, _, px_data_u = mu.initialiseRecMap(
            unions_indiv[0].nside_coverage,
            unions_indiv[0].nside_sparse,
            maps[0].dtype.names,
            pixels=unions_indiv[0].valid_pixels,
            dtypes=str(maps[0].dtype[0])
            )
        #populate the maps
        for i in range(ndtype):
            union[union.dtype.names[i]].update_values_pix(
                px_data_u,
                unions_indiv[i][px_data_u]
                )
    else:
        union = hsp.operations.sum_union(maps)
    return union


#######################################################
###############    START OF SCRIPT    #################
#######################################################


#make a directory for the combined maps if necessary
OUT_MAIN = cf.paths.out + 'combined/'
OUT_SYSTS = OUT_MAIN + 'systmaps/'
if not os.path.exists(OUT_SYSTS):
    os.system(f'mkdir -p {OUT_SYSTS}')

#list of the fields being analysed
fields = cf.fields

#lists of quantities and systematics that have been mapped at the desired resolution
add = lambda x,y : x + y
quants = [[os.path.basename(m) for m in glob.glob(f'{cf.paths.out}{fd}/*_{cf.nside_hi}.hsp') 
                                   + glob.glob(f'{cf.paths.out}{fd}/*_{cf.nside_hi}_*.hsp')]
            + ['systmaps/'+os.path.basename(m) for m in glob.glob(f'{cf.paths.out}{fd}/systmaps/*_{cf.nside_hi}.hsp') 
                                   + glob.glob(f'{cf.paths.out}{fd}/systmaps/*_{cf.nside_hi}_*.hsp')]
            for fd in fields
            ]
quants = sorted(np.unique(reduce(add, quants)))

for q in quants:
    #load the HealSparse maps from each field for this quantity and combine them
    union = combine_maps(q, fields)
    #save the map
    if union is not None:
        union.write(f'{OUT_MAIN}{q}', clobber=True)
