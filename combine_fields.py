#####################################################################################################
# - Combines maps from all individual fields specified in the config file into single maps (i.e. 
#   one map per quantity across all fields simultaneously).
#####################################################################################################

import os
import config
import healsparse as hsp
import glob
import map_utils as mu

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
    maps = [hsp.HealSparseMap.read(f'{cf.PATH_OUT}{fd}/{map_name}') for fd in fields]

    #check if these maps are single maps or recarrays of multiple maps
    ndtype = len(maps[0].dtype)
    if ndtype > 0:
        #cycle through the individual maps and combine from each field
        unions_indiv = [hsp.operations.sum_union([m[m.dtype.names[i]] for m in maps]) for i in range(ndtype)]
        #initialise a recarray for the output
        union, _, px_data_u = mu.initialiseRecMap(
            unions_indiv[0].nside_coverage,
            unions_indiv[0].nside_sparse,
            *unions_indiv[0].valid_pixels_pos(lonlat=True),
            maps[0].dtype.names,
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

cf = config.combineFields

#make a directory for the combined maps if necessary
OUT_MAIN = cf.PATH_OUT + 'combined/'
OUT_SYSTS = OUT_MAIN + 'systmaps/'
if not os.path.exists(OUT_SYSTS):
    os.system(f'mkdir -p {OUT_SYSTS}')

#list of the fields being analysed
fields = cf.get_global_fields()

#use the first field to identify the quantities that have been mapped (should be the same for all fields)
PATH_MAPS = cf.PATH_OUT + fields[0] + '/'
PATH_SYST = PATH_MAPS + 'systmaps/'
#lists of quantities and systematics that have been mapped at the desired resolution
quants = sorted(
    [os.path.basename(m) for m in glob.glob(f'{PATH_MAPS}*_{cf.nside_hi}.hsp') 
                                   + glob.glob(f'{PATH_MAPS}*_{cf.nside_hi}_*.hsp')]
)
quants.extend(
    sorted(
        [
            'systmaps/'+os.path.basename(m) for m in glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') 
                                   + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp')
        ]
    )
)

for q in quants:
    #load the HealSparse maps from each field for this quantity and combine them
    union = combine_maps(q, fields)
    #save the map
    union.write(f'{OUT_MAIN}{q}', clobber=True)
