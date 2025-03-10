##############################################################################
# - Combines maps from all individual fields specified in the config file into
#   single maps (i.e. one map per quantity across all fields simultaneously).
##############################################################################

import os
import sys
from configuration import PipelineConfig as PC
import healsparse as hsp
import glob
import map_utils as mu
from functools import reduce
import numpy as np

# SETTINGS #
config_file = sys.argv[1]
cf = PC(config_file, stage='combineFields')

#############################
#         FUNCTIONS         #
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

    # Load the maps
    maps = [hsp.HealSparseMap.read(f'{cf.paths.out}{fd}/{map_name}')
            for fd in fields
            if os.path.exists(f'{cf.paths.out}{fd}/{map_name}')
            ]
    # If map only exists for one field, return that as is
    if len(maps) == 1:
        return maps[0]
    # If no maps are found, return None
    elif len(maps) == 0:
        return None

    # Check if these maps are single maps or recarrays of multiple maps
    ndtype = len(maps[0].dtype)
    if ndtype > 0:
        # Cycle through the individual maps and combine from each field
        unions_indiv = [hsp.operations.sum_union([m[m.dtype.names[i]]
                                                  for m in maps])
                        for i in range(ndtype)]
        # Initialise a recarray for the output
        union, _, px_data_u = mu.initialiseRecMap(
            unions_indiv[0].nside_coverage,
            unions_indiv[0].nside_sparse,
            maps[0].dtype.names,
            pixels=unions_indiv[0].valid_pixels,
            dtypes=str(maps[0].dtype[0])
            )
        # Populate the maps
        for i in range(ndtype):
            union[union.dtype.names[i]].update_values_pix(
                px_data_u,
                unions_indiv[i][px_data_u]
                )
    else:
        union = hsp.operations.sum_union(maps)
    return union


def add(x, y):
    '''
    Simply takes two inputs and adds them together. For use with e.g.
    functools.reduce to aid in adding long sequences together.

    Parameters
    ----------
    x,y : int, float, or array-like
        The quantities to be added together.

    Returns
    -------
    z: int, float or array-like
        The sum of the two inputs.
    '''
    z = x + y
    return z

#######################################################
#                  START OF SCRIPT                    #
#######################################################


# Make a directory for the combined maps if necessary
OUT_MAIN = cf.paths.out + 'combined/'
OUT_SYSTS = OUT_MAIN + 'systmaps/'
if not os.path.exists(OUT_SYSTS):
    os.system(f'mkdir -p {OUT_SYSTS}')

# List of the fields being analysed
fields = cf.fields

# Lists of quantities and systematics mapped at the desired resolution
quants = [[os.path.basename(m)
           for m in glob.glob(f'{cf.paths.out}{fd}/*_nside{cf.nside_hi}*.hsp')]
          + ['systmaps/'+os.path.basename(m)
             for m in glob.glob(f'{cf.paths.out}{fd}/systmaps/'
                                f'*_nside{cf.nside_hi}*.hsp')]
          for fd in fields]
quants = sorted(np.unique(reduce(add, quants)))

for q in quants:
    # Load HealSparse maps from each field for this quantity and combine them
    union = combine_maps(q, fields)
    # Save the map
    if union is not None:
        union.write(f'{OUT_MAIN}{q}', clobber=True)
