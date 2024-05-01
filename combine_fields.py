#####################################################################################################
# - Combines maps from all individual fields specified in the config file into single maps (i.e. 
#   one map per quantity across all fields simultaneously).
#####################################################################################################

import os
import config
import healsparse
import glob


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
systs = sorted(
    [os.path.basename(m) for m in glob.glob(f'{PATH_SYST}*_{cf.nside_hi}.hsp') 
                                   + glob.glob(f'{PATH_SYST}*_{cf.nside_hi}_*.hsp')]
)


