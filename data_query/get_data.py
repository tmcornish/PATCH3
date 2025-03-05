from pdr3_queries import write_frames, write_fieldsearch
import os
import sys
from astropy.table import Table
import numpy as np
from configuration import PipelineConfig as PC
from output_utils import colour_string

# SETTINGS #
config_file = sys.argv[1]
cf = PC(config_file, stage='getData')

################################
#                              #
#  Download catalog-level data #
#                              #
################################

# Directory in which the SQL queries will be saved
PATH_SQL = 'sql_files/'
if not os.path.exists(PATH_SQL):
    os.system(f'mkdir -p {PATH_SQL}')

# Per-frame metadata
write_frames(cf.dr, f'{PATH_SQL}frames_wide.sql', submit=cf.submit,
             dir_out=cf.paths.data, do_download=cf.download)


# Retrieve the field names and boundary coordinates
bounds_file = './field_counts_and_boundaries.fits'
t_bounds = Table.read(bounds_file).to_pandas().set_index('field')

# WIDE fields
for fd in cf.get_subfields():
    print(colour_string(fd, 'purple'))
    fld = fd.encode('utf-8')
    # Get the number of objects
    Nobj = t_bounds.loc[fld]['countof']
    # Get the min and max Dec.
    dec_min = t_bounds.loc[fld]['dec_min']
    dec_max = t_bounds.loc[fld]['dec_max']

    # Calculate how many parts the catalogue should be split into
    nparts = int(np.ceil(Nobj / cf.Nmax))

    if nparts > 1:
        # Calculate the width of each part in declination
        ddec = (dec_max - dec_min) / nparts
        for i in range(nparts):
            # Determine the min and max Dec. for the current part
            decmin_now = dec_min + (i * ddec)
            decmax_now = decmin_now + ddec

            # Write and submit a query for the data in this sub-field
            write_fieldsearch(cf.dr, fd,
                              f'{PATH_SQL}run_{fd.upper()}_part{i+1:02d}.sql',
                              dir_out=cf.paths.data, do_photoz=cf.photoz,
                              mizuki_mstar_sfrs=cf.mizuki_mstar_sfrs,
                              submit=cf.submit, apply_cuts=cf.apply_cuts,
                              strict_cuts=cf.strict_cuts,
                              dec_range=[decmin_now, decmax_now],
                              do_download=cf.download, part=i+1)
    else:
        write_fieldsearch(cf.dr, fd, f'{PATH_SQL}run_{fd.upper()}.sql',
                          dir_out=cf.paths.data, do_photoz=cf.photoz,
                          mizuki_mstar_sfrs=cf.mizuki_mstar_sfrs,
                          submit=cf.submit, apply_cuts=cf.apply_cuts,
                          strict_cuts=cf.strict_cuts, do_download=cf.download)
