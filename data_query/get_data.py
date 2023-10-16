from pdr3_queries import write_frames, write_fieldsearch
import os
import contextlib
from astropy.table import Table
import numpy as np

################################
#                              #
#  Download catalog-level data #
#                              #
################################

#directory in which to save the data
PATH_DATA = '/mnt/extraspace/tmcornish/Datasets/HSC_DR3/'
#directory in which the SQL queries will be saved
PATH_SQL = 'sql_files/'
if not os.path.exists(PATH_SQL):
	os.system(f'mkdir -p {PATH_SQL}')


#Per-frame metadata
#write_frames('pdr3_wide', f'{PATH_SQL}frames_wide.sql', submit=True, dir_out=PATH_DATA, do_download=True)



#retrieve the field names and boundary coordinates
bounds_file = './field_counts_and_boundaries.fits'
t_bounds = Table.read(bounds_file).to_pandas().set_index('field')

#max number of sources allowed before splitting catalogues
Nmax = 20_000_000

#WIDE fields
for fld in t_bounds.index:
	fd = fld.decode('utf-8')
	#get the number of objects
	Nobj = t_bounds.loc[fld]['countof']
	#get the min and max Dec.
	dec_min = t_bounds.loc[fld]['dec_min']
	dec_max = t_bounds.loc[fld]['dec_max']

	#calculate how many parts the catalogue should be split into
	nparts = int(np.ceil(Nobj / Nmax))


	if nparts > 1:
		#calculate the width of each part in declination
		ddec = (dec_max - dec_min) / nparts
		for i in range(nparts):
			#determine the min and max Dec. for the current part
			decmin_now = dec_min + (i * ddec)
			decmax_now = decmin_now + ddec
			
			#write and submit a query for the data in this sub-field
			write_fieldsearch('pdr3_wide', fd, f'{PATH_SQL}run_{fd}_part{i+1}.sql', dir_out=PATH_DATA, do_photoz=True, submit=True, strict_cuts=True, dec_range=[decmin_now,decmax_now], do_download=True, part=i+1)
	else:
		write_fieldsearch('pdr3_wide', fd, f'{PATH_SQL}run_{fd}.sql', dir_out=PATH_DATA, do_photoz=True, submit=True, strict_cuts=True, do_download=True)




	


