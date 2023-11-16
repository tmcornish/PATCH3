#####################################################################################################
# - Applies various cuts to each of the catalogues downloaded from the HSC database. 
# - Combines all cleaned catalogues belonging to the same field.
#####################################################################################################

import os
import config
from astropy.table import Table, vstack
from astropy.io import fits
import numpy as np
import glob
from output_utils import colour_string, error_message

### SETTINGS ###
cf = config.cleanCats


###################
#### FUNCTIONS ####
###################

def basic_clean(t):
	'''
	Performs a basic clean of the input catalogue using pre-existing flags.

	Parameters
	----------
	t: astropy.table.Table
		Input catalogue.

	Returns
	-------
	t: astropy.table.Table
		Catalogue after applying the cuts.
	'''

	sel = np.ones(len(t), dtype=bool)
	#create an empty list to which column names with the 'is_null' suffix will be appended
	isnull_names = []

	#cycle through column names
	for key in t.colnames:
		if key.__contains__('isnull'):
			sel[t[key]] = 0
			isnull_names.append(key)
		else:
			#want to remove NaNs UNLESS they are photo-zs
			if not key.startswith('pz_'):
				sel[np.isnan(t[key])]

	t.remove_columns(isnull_names)
	t = t[sel]

	return t



def photom_cuts(t):
	'''
	Performs a photometric cuts on the catalogue.

	Parameters
	----------
	t: astropy.table.Table
		Input catalogue.

	Returns
	-------
	t: astropy.table.Table
		Catalogue after applying the cuts.
	'''

	#magnitude cut in primary band
	sel_maglim = np.ones(len(t), dtype=bool)
	maglim_mask = (t[f'{cf.band}_cmodel_mag'] - t[f'a_{cf.band}']) > cf.depth_cut
	sel_maglim[maglim_mask] = False

	#blending cut
	sel_blend = np.ones(len(t), dtype=bool)
	blend_mask = t[f'{cf.band}_blendedness_abs'] >= cf.blend_cut
	sel_blend[blend_mask] = False

	#S/N cut in primary band
	sel_sn_pri = np.ones(len(t), dtype=bool)
	sn_pri_mask = (t[f'{cf.band}_cmodel_flux'] / t[f'{cf.band}_cmodel_fluxerr']) < cf.sn_pri
	sel_sn_pri[sn_pri_mask] = False

	#S/N cut in all other bands
	sel_sn_sec = []
	for b in [x for x in cf.bands if x != cf.band]:
		sel_sn_now = np.ones(len(t), dtype=bool)
		sn_sec_mask = (t[f'{b}_cmodel_flux'] / t[f'{b}_cmodel_fluxerr']) < cf.sn_sec
		sel_sn_now[sn_sec_mask] = False
		sel_sn_sec.append(sel_sn_now)
	#ensure that source is above threshold in at least two of these bands
	sel_sn_sec = np.sum(sel_sn_sec, axis=0) >= 2

	#apply the cuts
	t = t[sel_maglim * sel_blend * sel_sn_pri * sel_sn_sec]

	return t
	

def gal_cut(t):
	'''
	Flags sources as stars or galaxies based on an 'extendedness' parameter.

	Parameters
	----------
	t: astropy.table.Table
		Input catalogue.

	Returns
	-------
	t_gals: astropy.table.Table
		Catalogue containing only the sources identified as galaxies.

	t_stars: astropy.table.Table
		Catalogue containing only the sources identified as stars.
	'''

	#identify stars via their 'extendedness' in the primary band
	star_mask = t[f'{cf.band}_extendedness_value'] == 0.

	#split the catalogue into stars and galaxies
	t_stars = t[star_mask]
	t_gals = t[~star_mask]

	return t_gals, t_stars


def write_output(t, fname):
	'''
	Flags sources as stars or galaxies based on an 'extendedness' parameter.

	Parameters
	----------
	t: astropy.table.Table
		Input catalogue.

	fname: str
		Filename to be given to the output file.
	'''

	#set up the header
	hdr = fits.Header()
	hdr['BAND'] = cf.band
	hdr['DEPTH'] = cf.depth_cut
	prm_hdu = fits.PrimaryHDU(header=hdr)

	#convert the catalogue into an HDU
	cat_hdu = fits.table_to_hdu(t)

	#write to file
	hdul = fits.HDUList([prm_hdu, cat_hdu])
	hdul.writeto(fname, overwrite=True)


#######################################################
###############    START OF SCRIPT    #################
#######################################################


#cycle through each of the fields
for fd in cf.fields:

	print(colour_string(fd, 'purple'))

	#create output directory for this field
	OUT = cf.PATH_OUT + fd
	print(OUT)
	if not os.path.exists(OUT):
		os.system(f'mkdir -p {OUT}')
	
	#see if the field has been split into multiple parts
	fname = f'{cf.PATH_DATA}{cf.prefix}{fd.upper()}{cf.suffix}.fits'
	print(fname)
	if os.path.exists(fname):
		data_all = Table.read(fname, format='fits')
		l_init = len(data_all)
		#apply basic clean
		data_all = basic_clean(data_all)
		#apply photometric cuts
		data_all = photom_cuts(data_all)
		l_final = len(data_all)

	else: 
		#see if catalogues exist for separate parts of the field
		parts = sorted(glob.glob(f'{cf.PATH_DATA}{cf.prefix}{fd.upper()}_part?{cf.suffix}.fits'))
		if len(parts) >= 1:
			#set up a list to contain data from all catalogues associated with this field
			data_all = []
			l_init = 0
			l_final = 0
			#cycle through each catalogue
			for cat in parts:
				data = Table.read(cat, format='fits')
				l_init += len(data)
				#apply basic clean
				data = basic_clean(data)
				#apply photometric cuts
				data = photom_cuts(data)
				data_all.append(data)
				l_final += len(data)
			#stack the data from each part
			data_all = vstack(data_all)

		else:
			error_message(cf.__name__, f'No catalogues found for field {fd.upper()}.')
			continue

	print(f'Began with {l_init} sources.')
	print(f'{l_final} sources remaining after cleaning.')

	#split catalogue into stars and galaxies
	data_gals, data_stars = gal_cut(data_all)
	print(f'{len(data_gals)} galaxies; {len(data_stars)} stars.')

	#write the catalogues to output files
	print('Writing outputs...')
	write_output(data_gals, f'{OUT}/{cf.cat_main}')
	write_output(data_stars, f'{OUT}/{cf.cat_stars}')






