#####################################################################################################
# - Applies various cuts to each of the catalogues downloaded from the HSC database. 
# - Combines all cleaned catalogues belonging to the same field.
#####################################################################################################

import os
import config
from astropy.table import Table, vstack, Column
from astropy.io import fits
import numpy as np
import glob
import h5py
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
	Writes an astropy Table to a file (type inferred from filename) with appropriate metadata.

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


def write_output_hdf(t, fname, colnames=None, group=None, mode='a'):
	'''
	Writes an astropy Table to a hdf5 file.

	Parameters
	----------
	t: astropy.table.Table
		Input catalogue.

	fname: str
		Filename to be given to the output file. If the file already exists, will try
		to append to existing data in the file.

	colnames: list or None
		List of columns to be included in the output file. If None, writes all columns
		to the output.

	group: str
		Name of the group to which the data will be written. If None, will write the data
		to the 'root' of the hierarchy in the output file.

	mode: str
		Mode in which to open the HDF file (e.g. 'w' for 'write', 'r' for 'read', etc.)
	'''

	#if colnames=None, set equal to the list of column names in t
	if colnames is None:
		colnames = t.colnames
	#if group=None, assign an empty string
	if group is None:
		group = ''

	#get the length of the Table
	N = len(t)
	#open the file
	with h5py.File(fname, mode) as hf:
		#cycle through the columns
		for col in colnames:
			#get the dtype of the column data
			dt = t[col].dtype
			#see if the current column already exists as a DataSet
			if col in hf[f'/{group}'].keys():
				#reshape the existing data and fill the empty entries with the new data
				dset = hf[f'/{group}/{col}']
				dset.resize((len(dset)+N,))
				dset[-N:] = t[col]
			else:
				#create the dataset if it doesn't exist
				dset = hf.create_dataset(f'{group}/{col}', shape=(N,), data=t[col], maxshape=(None,), dtype=dt)


def make_tomography_cat(z, zbins, fname):
	'''
	Makes an HDF file containing the information required by TXPipe's TXLensMaps stage.
	NOTE: This information includes data described as 'lens_weight', which is not relevant
	for this study. It is therefore assumed here that this quantity is 1 for all sources.

	Parameters
	----------
	z: array-like
		Redshifts for each object.

	zbins: array-like
		Edges of each redshift bin. Bins are left-inclusive and right-exclusive.

	fname: str
		Filename to be given to the output file. If the file already exists, will try
		to append to existing data in the file.
	'''

	#create a column for labeling sources according to their tomographic bin (set -1 by default)
	labels = np.full(len(z), -1, dtype='i1')
	#create an array to contain the counts in each bin
	counts = np.zeros(len(zbins)-1, dtype='i8')


	#cycle through the redshift bins
	for i in range(len(zbins)-1):
		zmask = (z >= zbins[i]) * (z < zbins[i+1])
		labels[zmask] = i
		counts[i] = zmask.sum()

	#create an array for containing the total counts in all bins
	counts_2d = np.array([counts.sum()])

	#create a column for 'lens_weight' and simply assign a value of 1 for every source
	lens_weight = np.ones(len(z), dtype='f8')

	#compile the relevant data into a list and assign them names
	data = [labels, counts, counts_2d, lens_weight]
	names = ['bin', 'counts', 'counts_2d', 'lens_weight']

	#write to the output file
	with h5py.File(fname, 'w') as hf:
		g = hf.create_group('tomography')
		for d, n in zip(data, names):
			dset = g.create_dataset(n, data=d, shape=(len(d),), dtype=d.dtype)




#######################################################
###############    START OF SCRIPT    #################
#######################################################


#cycle through each of the fields
for fd in cf.fields:

	print(colour_string(fd, 'purple'))

	#create output directory for this field
	OUT = cf.PATH_OUT + fd
	print(f'Output directory: {OUT}')
	if not os.path.exists(OUT):
		os.system(f'mkdir -p {OUT}')
	
	#filenames to be given to the HDF format output files
	hdf_basic = f'{OUT}/{cf.cat_basic}'
	hdf_full = f'{OUT}/{cf.cat_main}'
	#see if the field has been split into multiple parts
	fname = f'{cf.PATH_DATA}{cf.prefix}{fd.upper()}{cf.suffix}.fits'
	#initially enable 'write' mode for output files
	mode = 'w'
	if os.path.exists(fname):
		data_all = Table.read(fname, format='fits')
		l_init = len(data_all)
		#apply basic clean and write to HDF file
		data_all = basic_clean(data_all)
		write_output_hdf(data_all, hdf_basic)
		#apply photometric cuts and write to HDF file
		data_all = photom_cuts(data_all)
		write_output_hdf(data_all, hdf_full)
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
				#apply basic clean and write to HDF file
				data = basic_clean(data)
				write_output_hdf(data, hdf_basic, mode=mode)
				#apply photometric cuts and write to HDF file
				data = photom_cuts(data)
				write_output_hdf(data, hdf_full, mode=mode)
				data_all.append(data)
				l_final += len(data)
				try: 
					assert mode == 'a'
				except AssertionError:
					mode = 'a'
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
	#print('Writing outputs...')
	#write_output(data_gals, f'{OUT}/{cf.cat_main}')
	#write_output(data_stars, f'{OUT}/{cf.cat_stars}')

	#also produce a tomgraphy catalogue
	hdf_tomo = f'{OUT}/{cf.cat_tomo}'
	make_tomography_cat(data_gals[cf.zcol], cf.zbins, hdf_tomo)






