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
from output_utils import colour_string, write_output_hdf, h5py_dataset_iterator
from gen_utils import error_message

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
			#want to remove NaNs UNLESS they are photo-zs or r/i-band magnitude corrections
			if not key.startswith('pz_') and not key.startswith('corr_'):
				sel[np.isnan(t[key])] = 0
	
	#also remove any sources that are not primary detections
	sel *= t['isprimary']

	#if told to apply 'main' and 'strict'  cuts at this stage, do so
	if len(cf.remove_if_flagged) > 0:
		import flags as fl
		to_remove = fl.combine_flags(
						t,
						fl.get_flags(
							cf.band,
							cf.secondary_bands(),
							cf.remove_if_flagged
						),
						combine_type='or'
		)
		sel *= ~to_remove

	t.remove_columns(isnull_names)
	t = t[sel]

	#correct the r/i-band magnitudes to using the appropriate corrections
	for b in ['r', 'i']:
		corr_mask = ~np.isnan(t[f'corr_{b}mag'])
		t[f'{b}_cmodel_mag'][corr_mask] += t[f'corr_{b}mag'][corr_mask]

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




def make_tomography_cat(z, fname):
	'''
	Makes an HDF file containing the information required by TXPipe's TXLensMaps stage.
	NOTE: This information includes data described as 'lens_weight', which is not relevant
	for this study. It is therefore assumed here that this quantity is 1 for all sources.

	Parameters
	----------
	z: array-like
		Redshifts for each object.

	fname: str
		Filename to be given to the output file. If the file already exists, will try
		to append to existing data in the file.
	'''

	#create a column for labeling sources according to their tomographic bin (set -1 by default)
	labels = np.full(len(z), -1, dtype='i1')
	#create an array to contain the counts in each bin
	counts = np.zeros(len(cf.zbins)-1, dtype='i8')


	#cycle through the redshift bins
	for i in range(len(cf.zbins)-1):
		zmask = (z >= cf.zbins[i]) * (z < cf.zbins[i+1])
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
		g.attrs['nbin'] = len(cf.zbins) - 1
		for d, n in zip(data, names):
			dset = g.create_dataset(n, data=d, shape=(len(d),), dtype=d.dtype)


def flag_stars(t):
	'''
	Adds a column to the provided table flagging all stars.

	Parameters
	----------
	t: astropy.table.Table
		Input catalogue.

	'''

	#identify stars via their 'extendedness' in the primary band
	star_mask = t[f'{cf.band}_extendedness_value'] == 0.
	#add this as a column to the Table
	t['is_star'] = star_mask

#######################################################
###############    START OF SCRIPT    #################
#######################################################


#get a dictionary of all fields being analysed and their respective subfields
f_in_g = cf.fields_in_global()

#cycle through each global field
for g in f_in_g:
	print(colour_string(g.upper(), 'orange'))
	#cycle through each of the subfields
	for fd in f_in_g[g]:

		print(colour_string(fd, 'purple'))

		#create output directory for this field
		OUT = f'{cf.PATH_OUT}{g}/{fd}'
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
			print('Applying basic clean...')
			data_all = basic_clean(data_all)
			flag_stars(data_all)
			l_bc = len(data_all)
			write_output_hdf(data_all, hdf_basic, mode=mode, group='photometry')
			#apply photometric cuts and write to HDF file
			print('Applying photometric cuts...')
			data_all = photom_cuts(data_all)
			write_output_hdf(data_all, hdf_full, mode=mode, group='photometry')
			l_final = len(data_all)
		else: 
			#see if catalogues exist for separate parts of the field
			parts = sorted(glob.glob(f'{cf.PATH_DATA}{cf.prefix}{fd.upper()}_part??{cf.suffix}.fits'))
			if len(parts) >= 1:
				#set up a list to contain data from all catalogues associated with this field
				data_all = []
				l_init = 0		#counter for number of sources in raw data
				l_bc = 0		#counter for number of sources in basic-cleaned data
				l_final = 0		#counter for number of sources in final catalogue
				#cycle through each catalogue
				for i,cat in enumerate(parts):
					print(f'Cleaning part {i+1}...')
					data = Table.read(cat, format='fits')
					l_init += len(data)
					#apply basic clean and write to HDF file
					print('Applying basic clean...')
					data = basic_clean(data)
					flag_stars(data)
					l_bc += len(data)
					write_output_hdf(data, hdf_basic, mode=mode, group='photometry')
					#apply photometric cuts and write to HDF file
					print('Applying photometric cuts...')
					data = photom_cuts(data)
					write_output_hdf(data, hdf_full, mode=mode, group='photometry')
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

		print(colour_string(f'Began with {l_init} sources.', 'green'))
		print(colour_string(f'{l_bc} remained after basic cleaning.', 'green'))
		print(colour_string(f'{l_final} sources remaining after full cleaning.', 'green'))

		#split catalogue into stars and galaxies
		data_gals, data_stars = gal_cut(data_all)
		print(colour_string(f'{len(data_gals)} galaxies; {len(data_stars)} stars.', 'green'))

		#write the catalogues to output files
		print('Writing outputs...')
		hdf_stars = f'{OUT}/{cf.cat_stars}'
		write_output_hdf(data_gals, hdf_full, mode='w', group='photometry')
		write_output_hdf(data_stars, hdf_stars, mode='w', group='photometry')

		#also produce a tomgraphy catalogue
		hdf_tomo = f'{OUT}/{cf.cat_tomo}'
		make_tomography_cat(data_gals[cf.zcol], hdf_tomo)




print('Consolidating catalogues from subfields...')
cats = [cf.cat_basic, cf.cat_main, cf.cat_stars, cf.cat_tomo]
for g in f_in_g:
	print(colour_string(g.upper(), 'orange'))
	#cycle through the catalogue types
	for cat in cats:
		print(colour_string(cat, 'cyan')) 
		fname = f'{cf.PATH_OUT}{g}/{cat}'
		with h5py.File(fname, 'w') as fmain:
			for fd in f_in_g[g]:
				print(f'Adding data from subfield {fd}...')
				cat_now = f'{cf.PATH_OUT}{g}/{fd}/{cat}'
				with h5py.File(cat_now, 'r') as fnow:
					#check the current structure of the main file
					paths_current = [p for p,_ in h5py_dataset_iterator(fmain)]
					#iterate through each branch of the file tree
					for (path, dset) in h5py_dataset_iterator(fnow):
						N = len(dset)
						dt = dset.dtype
						if path in paths_current:
							dset_main = fmain[path]
							dset_main.resize((len(dset_main)+N,))
							dset_main[-N:] = dset[:]
						else:
							dset = fmain.create_dataset(path, shape=(N,), data=dset[:], maxshape=(None,), dtype=dt)
			#if the tomography catalogue, need to ensure that the nbin attribute is available for TXPipe
			if cat == cf.cat_tomo:
				fmain['tomography'].attrs['nbin'] = len(cf.zbins) - 1





