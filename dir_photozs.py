#####################################################################################################
# - Estimates the n(z) distributions for each tomographic bin using the direct reweighted estimation
#   method (DIR). Uses many-band photo-zs from COSMOS as the 'true' valus for calibration.
# - NOTE: Unlike other stages, this script has no dependence on the HSC field being analysed, as the
#   n(z) distributions are computed using specific catalogues.
# - TODO: 
#	- Adapt script so that the data can be automatically retrieved if they do not exist?
#####################################################################################################

import os
import sys
from configuration import PipelineConfig as PC
import h5py
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack
from sklearn.neighbors import NearestNeighbors
import scipy.spatial as spatial

### SETTINGS ###
config_file = sys.argv[1]
cf = PC(config_file, stage='dirPhotozs')


#######################################################
###############    START OF SCRIPT    #################
#######################################################

print('Loading data...')
########################
# see if a catalogue containing HSC phtoometry in the COSMOS field exists
if os.path.exists(cf.hsc_cosmos_cat):
	#load the data and convert to an astropy Table
	with h5py.File(cf.hsc_cosmos_cat, 'r') as hf:
		gp = hf['photometry']
		#columns to include in the Table
		cols_hsc = ['ra', 'dec']
		cols_hsc += [f'{b}_cmodel_mag' for b in cf.bands.all]
		cols_hsc += [f'pz_{s}_dnnz' for s in ['mean', 'mode', 'best', 'mc', 'err68_min', 'err68_max', 'err95_min', 'err95_max']]
		hsc_data = Table([gp[col][:] for col in cols_hsc], names=cols_hsc)
else:
	#tell user how to create the required catalogue and exit
	raise FileNotFoundError(
		f'''
		File {cf.hsc_cosmos_cat} not found.
		Download the HSC data using the SQL query in data_query/run_COSMOS_DUD.sql 
		then use the cleanCatalogues stage of the pipeline to create a cleaned 
		version. 
		hsc_cosmos_cat should then be set to the path of this cleaned version.
		'''
	)

# see if the COSMOS2020	catalogue exists at the specified destination
if os.path.exists(cf.cosmos_cat):
	#load the catalogue as an astropy Table
	cosmos_data = Table.read(cf.cosmos_cat)
else:
	raise FileNotFoundError(
		f'''
		File {cf.cosmos_cat} not found.
		Download the COSMOS2020 catalogue by running
		'wget https://cdsarc.cds.unistra.fr/ftp/J/ApJS/258/11/fits/COSMOS2020_CLASSIC_R1_v2.1_p3.fits'
		cosmos_cat should then be set to the path of this file.
		'''
	)

#apply cuts to the COSMOS data
cosmos_mask = (cosmos_data['lp_zMinChi2'] > 0.01) & (cosmos_data['lp_zMinChi2'] < 9) & \
			(cosmos_data['lp_type'] == 0.) & \
			(cosmos_data['lp_mass_best'] > 7.5) & \
			(np.maximum(cosmos_data['lp_zPDF_u68'] - cosmos_data['lp_zPDF'], cosmos_data['lp_zPDF'] - cosmos_data['lp_zPDF_l68']) < 0.05 * (1. + cosmos_data['lp_zMinChi2'])) & \
			(cosmos_data['lp_chi2_best'] < cosmos_data['lp_chis']) & \
			(cosmos_data['lp_chi2_best'] / cosmos_data['lp_NbFilt'] < 5) & \
			(np.isnan(cosmos_data['lp_zp_2']))
cosmos_data = cosmos_data[cosmos_mask]

print('Cross-matching HSC and COSMOS data...')
##############################################
hsc_skycoord = SkyCoord(hsc_data['ra']*u.deg, hsc_data['dec']*u.deg)
cosmos_skycoord = SkyCoord(cosmos_data['ALPHA_J2000'], cosmos_data['DELTA_J2000'])
#get nearest neighbours on the celestial sphere
cosmos_index, dist_2d, _ = hsc_skycoord.match_to_catalog_sky(cosmos_skycoord)
#identify everything within the tolerance
within_tol = dist_2d.degree*3600 < cf.cross_tol
cosmos_index_matched = cosmos_index[within_tol]
#mask the Tables and join
hsc_matched = hsc_data[within_tol]
cosmos_matched = cosmos_data[cosmos_index_matched]
cat_matched = hstack([hsc_matched, cosmos_matched])
#clear some memory
del cosmos_data, hsc_matched, cosmos_matched

print('Computing colour-space weights...')
##########################################
train_sample = np.array([np.array(cat_matched[f'{b}_cmodel_mag']) for b in cf.bands.all]).T
train_z = np.array(cat_matched['lp_zMinChi2'])
photoz_sample = np.array([hsc_data[f'{b}_cmodel_mag'] for b in cf.bands]).T
#find k nearest neighbours in 5D colour space
n_nbrs = NearestNeighbors(n_neighbors=cf.kNN, algorithm='kd_tree', metric='euclidean').fit(train_sample)
distances, _ = n_nbrs.kneighbors(train_sample)
#find the distance to the kth-nearest neighbour
distances = np.amax(distances, axis=1)
#find all sources in the full HSC data within the distance to the 20th nearest neighbour
tree_NN_lookup = spatial.cKDTree(photoz_sample, leafsize=cf.kNN*2)
num_photoz=np.array([len(tree_NN_lookup.query_ball_point(t,d+1E-6)) 
					for t,d in zip(train_sample,distances)])
#compute weights as the ratio of the number of photo-z neighbours (N_HSC=num_photoz) to the 
#number of COSMOS neighbours (N_COSMOS=kNN), normalised by the number of photo-z (HSC) objects
weights = (num_photoz / 20) * (len(train_sample) / len(photoz_sample))

print('Computing n(z) for each tomographic bin...')
###################################################
zbin_masks = cf.get_samples(cat_matched)
#bin edges and centres for the n(z) distributions
bins = np.arange(0., 7.+cf.dz, cf.dz) 
bin_centres = (bins[1:] + bins[:-1]) / 2
#calculate and ave outputs to file
with h5py.File(cf.nofz_files.nz_dir, 'w') as hf:
	hf.create_dataset('z', data=bin_centres)
	for i,zb in enumerate(zbin_masks):
		nz, _ = np.histogram(cat_matched['lp_zMinChi2'][zbin_masks[zb]], bins, weights=weights[zbin_masks[zb]], density=True)
		hf.create_dataset(f'nz_{i}', data=nz)
