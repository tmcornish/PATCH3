#####################################################################################################
# Config file for the pHSC3 pipeline.
#####################################################################################################


################
#### GLOBAL ####
################

class cf_global:

	#whether to run on glamdring or locally
	LOCAL = True

	#relevant directories (dependent on whether being run locally or on glamdring)
	if LOCAL:
		PATH_PIPE = '/home/cornisht/LSST_clustering/pHSC3/'			#path to pipeline scripts
		PATH_DATA = '/home/cornisht/LSST_clustering/Data/HSC_DR3/'	#where the raw data are stored
		PATH_OUT =  f'{PATH_PIPE}out/'								#main directory for outputs
	else:
		PATH_PIPE = '/mnt/zfsusers/tcornish/pHSC3/'					#path to pipeline scripts
		PATH_DATA = '/mnt/extraspace/tmcornish/Datasets/HSC_DR3/'	#where the raw data are stored
		PATH_OUT = '/mnt/extraspace/tmcornish/pHSC3_out/'			#main directory for outputs
	#directory for figures
	PATH_PLOTS = PATH_OUT + 'figures/'
	

	#data release
	dr = 'pdr3_wide'

	#lists detailing which sub-fields belong to which equatorial field
	equatora = [f'equator{i:02d}' for i in [21,22,23,0,1,2]]
	equatorb = [f'equator{i:02d}' for i in [8,9,10,11,12,13,14,15]]
	
	#fields for which the pipeline is to be run
	#fields = ['aegis']
	fields = ['hectomap']
	#fields = equatora
	#fields = equatorb
	#fields = ['aegis', 'hectomap'] + equatora + equatorb

	#file containing the metadata
	metafile = f'{PATH_DATA}PDR3_WIDE_frames.fits'
	#basename to be given to the split metadata files
	metasplit = 'frames_metadata.fits'

	#main photometric band
	band = 'i'
	#list of all photometric bands
	bands = ['g', 'r', 'i', 'z', 'y']
	#dictionary of alternative names for certain bands
	bands_alt = {
		'i' : ['i2'],
		'r' : ['r2']
	}

	#S/N thresholds in primary band and other bands
	sn_pri = 10.
	sn_sec = 5.
	#depth limit in primary band
	depth_cut = 24.5

	#names to be given to the relevant catalogues
	cat_basic = 'basicclean_catalogue.hdf5'
	cat_main = 'clean_catalogue.hdf5'
	cat_stars = 'star_catalogue.hdf5'
	cat_tomo = 'tomography_catalogue.hdf5'

	#NSIDE parameter for the low- and high-resolution components of the maps
	nside_lo = 32
	nside_hi = 512
	#low-resolution NSIDE parameter to use for splitting the data
	nside_cover = 8
	
	#threshold below which pixels in the survey mask will be considered masked
	weight_thresh = 0.5

	#basenames for the various maps
	dustmaps = 'dustmaps.hsp'
	bo_mask = 'bo_mask.hsp'
	masked_frac = 'masked_fraction.hsp'
	survey_mask = 'survey_mask.hsp'
	star_map = 'star_counts.hsp'
	depth_map = 'depth_map.hsp'

	#basenames for the count and density maps
	ngal_maps = 'ngal_maps.hsp'
	deltag_maps = 'deltag_maps.hsp'

	#redshift column to use for tomography
	zcol = 'pz_best_dnnz'
	#redshift bin edges
	zbins = [0.3, 0.6, 0.9, 1.2, 1.5]


	@classmethod
	def get_global_fields(cls):
		fields_global = []
		if 'hectomap' in cls.fields:
			fields_global.append('hectomap')
		if 'aegis' in cls.fields:
			fields_global.append('aegis')
		if any(x in cls.equatora for x in cls.fields):
			fields_global.append('equatora')
		if any(x in cls.equatorb for x in cls.fields):
			fields_global.append('equatorb')

		return fields_global


##################
#### get_data ####
##################

class getData(cf_global):

	name = 'getData'

	#submit/download data requests
	submit = True
	download = True
	#include photo-z information
	photoz = True
	#apply cuts based on existing flags in the catalogues
	strict_cuts = True

	#maximum number of sources allowed before catalogues split for field
	Nmax = 20_000_000

	@classmethod
	def suffix(cls):
		suff = '_forced'
		if not cls.strict_cuts:
			suff += '_shearcat'
		return suff


########################
#### split_metadata ####
########################

class splitMetadata(cf_global):

	name = 'splitMetadata'

	#boundaries of each global field, ordered [RA_min, RA_max, DEC_min, DEC_max]
	bounds = {
		'aegis' : [212., 216., 51.6, 53.6],
		'equatora' : [326.25, 41.25, -8., 8.],
		'equatorb' : [125., 227.5, -4., 7.],
		'hectomap' : [195., 255., 41.5, 45.]
	}

	#whether to further split the metadata by filter (can help Decasu avoid memory issues)
	split_by_band = True


##########################
#### clean_catalogues ####
##########################

class cleanCats(cf_global):

	name = 'cleanCats'

	#parts of the naming structure for raw data files
	prefix = f'{cf_global.dr.upper()}_'
	suffix = getData.suffix()

	#blending cut (maximum allowed flux estimated to be from blending)
	blend_cut = 10. ** (-0.375)

	@classmethod
	def fields_in_global(cls):
		fields_global = cls.get_global_fields()
		f_in_g = {}
		for g in fields_global:
			if g == 'equatora':
				f_in_g[g] = [f for f in cls.equatora if f in cls.fields]
			elif g == 'equatorb':
				f_in_g[g] = [f for f in cls.equatorb if f in cls.fields]
			else:
				f_in_g[g] = [g]
		return f_in_g



##################################
#### make_maps_from_catalogue ####
##################################

class makeMapsFromCat(cf_global):

	name = 'makeMapsFromCat'

	#NSIDE for the upgraded-resolution version of the bright object mask
	nside_mask = 16384

	#column names for flags identifying sources near bright objects
	bo_flags = [f'{cf_global.band}_mask_brightstar_ghost15',
				f'{cf_global.band}_mask_brightstar_halo',
				f'{cf_global.band}_mask_brightstar_blooming']


#################################
#### make_maps_from_metadata ####
#################################

class makeMapsFromMetadata(splitMetadata):

	name = 'makeMapsFromMetadata'

	#decasu config file
	configfile = f'{cf_global.PATH_PIPE}decasu_config_hpix_hsc_dr3.yaml'
	#number of cores to use if running locally (if running on glamdring need to specify that elsewhere)
	ncores = 18
	


##########################
#### make_galaxy_maps ####
##########################

class makeGalaxyMaps(cf_global):

	name = 'makeGalaxyMaps'


###############################
#### compute_power_spectra ####
###############################

class computePowerSpectra(cf_global):

	name = 'computePowerSpectra'

	#systematics maps to deproject
	systs = [
		]
	
	#TODO: add option to incorporate all decasu output maps in systs
	
	#output file for power spectrum information
	outfile = f'power_spectra_info_{cf_global.nside_hi}.hdf5'

	

