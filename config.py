#####################################################################################################
# Config file for the pHSC3 pipeline.
#####################################################################################################


################
#### GLOBAL ####
################

class cf_global:

	LOCAL = True	#whether to run on locally or remotely
	NERSC = False	#if remotely: NERSC or glamdring

	#relevant directories (dependent on whether being run locally or on glamdring)
	if LOCAL:
		PATH_PIPE = '/home/cornisht/LSST_clustering/pHSC3/'			#path to pipeline scripts
		PATH_DATA = '/home/cornisht/LSST_clustering/Data/HSC_DR3/'	#where the raw data are stored
		PATH_OUT =  f'{PATH_PIPE}out/'								#main directory for outputs
	else:
		if NERSC:
			PATH_PIPE = '/globals/homes/t/tcornish/pHSC3/'			#path to pipeline scripts
			PATH_DATA = '/pscratch/sd/t/tcornish/HSC_DR3/'			#where the raw data are stored
			PATH_OUT = '/pscratch/sd/t/tcornish/pHSC3_out/'			#main directory for outputs
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
	#fields = ['combined']
	#fields = ['hectomap'] + equatora + equatorb
	#fields = ['combined', 'hectomap'] + equatora + equatorb


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
	nside_hi = 1024
	#low-resolution NSIDE parameter to use for splitting the data
	nside_cover = 8
	
	#threshold below which pixels in the survey mask will be considered masked
	weight_thresh = 0.5
	#whether to smooth the mask
	smooth_mask = True
	#width (in arcmin) of the smoothing kernel
	smooth_fwhm = 60.
	#threshold to apply to smoothed mask
	smooth_thresh = 0.4

	#function for generating filenames for band-specific dust extinction maps
	def dustmap_names(bands, nside_hi):
		return [f'dustmaps_{b}_{nside_hi}.hsp' for b in bands]
	#basenames for the various maps
	dustmaps = dustmap_names(bands, nside_hi)
	bo_mask = f'bo_mask_{nside_hi}.hsp'
	masked_frac = f'masked_fraction_{nside_hi}.hsp'
	survey_mask = f'survey_mask_{nside_hi}.hsp'
	star_map = f'star_counts_{nside_hi}.hsp'
	depth_map = f'depth_map_{nside_hi}.hsp'

	#basenames for the count and density maps
	ngal_maps = f'ngal_maps_{nside_hi}.hsp'
	deltag_maps = f'deltag_maps_{nside_hi}.hsp'

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
		if 'combined' in cls.fields:
			fields_global.append('combined')

		return fields_global


	@classmethod
	def get_bin_pairings(cls):
		'''
		Returns pairs of IDs for each tomographic bin being analysed. Also
		returns comma-separated string versions of the ID pairs.
		'''
		import itertools

		nbins = len(cls.zbins) - 1
		l = list(range(nbins))
		pairings = [i for i in itertools.product(l,l) if tuple(reversed(i)) >= i]
		pairings_s = [f'{p[0]},{p[1]}' for p in pairings]
		return pairings, pairings_s


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
	apply_cuts = False
	strict_cuts = False

	#maximum number of sources allowed before catalogues split for field
	Nmax = 5_000_000

	@classmethod
	def suffix(cls):
		if cls.strict_cuts:
			suff = ''
		else:
			suff = '_shearcat'
		suff += '_forced'
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

	#whether or not to initially make the map at high resolution and degrade
	highres_first = False
	#NSIDE for the upgraded-resolution version of the bright object mask
	nside_mask = 16384

	#column names for flags identifying sources near bright objects
	bo_flags = [f'{cf_global.band}_mask_brightstar_ghost15',
				f'{cf_global.band}_mask_brightstar_halo',
				f'{cf_global.band}_mask_brightstar_blooming']
	
	#whether or not to include additional flags in the mask
	incl_main = True
	incl_strict = True


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

	#whether to smooth the map prior to removing pixels


########################
#### combine_fields ####
########################

class combineFields(cf_global):

	name = 'combineFields'


#########################
#### pca_systematics ####
#########################

class pcaSystematics(cf_global):

	name = 'pcaSystematics'

	#plot eigenvalues of the principal components
	plot_eigen = True
	#fraction of total variance to keep with principal components
	var_thresh = 0.98


############################
#### theory_predictions ####
############################

class theoryPredictions(cf_global):

	name = 'theoryPredictions'

	#column in the catalogues containing the random MC draws from the redshift distribution
	z_mc_col = 'pz_mc_dnnz'
	#number of bins to use in the n(z) histograms
	nbins_nofz = 50

	#fiducial cosmology parameters
	cosmo_fiducial = {
		'Omega_c' : 0.27,
		'Omega_b' : 0.045,
		'h'       : 0.67,
		'sigma8'  : 0.83,
		'n_s'     : 0.96
	}

	#base name of the files to which theory power spectra will be saved
	theory_out = 'theory_cells.hdf5'

###############################
#### compute_power_spectra ####
###############################

class computePowerSpectra(cf_global):

	name = 'computePowerSpectra'

	#systematics maps to deproject (add 'all' to the list of wanting to deproject all
	#maps in the systematics directory)
	systs = [
		#'all'
		]
	#(optional) maximum number of systematics to deproject - uses all provided if set to None
	Nsyst_max = None

	#approximately logarithmically-spaced bandpowers used in Nicola+19
	use_N19_bps = False
	bpw_edges = [100, 200, 300, 400, 600, 800, 1000, 1400, 1800, 2200, 3000,
			 3800, 4600, 6200, 7800, 9400, 12600, 15800]
	#Number of bandpowers to use if not using edges from Nicola+19
	nbpws = 18
	#minimum ell (i.e. largest scale) to use
	ell_min = 1
	#whether to use linear or log spacing for the bandpowers
	log_spacing = False
	
	#output file for power spectrum information
	outfile = f'power_spectra_info_{cf_global.nside_hi}_nodegrade_smoothedmask.hdf5'

	#output files for the NmtWorkspace and NmtCovarianveWorkspace
	wsp_file = f'workspace_{cf_global.nside_hi}_nodegrade_smoothedmask.fits'
	covwsp_file = f'covworkspace_{cf_global.nside_hi}_nodegrade_smoothedmask.fits'
	#cache file for keeping track of which systematics have been deprojected previously
	deproj_file = f'deprojected_{cf_global.nside_hi}_nodegrade_smoothedmask.txt'

	#apply a multiplicative correction to delta_g due to star contamination
	correct_for_stars = True
	#fiducial estinate for the fraction of stars making it into the final sample (from Nicola+19)
	Fs_fiducial = 0.02 

	#create lightweight NmtFields (cannot calculate deproj. bias, but saves memory)
	lite = False



###############################
#### amend_power_spectra ####
###############################

class amendPowerSpectra(computePowerSpectra):

	name = 'amendPowerSpectra'


###############################
#### plot_power_spectra ####
###############################

class plotPowerSpectra(computePowerSpectra):

	name = 'plotPowerSpectra'

	#whether to multiply the C_ells by l(l+1)/(2*pi) for the figure
	normalise = False

	#show the C_ells pre-debiasing
	show_pre_debias = False
	#show the C_ells without deprojection
	show_no_deproj = True
	#show the theory predictions
	show_theory = True
	#make a figure showing all C_ells simultaneously
	make_combined = False


	
############################
#### make_txpipe_inputs ####
############################

class makeTXPipeInputs(cf_global):

	name = 'makeTXPipeInputs'


