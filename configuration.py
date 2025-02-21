#####################################################################################################
# Configuration class for the various stages of the GC-HSC3 pipeline.
#####################################################################################################

import yaml

class DictAsMember(dict):
	'''
	Class for converting dictionary entries into object
	attributes. Allows for nested dictionaries.
	'''
	def __getattr__(self, name):
		value = self[name]
		if isinstance(value, dict):
			value = DictAsMember(value)
		return value
	

class PipelineConfig():
	''' 
	A class for identifying and defining the settings for the pipeline,
	both globally and at the specified stage.
	'''

	def __init__(self, config_file, stage=None):
		'''
		Loads the contents of a YAML file and uses it to define
		properties representing settings for the pipeline.
		'''
		
		with open(config_file) as f:
			config_all = yaml.safe_load(f)
			self.config_dict = config_all['global']
			if stage is not None:
				self.config_dict = {**self.config_dict, **config_all[stage]}
		
		#define a new property containing all bands
		self.config_dict['bands']['all'] = [self.bands.primary] + self.bands.secondary

		#get the number of samples and set it as a property
		self.config_dict['nsamples'] = len(self.samples)

		#store the path of the config file
		self.config_dict['config_file'] = config_file
		#set the map and catalogue names
		self._set_output_names()
		#identify the machine or cluster on which this is being run
		self._set_platform()
		#identify the various flags to be applied when removing sources
		self._set_flags()		


	def __getattr__(self, name):
		'''
		Enables retrieval of pipeline settings as object attributes
		rather than dictionary entries.
		'''
		value = self.config_dict[name]
		if isinstance(value, dict):
			value = DictAsMember(value)
		return value


	def _set_output_names(self):
		'''
		Appends the base names for each map file with the specified suffix
		and (where appropriate) the NSIDE used for the analysis.
		'''
		#hdf5 catalogues
		for key in self.cats:
			self.config_dict['cats'][key] = self.cats[key] + '.hdf5'
		
		#healsparse maps
		for key in self.maps:
			#if dustmaps, replace name with list of names for all bands
			if key == 'dustmaps':
				self.config_dict['maps'][key] = [self.maps[key] 
									 			+ f'_{b}_nside{self.nside_hi}{self.suffix}.hsp'
												for b in self.bands.all]
			else:
				self.config_dict['maps'][key] = self.maps[key] + f'_nside{self.nside_hi}{self.suffix}.hsp'
		
		#data files
		for key in self.data_files:
			self.config_dict['data_files'][key] = self.paths.data + self.data_files[key] + '.fits'

		#n(z) hdf5 files
		for key in self.nofz_files:
			self.config_dict['nofz_files'][key] = self.nofz_files[key] + self.suffix + '.hdf5'
		
		#power spectra hdf5 files
		for key in self.cell_files:
			self.config_dict['cell_files'][key] = self.cell_files[key] + f'_nside{self.nside_hi}{self.suffix}.hdf5'
		
		#cache files
		for key in self.cache_files.workspaces:
			self.config_dict['cache_files']['workspaces'][key] = self.cache_files.workspaces[key] + \
																f'_nside{self.nside_hi}{self.suffix}.fits'
		for key in self.cache_files.deproj:
			self.config_dict['cache_files']['deproj'][key] = self.cache_files.deproj[key] + \
																f'_nside{self.nside_hi}{self.suffix}.txt'
		for key in self.cache_files.hods:
			self.config_dict['cache_files']['hods'][key] = self.cache_files.hods[key] + \
																f'_nside{self.nside_hi}{self.suffix}.txt'
		
		#output SACC files
		for key in self.sacc_files:
			self.config_dict['sacc_files'][key] = self.sacc_files[key] + \
													f'_nside{self.nside_hi}{self.suffix}.fits'
		
		#auxiliary files (assumed to be stored in the same directory as the config file)
		conf_dir = '/'.join(self.config_file.split('/')[:-1]) + '/'
		for key in self.auxfiles:
			self.config_dict['auxfiles'][key] = conf_dir + self.auxfiles[key]


	def _set_platform(self):
		'''
		Identifies the machine or cluster on which the pipeline is being run, and sets it as
		a property of the class.
		'''
		import platform as pf
		#get the name of the node on which this is being run
		node = pf.node()
		if node.startswith('comp'):
			node = 'glamdring'
		elif node.startswith('nid'):
			node = 'nersc'
		else:
			node = 'local'
		
		#set the name of the node as a property of the class
		self.config_dict['platform'] = node
	

	def _set_flags(self):
		'''
		Retrieves the column names for all flags to be used when excluding sources from analysis,
		using names listed in the file specified as self.auxfiles.flags.
		'''
		self.config_dict['flags'] = {}
		with open(self.auxfiles.flags) as f:
			flags_dict = yaml.safe_load(f)

		#brightstar flags are assumed to be associated with the primary band unless they begin with the band name
		self.config_dict['flags']['brightstar'] = [f'{self.bands.primary}{fl}' if fl.startswith('_') else fl 
											 		for fl in flags_dict['brightstar']]
		
		#'main' flags are associated with the primary band
		self.config_dict['flags']['main'] = [f'{self.bands.primary}{fl}' for fl in flags_dict['main']]

		#'strict' flags are associated with all bands
		import itertools
		self.config_dict['flags']['strict'] = list(
			itertools.chain.from_iterable(
				[f'{b}{fl}' for fl in flags_dict['strict']]
				for b in self.bands.all
			)
		)


	def get_subfields(self):
		'''
		Identifies which subfields belong to the fields specified in the config file.
		'''
		subfields = []
		if 'hectomap' in self.fields:
			subfields.append('hectomap')
		if 'spring' in self.fields:
			subfields.extend([f'equator{i:02d}' for i in [21,22,23,0,1,2]])
		if 'autumn' in self.fields:
			subfields.extend([f'equator{i:02d}' for i in [8,9,10,11,12,13,14,15]])
		if 'cosmos' in self.fields:
			subfields.append('cosmos')

		return subfields


	def get_samples(self, cat):
		'''
		Returns masks to apply to the input catalogue in order to select sources
		belonging to each sample.

		Parameters
		----------
		cat: h5py._hl.files.File or h5py._hl.group.Group
			The input catalogue or Group. Must contain the relevant Datasets
			needed to define the cuts.
		
		Returns
		-------
		sample_masks: dict[numpy.array]
			Dictionary of boolean arrays identifying which sources in the catalogue
			belong to each sample.
		'''
		import numpy as np

		sample_masks = {}
		for s in self.samples:
			samp = self.samples[s]
			#replace any references to key columns with their actual names
			for key in self.key_cols:
				samp = samp.replace(key, f'cat["{self.key_cols[key]}"][:]')
			#split the expression at any semicolons and evaluate each term
			samp = samp.split(';')
			
			samp_mask = np.ones_like(cat[next(iter(cat.keys()))][:], dtype=bool)
			for ex in samp:
				samp_mask *= eval(ex)
			sample_masks[s] = samp_mask
		return sample_masks
	

	@staticmethod
	def get_field_boundaries(field):
		'''
		Given the name of an HSC field, will return the approximate corner coordinates of
		a rectangular boundary encompassing the field. These coordinates are listed as
		[RA_min, RA_max, Dec_min, Dec_max].

		Parameters
		----------
		field: str
			Name of the field whose boundaries are to be returned. Must be either 'hectomap',
			'spring', 'autumn', or 'aegis'.
		
		Returns
		-------
		bounds: list[float]
			Coordinates defining the boundary of the field, given as [RA_min, RA_max, Dec_min,
			Dec_max]. In the event that the field crosses RA=0, RA_min will lie westward of 
			this longitude, and RA_max will lie eastward (in this situation, RA_min > RA_max).
		'''
		if field == 'hectomap':
			bounds = [195., 255., 41.5, 45.]
		elif field == 'spring':
			bounds = [326.25, 41.25, -8., 8.]
		elif field == 'autumn':
			bounds = [125., 227.5, -4., 7.]
		elif field == 'aegis':
			bounds = [212., 216., 51.6, 53.6]
		else:
			raise ValueError('field must be either "hectomap", "spring", "autumn", or "aegis".')
		return bounds
	

	@staticmethod
	def combine_flags(dataset, flagcols, combine_type='or'):
		'''
		Combines flag columns with the specified operation ('and' or 'or').

		Parameters
		----------
		dataset: HDF5 Dataset or Group
			Dataset containing the relevant flags.

		flagcols: list
			List of all columns containing the flags to be combined.

		combine_type: str
			Operation to be used to combine the flags ('and' or 'or').

		Returns
		-------
		flags_comb: array-like
			The result of combining the flags.
		'''

		#see if valid combine_type has been chosen
		combine_type = combine_type.lower()
		if combine_type not in ['or', 'and']:
			raise TypeError('combine_type must be either "and" or "or" (capital letters allowed).')

		if combine_type == 'or':
			combine = lambda x,y : x + y
		else:
			combine = lambda x,y : x * y

		flags_comb = [dataset[fl][:] for fl in flagcols]
		from functools import reduce
		flags_comb = reduce(combine, flags_comb)

		return flags_comb