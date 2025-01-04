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

		#set the map and catalogue names
		self._set_map_names()
		self._set_cat_names()
		
	def __getattr__(self, name):
		'''
		Enables retrieval of pipeline settings as object attributes
		rather than dictionary entries.
		'''
		value = self.config_dict[name]
		if isinstance(value, dict):
			value = DictAsMember(value)
		return value

	def _set_map_names(self):
		'''
		Appends the base names for each map file with the specified suffix
		and NSIDE used for the analysis.
		'''
		for key in self.maps:
			#if dustmaps, replace name with list of names for all bands
			if key == 'dustmaps':
				self.config_dict['maps'][key] = [self.maps[key] 
									 			+ f'_{b}_nside{self.nside_hi}{self.suffix}.hsp'
												for b in self.bands.all]
			else:
				self.config_dict['maps'][key] = self.maps[key] + f'_nside{self.nside_hi}{self.suffix}.hsp'

	def _set_cat_names(self):
		'''
		Appends the base names for each catalogue file with the specified suffix.
		'''
		for key in self.cats:
			self.config_dict['cats'][key] = self.cats[key] + self.suffix + '.hdf5'

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
		sample_masks: list[numpy.array]
			List of boolean arrays identifying which sources in the catalogue
			belong to each sample.
		'''
		import numpy as np

		sample_masks = []
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
			sample_masks.append(samp_mask)
		return sample_masks