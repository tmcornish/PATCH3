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
	

class cf_global(object):

	def __init__(self, config_file):
		'''
		Loads the contents of a YAML file and uses it to define
		properties representing settings for the pipeline.
		'''
		
		with open(config_file) as f:
			self.config_dict = yaml.safe_load(f)['global']
		
	def __getattr__(self, name):
		'''
		Enables retrieval of pipeline settings as object attributes
		rather than dictionary entries.
		'''
		value = self.config_dict[name]
		if isinstance(value, dict):
			value = DictAsMember(value)
		return value
