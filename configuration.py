##############################################################################
# Configuration class for the various stages of the PATCH3 pipeline.
##############################################################################

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

        # Add run name to name output directory
        self.config_dict['paths']['out'] += self.run_name + '/'
        # Handle case where all fields are to be run
        if 'all' in self.config_dict['fields']:
            self.config_dict['fields'] = [
                'hectomap',
                'spring',
                'autumn'
            ]

        # Define a new property containing all bands
        self.config_dict['bands']['all'] = [self.bands.primary] + \
            self.bands.secondary

        # Identify the machine or cluster on which this is being run
        self._set_platform()

    def __getattr__(self, name):
        '''
        Enables retrieval of pipeline settings as object attributes
        rather than dictionary entries.
        '''
        value = self.config_dict[name]
        if isinstance(value, dict):
            value = DictAsMember(value)
        return value

    def _set_platform(self):
        '''
        Identifies the machine or cluster on which the pipeline is being run,
        and sets it as a property of the class.
        '''
        import platform as pf
        # Get the name of the node on which this is being run
        node = pf.node()
        if node.startswith('cx3') or node.startswith('login-'):
            node = 'imperial-cx3'
        else:
            node = 'local'

        # Set the name of the node as a property of the class
        self.config_dict['platform'] = node

    def get_subfields(self):
        '''
        Identifies which subfields belong to the fields specified in the
        config file.
        '''
        subfields = []
        if 'hectomap' in self.fields:
            subfields.append('hectomap')
        if 'spring' in self.fields:
            subfields.extend([f'equator{i:02d}'
                              for i in [21, 22, 23, 0, 1, 2]])
        if 'autumn' in self.fields:
            subfields.extend([f'equator{i:02d}'
                              for i in [8, 9, 10, 11, 12, 13, 14, 15]])
        if 'cosmos' in self.fields:
            subfields.append('cosmos')
        if 'aegis' in self.fields:
            subfields.append('aegis')

        return subfields

    @staticmethod
    def get_field_boundaries(field):
        '''
        Given the name of an HSC field, will return the approximate corner
        coordinates of a rectangular boundary encompassing the field. These
        coordinates are listed as [RA_min, RA_max, Dec_min, Dec_max].

        Parameters
        ----------
        field: str
            Name of the field whose boundaries are to be returned. Must be
            either 'hectomap', 'spring', 'autumn', or 'aegis'.

        Returns
        -------
        bounds: list[float]
            Coordinates defining the boundary of the field, given as [RA_min,
            RA_max, Dec_min, Dec_max]. In the event that the field crosses
            RA=0, RA_min will lie westward of this longitude, and RA_max will
            lie eastward (in this situation, RA_min > RA_max).
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
            raise ValueError('field must be either "hectomap", "spring", '
                             '"autumn", or "aegis".')
        return bounds
