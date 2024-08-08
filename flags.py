#####################################################################################################
# File containing column names for the various flags in the HSC data,
# and means of retrieving them for analysis.
#####################################################################################################

#only sources with these flags set to False must be included in analysis
flags_keep = [
	'isprimary'
]

#flags for masking regions contaminated by bright stars
flags_brightstar = [
	'_mask_brightstar_ghost15',
	'_mask_brightstar_halo',
	'_mask_brightstar_blooming',
	#'_mask_brightstar_dip'
]
flag_channelstop = 'y_mask_brightstar_channel_stop'

#flags associated with the primary photometric band
flags_main = [
	'_cmodel_flag_badcentroid',
	'_sdsscentroid_flag',
	'_pixelflags_edge',
	'_pixelflags_interpolatedcenter',
	'_pixelflags_saturatedcenter',
	'_pixelflags_crcenter',
	'_pixelflags_bad',
	'_pixelflags_suspectcenter',
	'_pixelflags_clipped',
	'_deblend_skipped'
]

#flags associated with all photometric bands (considered 'strict' cuts)
flags_strict = [
	'_sdsscentroid_flag',
	'_cmodel_flag',
	'_psfflux_flag',
	'_pixelflags_edge',
	'_pixelflags_interpolatedcenter',
	'_pixelflags_saturatedcenter',
	'_pixelflags_crcenter',
	'_pixelflags_bad'
]


def get_flags(b_primary, b_secondary=[], types=['main'], incl_channelstop=False):
	'''
	Returns column names of the requested flags.

	Parameters
	----------
	b_primary: str
		The primary photometric band.
	
	b_secondary : list
		List of other bands included in the analysis.
	
	types: list
		Which types of flags to be returned.

	incl_channelstop: bool
		If True and brightstar flags requested, will include the channel stop
		flag in the y-band.
	
	Returns
	-------
	flags_all: list
		List of all requested flag column names.
	'''

	#list for containing flag column names
	flags_all = []
	#cycle through requested flag types
	for ty in types:

		#check if type is valid
		if ty not in ['main', 'strict', 'brightstar']:
			raise TypeError('Argument "types" must contain only combinations of "main", "strict", or "brightstar".')
		
		if ty == 'main':
			flags_all += [f'{b_primary}{fl}' for fl in flags_main]
		
		elif ty == 'strict':
			import itertools
			flags_all += list(
				itertools.chain.from_iterable(
					[f'{bs}{fl}' for fl in flags_strict]
					for bs in ([b_primary] + b_secondary)
				)
			)
		
		else:
			flags_all += [f'{b_primary}{fl}' for fl in flags_brightstar]
			if incl_channelstop:
				flags_all += [flag_channelstop]
	
	return list(set(flags_all))



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