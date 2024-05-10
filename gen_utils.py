############################################################################################################
# Module containing general-purpose functions.
###########################################################################################################

import numpy as np

def round_sigfigs(num, sf):
	'''
	Rounds a number to a given number of significant figures

	Parameters
	----------
	num: float
		The number to be rounded
		
	sf: int
		The number of significant figures to which num will be rounded.

	Returns
	-------
	num_rounded: float
		The input number rounded to the desired number of significant figures.
	'''
	if num != 0.:
		i = -int(np.floor(np.log10(abs(num))) - (sf - 1))		#the number of decimal places to round to
		num_rounded = round(num, i)
	else:
		num_rounded = 0.
	return num_rounded


def error_message(module, message):
	'''
	Prints a nicely formatted error message. For use within other modules as a means of identifying
	issues. 

	Parameters
	----------
	module: str
		The name of the module being debugged.

	message: str
		The body of the error message to print.
	'''
	from output_utils import colour_string
	err_str = [
		colour_string(f'{module}\n', 'cyan'),
		colour_string('Error: '),
		colour_string(message, 'white')]
	print(''.join(err_str))


def get_ndim(x):
	'''
	Takes a variable and determines its the number of dimensions, e.g. 0 for single number (scalar), 
	1 for a list or 1d array, 2 for a 2d array, etc. Note: only really effective with int, float, 
	str, list, tuple, and other array-like types.

	Parameters
	----------
	x: array-like
		The variable for which the rank is to be determined.

	Returns
	-------
	ndim: int
		Number of dimensions in x.
	'''
	xtype = type(x)

	#if the variable is a single integer, float or string, set ndim = 0
	if xtype in [int, float, str, np.float32, np.float64, np.int32, np.int64]:
		ndim = 0
	#otherwise, attempt to convert to an array then find the number of dimensions
	else:
		x_array = np.array(x)
		ndim = len(x_array.shape)
		#if the length is still 0, return None instead and print an error message
		if ndim == 0:
			ndim = None
			error_message('general.get_ndim', f'could not find ndim for data type {xtype.__name__}')

	return ndim