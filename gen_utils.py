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

