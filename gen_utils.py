############################################################################################################
# Module containing general-purpose functions.
###########################################################################################################

import numpy as np


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

