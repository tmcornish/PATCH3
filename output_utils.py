############################################################################################################
# Module containing functions for formatting outputs from scripts.
###########################################################################################################


def colour_string(s, c='red'):
	'''
	Reformats a string so that it can be printed to the Terminal in colour (against a black background).

	Parameters
	----------
	s: str
		The string to be printed in colour.
	
	c: str
		The desired colour (must be one of the seven available choices; see below).

	Returns
	-------
	s_new: str
		The input string, reformatted to be displayed with the desired colour in the Terminal.
	'''

	#list of possible colour symbols
	colours = ['red', 'green', 'orange', 'blue', 'purple', 'cyan', 'white']
	#corresponding codes
	codes = ['0;31;40', '0;32;40', '0;33;40', '0;34;40', '0;35;40', '0;36;40', '0;37;40']
	#identify the code corresponding to the colour selected in the argument
	try:
		code_sel = codes[colours.index(c)]
		#use the relevant code to ensure the string is printed to the terminal in the chosen colour
		s_new = '\x1b[%sm%s\x1b[0m'%(code_sel, s)
		#return the new string
		return s_new
	except ValueError:
		#if the use did not select an available colour, print an explanation and return the original string
		print('colour_string: Selected colour not available. Available colours are:\nred\ngreen\norange\nblue\npurple\nolive\nwhite\n')
		return s


def string_important(s):
	'''
	Prints the provided string, using large numbers of '#'s to make it easy to spot when running the code.

	Parameters
	----------
	s: str
		String to print.

	Returns
	-------
	Input string reformatted to be easily visible in the Terminal.
	'''
	N_str = len(s)
	N_pad = 8
	N_total = N_str + 2 * (N_pad + 1) 
	pad_newline = '#' * N_total
	pad_textline = '#' * N_pad
	textline = ' '.join([pad_textline, s, pad_textline])	#line containing the important text with padding
	return '\n'.join([pad_newline, textline, pad_newline])


def array_to_fits(data, filename, CTYPE1='RA', CTYPE2='DEC', CRPIX=[1,1], CRVAL=[0,0], CDELT=[1,1]):
	'''
	Takes an array and details describing a coordinate system and creates a FITS file.

	Parameters
	----------
	data: array
		The array containing the data.
	
	filename: str
		Output file name.
	
	CTYPE1: str
		Name of the variable along axis 1.
	
	CTYPE2: str
		Name of the variable along axis 2.
	
	CRPIX: array-like
		Reference pixels in [X,Y].
	
	CRVAL: array-like
		Reference pixel values in [X,Y].
	
	CDELT: array-like
		Pixel scales in [X,Y].
	'''
	#create a PrimaryHDU from the chi^2 grid
	hdu = fits.PrimaryHDU(data)
	#update the relevant parameters in the header
	hdu.header['CTYPE1'] = CTYPE1
	hdu.header['CTYPE2'] = CTYPE2
	hdu.header['CRPIX1'] = CRPIX[0]
	hdu.header['CRPIX2'] = CRPIX[1]
	hdu.header['CRVAL1'] = CRVAL[0]
	hdu.header['CRVAL2'] = CRVAL[1]
	hdu.header['CDELT1'] = CDELT[0]
	hdu.header['CDELT2'] = CDELT[1]
	#create an HDUList
	hdul = fits.HDUList([hdu])
	#write this to the file
	hdul.writeto(filename, overwrite=True)


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
	err_str = [
		colour_string(f'{module}\n', 'cyan'),
		colour_string('Error: '),
		colour_string(message, 'white')]
	print(''.join(err_str))
