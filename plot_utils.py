############################################################################################################
# Module containing functions and variables for formatting plots.
###########################################################################################################


def check_for_latex():
	'''
	Checks if LaTeX is installed on machine; returns True if so and False otherwise.

	Returns
	-------
	installed: bool
		True if LaTeX is installed and False otherwise.
	'''
	import shutil
	
	if shutil.which('latex'):
		installed = True
	else:
		installed = False
	
	return installed


#dictionary containing custom formatting for plots
styledict = {
	'figure.figsize' : (8., 6.),
	'legend.fontsize' : 14,
	'legend.shadow' : False,
	'legend.framealpha' : 0.9,
	'xtick.labelsize' : 22,
	'ytick.labelsize' : 22,
	'axes.labelsize' : 24,
	'axes.linewidth' : 2.,
	'image.origin' : 'lower',
	'xtick.minor.visible' : True,
	'xtick.major.size' : 7,
	'xtick.minor.size' : 4,
	'xtick.major.width' : 2.,
	'xtick.minor.width' : 1.5,
	'xtick.top' : True,
	'xtick.major.top' : True,
	'xtick.minor.top' : True,
	'xtick.direction' : 'in',
	'ytick.minor.visible' : True,
	'ytick.major.size' : 7,
	'ytick.minor.size' : 4,
	'ytick.major.width' : 2.,
	'ytick.minor.width' : 1.5,
	'ytick.right' : True,
	'ytick.major.right' : True,
	'ytick.minor.right' : True,
	'font.size' : 22,
	'font.weight' : 'bold',
	'ytick.direction' : 'in',
	'text.usetex' : check_for_latex(),				#enables the use of LaTeX style fonts and symbols
	'mathtext.fontset' : 'stix',
	'font.family' : 'STIXGeneral',
	'axes.formatter.useoffset' : False,
}

#colours
red = '#eb0505'
dark_red = '#ab0000'
ruby = '#C90058'
crimson = '#AF0404'
coral = '#FF4848'
magenta = '#C3027D'
orange = '#ED5A01'
green = '#0A8600'
light_green = '#11C503'
teal = '#00A091'
cyan = '#00d0f0'
blue = '#0066ff'
light_blue = '#00C2F2'
dark_blue = '#004ab8'
purple = '#6D04C4'
lilac = '#EB89FF'
plum = '#862388'
pink = '#E40ACA'
baby_pink = '#FF89FD'
fuchsia = '#E102B5'
grey = '#969696'

#obtain the figure size in inches
x_size, y_size = 8., 8.

#formatting for any arrows to be added to the plot for representing upper/lower limits
ax_frac = 1/40.				#the fraction of the y axis that the total length of a vertical arrow should occupy
al = ax_frac * y_size		#the length of each arrow in inches (ew, but sadly metric isn't allowed)
scale = 1./al				#'scale' parameter used for defining the length of each arrow in a quiver
scale_units = 'inches'		#units of the scale parameter
aw = 0.0175 * al			#the width of each arrow shaft in inches
hw = 4.						#width of the arrowheads in units of shaft width
hl = 3.						#length of the arrowheads in units of shaft width
hal = 2.5					#length of the arrowheads at the point where they intersect the shaft 
							#(e.g. hal = hl gives a triangular head, hal < hl gives a more pointed head)
#put all of the above arrow settings into a dictionary
arrow_settings = {
	'scale' : scale,
	'scale_units' : scale_units,
	'width' : aw,
	'headwidth' : hw,
	'headlength' : hl,
	'headaxislength' : hal
}

#default markers to cycle through when plotting multiple datasets on a set of axes
cycle_mkr = ['s', '^', 'v', '*', (8,1,0), 'x', '<', '>', 'p']
#default colours to cycle through when plotting multiple datasets
cycle_clr = [plum, lilac, blue, teal, green, orange, red, dark_red]



###################
#### FUNCTIONS ####
###################

def scale_RGB_colour(rgb, scale_l=1., scale_s=1.):
	'''
	Takes any RGB code and scales its 'lightness' and saturation (according to the hls colour 
	model) by factors of scale_l and scale_s.

	Parameters
	----------
	rgb: array-like
		The RGB colour specifications.
	
	scale_l: float
		The factor by which to scale the lightness.
	
	scale_s: float
		The factor by which to scale the saturation.

	Returns
	-------
	rgb: array-like
		The RGB colour specifications of the scaled colour.
	'''
	import colorsys
	#convert the rgb to hls (hue, lightness, saturation)
	h, l, s = colorsys.rgb_to_hls(*rgb)
	#scale the lightness ad saturation and ensure the results are between 0 and 1
	l_new = max(0, min(1, l * scale_l))
	s_new = max(0, min(1, s * scale_s))
	#convert back to rgb and return the result
	return colorsys.hls_to_rgb(h, l_new, s_new)



def plot_correlation_matrix(S, **kwargs):
	'''
	Given a Sacc object, computes the correlation matrix from the covariance
	matrix and creates a figure displaying it.

	Parameters
	----------
	S: sacc.sacc.Sacc object or str
		The Sacc object containing the covariance matrix to be converted into
		a correlation matrix. If a string, must be the path to a Sacc file.
	
	Returns
	-------
	fig: matplotlib.figure.Figure
		Figure in which the correlation matrix is displayed.
	
	ax: matplotlib.axes._axes.Axes
		Axes on which the correlation matrix is displayed.
	'''

	import matplotlib as mpl
	import matplotlib.pyplot as plt
	import sacc
	import numpy as np

	plt.style.use(styledict)
	#check if input is a string
	if type(S) == str:
		S = sacc.Sacc.load_fits(S)
	#retrieve the covariance matrix
	cov = S.covariance.covmat
	#convert to correlation matrix
	Dinv = np.diag((1. / np.sqrt(np.diag(cov))))
	corr = Dinv @ cov @ Dinv

	#display the correlation matrix on a figure
	f, ax = plt.subplots()
	im = ax.imshow(corr, origin='upper', **kwargs)
	cbar = f.colorbar(im, ax=ax)
	cbar.set_label(r'$\rho$')

	#adjust the axis labels so that they display the power spectra pairings
	#(assumes the last character in each C_ell label describes the bin)
	pairs = [f'{i[-1]}{j[-1]}' for i,j in S.get_tracer_combinations()]
	npairs = len(pairs)
	#determine the width of each sub-grid in the correlation matrix
	wsg = corr.shape[0] / npairs
	#tick positions
	tick_pos = np.arange(0.5 * wsg, (npairs+0.5)*wsg, wsg)
	#set the positions and labels
	ax.set_xticks(tick_pos, labels=pairs)
	ax.set_yticks(tick_pos, labels=pairs)
	ax.minorticks_off()
	ax.tick_params(direction='out')

	#return the figure and axes
	return f, ax


def plot_map(mp, field, vals_unseen=None, unseen_thresh=None, title='', **kwargs):
	'''
	Given a map (HealPIX or HealSparse format), will return a figure displaying
	that map. 

	Parameters
	----------
	mp: numpy.ndarray or healsparse.HealSparseMap.HealSparseMap
		The map to be displayed.
	
	field: str
		name of the HSC field being mapped. Must be either "hectomap",
		"spring" or "autumn".
	
	vals_unseen: list or None
		List of values to be treated as equivalent to hp.UNSEEN for display purposes.
		If None, all pixels will be displayed with their exact value.
	
	unseen_thresh: float or None
		Value below which pixels will be counted as UNSEEN. If None, all pixels will 
		be displayed with their exact value.
	
	Returns
	-------
	fig: matplotlib.figure.Figure
		Figure in which the map is displayed.
	'''
	import healpy as hp
	import healsparse as hsp
	import matplotlib.pyplot as plt
	plt.style.use(styledict)

	#approximate field centers and corresponding suitable image dimensions
	if field == 'hectomap':
		ra_mean = 231.71
		dec_mean = 43.34
		xsize = 500
		ysize = 100
	elif field == 'spring':
		ra_mean = 5
		dec_mean = 0
		xsize = 1500
		ysize = 300
	elif field == 'autumn':
		ra_mean = 177.16
		dec_mean = 1.05
		xsize = 2500
		ysize = 300
	else:
		raise ValueError("Argument 'field' must be one of {'hectomap', 'spring', 'autumn'}.")

	#get the resolution of the map
	if type(mp) == hsp.HealSparseMap:
		nside = mp.nside_sparse
		mp = mp.generate_healpix_map(nest=False)
	else:
		nside = hp.npix2nside(len(mp))
	reso = hp.nside2resol(nside, arcmin=True)
	#make a copy of the input map to avoid editing in place
	mp = mp.copy()
	
	#create the figure
	fig = plt.figure(figsize=(18, 5))
	if vals_unseen is not None:
		mp = mp.copy()
		for val in vals_unseen:
			mp[mp == val] = hp.UNSEEN
	if unseen_thresh is not None:
		mp = mp.copy()
		mp[mp <= unseen_thresh] = hp.UNSEEN
	hp.gnomview(mp, rot=[ra_mean, dec_mean, 0], xsize=xsize, ysize=ysize, reso=reso, notext=True, fig=fig, title=title, **kwargs)