############################################################################################################
# Module containing convenience functions for calculating angular power spectra.
############################################################################################################

import pymaster as nmt
import numpy as np


def get_bin_pairings(nbins, auto_only=False):
	'''
	Returns pairs of IDs for each tomographic bin being analysed. Also
	returns comma-separated string versions of the ID pairs.

	Parameters
	----------
	nbins: int
		Number of bins being considered.
	
	auto_only: bool
		If True, will only return each bin paired with itself.
	
	Returns
	-------
	pairings: list[tuple]
		List of possible bin pairings.
	
	pairings_s: list[str]
		List of comma-separated string versions of the ID pairs.
	'''
	import itertools

	l = list(range(nbins))
	if auto_only:
		pairings = [(i,i) for i in l]
	else:
		pairings = [i for i in itertools.product(l,l) if tuple(reversed(i)) >= i]
	pairings_s = [f'{p[0]},{p[1]}' for p in pairings]
	return pairings, pairings_s


def get_data_from_sacc(s, auto_only=False):
	'''
	Retrieves the relevant covariance information, which depends on the
	user settings in config.py.

	Parameters
	----------
	s: sacc.sacc.Sacc
		Sacc object containing the angular power spectrum information.
	
	auto_only: bool
		If True, retrieves information for the autocorrelations only.

	Returns
	-------
	ells: np.ndarray
		Array of multipoles at which the power spectra are defined (i.e. effective
		multipoles of each bandpower).
	
	cells: np.ndarray
		Array of arrays containing the desired power spectra.
	
	cov: np.ndarray
		Matrix of the relevant convariances.
	'''

	#get the tracer combinations
	combos = s.get_tracer_combinations()
	#get the number of tomographic bins
	nbins = len(np.unique(combos))
	#get the full covariance matrix from the Sacc
	cov_full = s.covariance.covmat

	if auto_only:
		#cycle through the redshift bins
		ells, cells, inds = [], [], []
		for i in range(nbins):
			ells_now, cells_now, inds_now = s.get_ell_cl('cl_00', f'bin_{i}', f'bin_{i}', return_ind=True)
			ells.append(ells_now)
			cells.append(cells_now)
			inds.extend(list(inds_now))
		#use the returned indices to retrieve the relevant covariance information
		cov = cov_full[inds][:,inds]
	
	else:
		#get data for all power spectra
		ell_cl = [s.get_ell_cl('cl_00', i, j) for i, j in combos]
		cov = cov_full
		#list the ells and cells for each pairing
		ells = [ell_cl[i][0] for i in range(len(ell_cl))]
		cells = [ell_cl[i][1] for i in range(len(ell_cl))]
	
	return ells, cells, cov
	

def select_from_sacc(s, tracer_combos, data_type):
	'''
	Given a Sacc object and a set of tracer combinations, will return a new Sacc object
	containing only the information for those tracer combinations.

	Parameters
	----------
	s: sacc.sacc.Sacc or str
		The Sacc object containing the information for many tracers. If a string, must be 
		the path to a Sacc file.
	
	tracer_combos: list[tuple]
		List of tuples, with each tuple containing a pair of tracer names.
	
	data_type: str
		Data type for which the information is to be extracted. E.g. 'cl_00' or
		'galaxy_density_cl'. Use print(sacc.standard_types) to see list of possible
		values.
	
	Returns
	-------
	s_new: sacc.sacc.Sacc
		Sacc object containing only the desired information.
	'''
	import sacc

	#check if input is a string
	if type(s) == str:
		s = sacc.Sacc.load_fits(s)

	#get the unique tracer names
	tc_unique = np.unique(tracer_combos)
	#set up a new Sacc object and add tracers
	s_new = sacc.Sacc()
	for tc in tc_unique:
		s_new.add_tracer_object(s.get_tracer(tc))
	#now add ell and C_ell info for each desired combination
	inds_all = []
	for tc in tracer_combos:
		ells, cells, inds = s.get_ell_cl(data_type, *tc, return_ind=True)
		#get the window functions
		wins = s.get_bandpower_windows(inds)
		#add the ell_cl info
		s_new.add_ell_cl(data_type, *tc, ells, cells, window=wins)
		#add the indices to the list
		inds_all.extend(list(inds))
	#add the covariance info
	if s.covariance is not None:
		s_new.covariance = sacc.covariance.FullCovariance(s.covariance.covmat[inds_all][:,inds_all])
	else:
		s_new.covariance = None

	return s_new
	


def compute_covariance(w, cw, f_i1, f_i2, f_j1=None, f_j2=None, f_sky=None, return_cl_coupled=False, return_cl_guess=False):
	'''
	Computes the Gaussian covariance for a pair of angular power spectra,
	given the fields used to compute said spectra. 
	NOTE: this currently assumes the same mask is used for all fields
	involved.

	Parameters
	----------
	w: pymaster.NmtWorkspace
		NaMaster workspace defined with the mask applied to all invovled fields.
	
	cw: pymaster.NmtCovarianceWorkspace
		Covariance workspace defined with the mask applied to all invovled 
		fields.
	
	f_i1, f_i2: pymaster.NmtField, pymaster.NmtField
		Fields contributing to the first power spectrum.
	
	f_j1, f_j2: pymaster.NmtField, pymaster.NmtField
		Fields contributing to the second power spectrum. If None, will be set
		to f_i1 and f_i2.
	
	f_sky: float
		Estimate of the observed sky fraction (assumed the same for all fields).
		If None, will calculate using the mask.
	
	return_cl_coupled: bool
		Whether to return the coupled C_ells computed as part of the covariance
		calculation.
	
	return_cl_guess: bool
		Whether to return the best-guess coupled C_ells computed as part of the 
		covariance calculation.
	
	Returns
	-------
	covar: numpy.ndarray
		Covariance matrix for the pair of angular power spectra.
	
	err_cell: numpy.ndarray
		1-sigma uncertainty for the cross-correlation in each bandpower.
	
	cl_coupled_X1Y2: list[numpy.ndarray]
		(Optional) Coupled C_ells for each field pair.
	
	cl_guess_X1Y2: list[numpy.ndarray]
		(Optional) Best-guess coupled C_ells for each field pair.
	'''
	#see if additional fields have been provided for the second power spectrum
	if f_j1 is None:
		f_j1 = f_i1
	if f_j2 is None:
		f_j2 = f_i2
	
	#if no sky fraction estimate is provided, compute from mask
	if f_sky is None:
		f_sky = np.mean(f_i1.get_mask() ** 2.)
	
	#compute coupled c_ells for each possible combination of i and j
	cl_coupled_i1j1 = nmt.compute_coupled_cell(f_i1, f_j1)
	cl_coupled_i1j2 = nmt.compute_coupled_cell(f_i1, f_j2)
	cl_coupled_i2j1 = nmt.compute_coupled_cell(f_i2, f_j1)
	cl_coupled_i2j2 = nmt.compute_coupled_cell(f_i2, f_j2)
	#use these along with the mask to get a guess of the true C_ell
	cl_guess_i1j1 = cl_coupled_i1j1 / f_sky
	cl_guess_i1j2 = cl_coupled_i1j2 / f_sky
	cl_guess_i2j1 = cl_coupled_i2j1 / f_sky
	cl_guess_i2j2 = cl_coupled_i2j2 / f_sky


	covar = nmt.gaussian_covariance(cw, 
									0, 0, 0, 0,			#spin of each field
									[cl_guess_i1j1[0]],	
									[cl_guess_i1j2[0]],
									[cl_guess_i2j1[0]],
									[cl_guess_i2j2[0]],
									w)
	#errorbars for each bandpower
	err_cell = np.diag(covar) ** 0.5

	to_return = [covar, err_cell]
	if return_cl_coupled:
		to_return.append([cl_coupled_i1j1, cl_coupled_i1j2, cl_coupled_i2j1, cl_coupled_i2j2])
	if return_cl_guess:
		to_return.append([cl_guess_i1j1, cl_guess_i1j2, cl_guess_i2j1, cl_guess_i2j2])


	return (*to_return,)