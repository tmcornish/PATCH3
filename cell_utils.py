############################################################################################################
# Module containing convenience functions for calculating angular power spectra.
############################################################################################################

import pymaster as nmt
import numpy as np

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