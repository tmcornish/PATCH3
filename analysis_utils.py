############################################################################################################
# Module containing generally useful functions for data analysis.
###########################################################################################################

import numpy as np

def fit_polynomial_analytic(x, y, cov, degree=1, return_cov=False, return_chi2nu=False):
	'''
	Uses linear algebra to analytically fit a polynomial of degree n to data.

	Parameters
	----------
	x: array-like
		Values of the independent variable.
	
	y: array-like
		Values of the dependent variable.
	
	cov: array-like
		Covariance matrix of the data.
	
	degree: int
		Degree of the polynomial to be fitted (default 1).
	
	return_cov: bool
		Whether to return the covariance matrix of the best-fit coefficients
		(default False).
	
	return_chi2nu: bool
		Whether to return the reduced chi^2 value.
		
	Returns
	-------
	BF: numpy.ndarray
		Best-fit coefficients.
	
	cov_bf: numpy.ndarray
		Covariance matrix of the best-fit coefficients (only returned if 
		return_cov is True).
	
	chi2_red: float
		Reduced chi^2 value (only returned if return_chi2nu is True).
	'''

	#construct an array where each column is x^0, x^1, ... x^n
	n = degree + 1
	A = np.vander(x, n, increasing=True)
	#invert the data covariance matrix
	icov = np.linalg.inv(cov)

	#set up the matrices for solving the linear system of equations
	M1 = A.T @ icov @ A
	M2 = A.T @ icov @ y

	#matrix containing the best-fit coefficients in ascending order of degree in x
	BF = np.linalg.inv(M1) @ M2
	to_return = [BF]

	if return_cov:
		#covariance matrix for the best-fit coefficients
		cov_bf = np.linalg.inv(M1)
		to_return.append(cov_bf)
	
	if return_chi2nu:
		#compute reduced chi^2
		M3 = y - A @ BF
		chi2 = M3.T @ icov @ M3
		chi2_red = chi2 / (len(y) - n)
		to_return.append(chi2_red)
	
	return (*to_return,)
