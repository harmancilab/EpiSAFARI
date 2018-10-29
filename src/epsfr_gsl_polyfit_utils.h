#ifndef _POLIFITGSL_H
#define _POLIFITGSL_H

#include <vector>
using namespace std;

#ifdef __unix__
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_eigen.h>
	#include <gsl/gsl_sort.h>
	//#include <gsl/gsl_wavelet.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_errno.h>
	#include <gsl/gsl_fft_complex.h>
	#include <gsl/gsl_multifit.h>

	#define REAL(z,i) ((z)[2*(i)])
	#define IMAG(z,i) ((z)[2*(i)+1])
#endif

bool bsplinefit(int n_data_pts, 
				int n_spline_coeffs, 
				int spline_order,
				double *dx, double *dy, 
				double* reconst_y, 
				double *store);

#endif