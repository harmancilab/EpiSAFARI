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

void select_points_of_interest_per_RD_signal_profile(double* signal_profile, 
													int start_i, int end_i, 
													int max_dist_between_cons_pts, 
													vector<double>* x_vec, vector<double>* y_vec);

void bspline_encode_mapped_read_profile(char* signal_dir,
	char* chr_id,
	int l_frag,
	int n_spline_coeff,
	int bspline_order,
	int min_n_pts_2_encode,
	int max_dist_between_cons_pts,
	double max_max_err,
	double max_avg_err,
	int l_win);

#endif