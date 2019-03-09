#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "epsfr_utils.h"
#include "epsfr_rng.h"
#include "epsfr_seed_manager.h"
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include "epsfr_gsl_polyfit_utils.h"

#include <vector>
#include <algorithm>
using namespace std;

bool __DUMP_POLYFIT_MESSAGES__ = false;

bool sort_derivs(double* deriv_entry1, double* deriv_entry2)
{
	return(deriv_entry1[1] > deriv_entry2[1]);

	//if (fabs(deriv_entry1[1] - deriv_entry2[1]) > 0.01)
	//{
	//	return(deriv_entry1[1] > deriv_entry2[1]);
	//}
	//else
	//{
	//	return(deriv_entry1[2] > deriv_entry2[2]);
	//}	
}

vector<double>* brkpt_select_uniform(double* dx, double* dy, int n_data_pts, int max_n_internal_brkpts)
{
	vector<double>* internal_brkpts = new vector<double>();
	double delta = (dx[n_data_pts - 1] + 0.01 - dx[0]) / (max_n_internal_brkpts + 2 - 1);
	double cur_brkpt = dx[0] + delta;
	for (int i = 0; i < max_n_internal_brkpts; i++)
	{
		internal_brkpts->push_back(cur_brkpt);
		cur_brkpt += delta;
	} // i loop.

	return(internal_brkpts);
}

// Randomly place the breakpoints.
vector<double>* brkpt_select_random(t_rng* rng, double* dx, double* dy, int n_data_pts, int max_n_internal_brkpts)
{
	vector<double>* internal_brkpts = new vector<double>();
	for (int i_brk = 0; i_brk < max_n_internal_brkpts; i_brk++)
	{
		internal_brkpts->push_back(dx[0] + rng->random_double_ran3() * (dx[n_data_pts-1] - dx[0] - 0.1));
	} // i_brk loop.

	sort(internal_brkpts->begin(), internal_brkpts->end());

	return(internal_brkpts);
}

vector<double>* brkpt_select_per_hill_derivative(double* dx, double* dy, int n_data_pts, int max_n_internal_brkpts)
{
	vector<double*>* derivs = new vector<double*>();

	// Compute the derivative of the signal.
	int i = 1;
	double cur_hill_deriv_sign = (dy[i] - dy[i - 1]);
	int cur_hill_start_i = 0;
	while ( i < n_data_pts)
	{
		double cur_pos_deriv_sign = (dy[i] - dy[i - 1]);

		// If the current sign is 0, reset it on the region.
		if (cur_hill_deriv_sign == 0)
		{
			cur_hill_deriv_sign = (dy[i] - dy[i - 1]);
		}

		if (cur_hill_deriv_sign * cur_pos_deriv_sign < 0)
		{
			double* cur_deriv = new double[3];
			cur_deriv[0] = (dx[i-1] + dx[cur_hill_start_i]) / 2; // Set the midpoint of the hill.
			cur_deriv[1] = fabs((dy[i-1] - dy[cur_hill_start_i]) / (dx[i-1] - dx[cur_hill_start_i]));
			cur_deriv[2] = (dy[i] + dy[i - 1]) / 2;

			if (__DUMP_POLYFIT_MESSAGES__)
			{
				fprintf(stderr, "Hill: %.1f-%.1f; Sign: %1.f\n", dx[cur_hill_start_i], dx[i - 1], cur_hill_deriv_sign);
			}

			derivs->push_back(cur_deriv);

			// Update the start and sign.
			cur_hill_start_i = i - 1;
			cur_hill_deriv_sign = (dy[i] - dy[i - 1]);
		}

		i++;
	} // i loop.

	// Add the last hill.
	double* cur_deriv = new double[3];
	cur_deriv[0] = (dx[i - 1] + dx[cur_hill_start_i]) / 2; // Set the midpoint of the hill.
	cur_deriv[1] = fabs((dy[i - 1] - dy[cur_hill_start_i]) / (dx[i - 1] - dx[cur_hill_start_i]));
	cur_deriv[2] = (dy[i] + dy[i - 1]) / 2;

	if (__DUMP_POLYFIT_MESSAGES__)
	{
		fprintf(stderr, "Hill: %.1f-%.1f\n", dx[cur_hill_start_i], dx[i - 1]);
	}

	derivs->push_back(cur_deriv);

	// Sort derivatives.
	sort(derivs->begin(), derivs->end(), sort_derivs);

	// Set the breakpoints as the points in the middle of regions with top derivatives.
	vector<double>* brkpts = new vector<double>();
	for (int i = 0; i < (int)derivs->size(); i++)
	{
		brkpts->push_back(derivs->at(i)[0]);

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(stderr, "Top brkpt @ %d: %.3f\n", i, derivs->at(i)[0]);
		}

		// Break, if we have n breakpoints.
		if ((int)brkpts->size() == max_n_internal_brkpts)
		{
			break;
		}
	} // i loop.

	// Add the remaining breakpoints uniformly.
	double delta = (dx[n_data_pts - 1] + 0.01 - dx[0]) / (max_n_internal_brkpts - (int)brkpts->size() + 2 - 1);
	double cur_brkpt = dx[0] + delta;

	for (int i = (int)brkpts->size(); i < max_n_internal_brkpts; i++)
	{
		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(stderr, "Uni brkpt @ %d: %.3f\n", i, cur_brkpt);
		}

		brkpts->push_back(cur_brkpt);
		cur_brkpt += delta;
	} // i loop.

	// Free memory.
	for (int i = 0; i < (int)derivs->size(); i++)
	{
		delete[] derivs->at(i);
	} // i loop.
	delete derivs;

	// Set the brkpt's at those locations.
	sort(brkpts->begin(), brkpts->end());

	return(brkpts);
}

//vector<double>* brkpt_select_per_window_derivative(double* dx, double* dy, int n_data_pts, int max_n_internal_brkpts)
//{
//	vector<double>* brkpts = new vector<double>();
//
//	double l_win = (dx[n_data_pts - 1] + 0.01 - dx[0]) / (max_n_internal_brkpts + 2 - 1);
//	double cur_win_start = dx[0];
//	for (int i_brkpt = 0; i_brkpt < max_n_internal_brkpts; i_brkpt++)
//	{
//		// find the highest derivative in the current window.
//		vector<double*>* derivs = new vector<double*>();
//		for (int i = 0; i < n_data_pts; i++)
//		{
//			// Make sure we are in the window.
//			if (dx[i] > cur_win_start &&
//				dx[i] < cur_win_start + l_win)
//			{
//				double* cur_deriv = new double[3];
//				cur_deriv[0] = (dx[i] + dx[i - 1]) / 2; // Set the midpoint.
//				cur_deriv[1] = fabs((dy[i] - dy[i - 1]) / (dx[i] - dx[i - 1]));
//				cur_deriv[2] = (dy[i] + dy[i - 1]) / 2;
//
//				derivs->push_back(cur_deriv);
//			}
//		} // i loop.
//
//		sort(derivs->begin(), derivs->end(), sort_derivs);
//
//		brkpts->push_back(derivs->at(0)[0]);
//
//		if (__DUMP_POLYFIT_MESSAGES__)
//		{
//			fprintf(stderr, "Window: %.1f-%.1f: Top deriv @ %.1f\n", 
//					cur_win_start, 
//					cur_win_start + l_win, 
//					derivs->at(0)[0]);
//		}
//
//		for (int i = 0; i < derivs->size(); i++)
//		{
//			delete[] derivs->at(i);
//		} // i loop.
//
//		delete derivs;
//
//		cur_win_start += l_win;
//	} // i_brkpt loop.
//
//	sort(brkpts->begin(), brkpts->end());
//
//	if (brkpts->size() != max_n_internal_brkpts)
//	{
//		fprintf(stderr, "Generated %d breakpoints (%d)\n", brkpts->size(), max_n_internal_brkpts);
//		exit(0);
//	}
//
//	return(brkpts);
//}
//
//vector<double>* brkpt_select_per_window_extrema(double* dx, double* dy, int n_data_pts, int max_n_internal_brkpts)
//{
//	vector<double>* brkpts = new vector<double>();
//
//	double l_win = (dx[n_data_pts - 1] + 0.01 - dx[0]) / (max_n_internal_brkpts + 2 - 1);
//	double cur_win_start = dx[0];
//	for (int i_brkpt = 0; i_brkpt < max_n_internal_brkpts; i_brkpt++)
//	{
//		// find the highest derivative in the current window.
//		vector<double*>* maxima = new vector<double*>();
//		for (int i = 0; i < n_data_pts; i++)
//		{
//			if (dx[i] > cur_win_start &&
//				dx[i] < cur_win_start + l_win)
//			{
//				double* cur_val = new double[3];
//				cur_val[0] = dx[i]; // Set the midpoint.
//				cur_val[1] = dy[i];
//				cur_val[2] = dy[i];
//
//				maxima->push_back(cur_val);
//			}
//		} // i loop.
//
//		sort(maxima->begin(), maxima->end(), sort_derivs);
//
//		brkpts->push_back(maxima->at(0)[0]);
//
//		if (__DUMP_POLYFIT_MESSAGES__)
//		{
//			fprintf(stderr, "Window: %.1f-%.1f: Top deriv @ %.1f\n",
//				cur_win_start,
//				cur_win_start + l_win,
//				maxima->at(0)[0]);
//		}
//
//		for (int i = 0; i < maxima->size(); i++)
//		{
//			delete[] maxima->at(i);
//		} // i loop.
//
//		delete maxima;
//
//		cur_win_start += l_win;
//	} // i_brkpt loop.
//
//	sort(brkpts->begin(), brkpts->end());
//
//	if (brkpts->size() != max_n_internal_brkpts)
//	{
//		fprintf(stderr, "Generated %d breakpoints (%d)\n", brkpts->size(), max_n_internal_brkpts);
//		exit(0);
//	}
//
//	return(brkpts);
//}

vector<double>* brkpt_select_per_vicinity_derivative(double* dx, double* dy, int n_data_pts, int max_n_internal_brkpts)
{
	vector<double*>* derivs = new vector<double*>();

	//int l_vic = (dx[n_data_pts-1] - dx[0]) / (max_n_internal_brkpts + 1);
	int l_vic = 100;

	// Compute the derivative of the signal.
	for (int i = 1; i < n_data_pts; i++)
	{
		int j_post = i;
		while (j_post < (n_data_pts-1) &&
			dx[j_post] < dx[i] + l_vic / 2)
		{
			j_post++;
		} // j loop.

		int j_pre = i;
		while (j_pre > 0 &&
			dx[j_pre] > dx[i] - l_vic / 2)
		{
			j_pre--;
		} // j loop.

		if(dx[j_post] >= dx[j_pre] + l_vic)
		{
			double* cur_deriv = new double[3];
			cur_deriv[0] = dx[i];
			//cur_deriv[1] = fabs((dy[j_post] - dy[j_pre]) / (dx[j_post] - dx[j_pre]));
			cur_deriv[1] = fabs((dy[j_post] - dy[j_pre]) / l_vic);
			cur_deriv[2] = (dy[j_post] + dy[j_pre]) / 2;

			derivs->push_back(cur_deriv);
		}

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(stderr, "@ %.1f: %.1f-%.1f: %.2f\n", dx[i], dx[j_pre], dx[j_post], (dy[j_post] - dy[j_pre]) / l_vic);
		}
	} // i loop.

	sort(derivs->begin(), derivs->end(), sort_derivs);

	// Set the breakpoints as the points in the middle of regions with top derivatives.
	vector<double>* brkpts = new vector<double>();
	for (int i = 0; i < (int)derivs->size(); i++)
	{
		brkpts->push_back(derivs->at(i)[0]);

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(stderr, "Top brkpt @ %d: %.3f\n", i, derivs->at(i)[0]);
		}

		// Break, if we have n breakpoints.
		if ((int)brkpts->size() == max_n_internal_brkpts)
		{
			break;
		}
	} // i loop.

	  // Set the brkpt's at those locations.
	sort(brkpts->begin(), brkpts->end());

	if ((int)brkpts->size() != max_n_internal_brkpts)
	{
		fprintf(stderr, "Too few breakpoints: %d, %d\n", (int)brkpts->size(), max_n_internal_brkpts);
		exit(0);
	}

	return(brkpts);
}


vector<double>* brkpt_select_per_derivative(double* dx, double* dy, int n_data_pts, int max_n_internal_brkpts)
{
	vector<double*>* derivs = new vector<double*>();

	// Compute the derivative of the signal.
	for (int i = 1; i < n_data_pts; i++)
	{
		double* cur_deriv = new double[3];
		cur_deriv[0] = (dx[i] + dx[i - 1]) / 2; // Set the midpoint.
		cur_deriv[1] = fabs((dy[i] - dy[i - 1]) / (dx[i] - dx[i - 1]));
		cur_deriv[2] = (dy[i] + dy[i - 1]) / 2;

		derivs->push_back(cur_deriv);
	} // i loop.

	sort(derivs->begin(), derivs->end(), sort_derivs);

	// Set the breakpoints as the points in the middle of regions with top derivatives.
	vector<double>* brkpts = new vector<double>();
	for (int i = 0; i < (int)derivs->size(); i++)
	{
		brkpts->push_back(derivs->at(i)[0]);

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(stderr, "Top brkpt @ %d: %.3f\n", i, derivs->at(i)[0]);
		}

		// Break, if we have n breakpoints.
		if ((int)brkpts->size() == max_n_internal_brkpts)
		{
			break;
		}
	} // i loop.

	// Set the brkpt's at those locations.
	sort(brkpts->begin(), brkpts->end());

	if ((int)brkpts->size() != max_n_internal_brkpts)
	{
		fprintf(stderr, "Too few breakpoints: %d, %d\n", (int)brkpts->size(), max_n_internal_brkpts);
		exit(0);
	}

	return(brkpts);
}

bool bsplinefit_nonuniform_per_breakpoint_type(int n_data_pts,
												int max_n_internal_breakpoints, // Maximum # of brkpts w/o 
												int breakpoint_type,
												int spline_order,
												double *dx, double *dy,
												double* reconst_y,
												double *store,
												t_rng* rng,
												bool dump_bsplines)
{
	//if (__DUMP_POLYFIT_MESSAGES__)
	//	fprintf(stderr, "Fitting %d points with B-spline with %d coeffs of order %d (%d points per fit).\n", n_data_pts, n_spline_coeffs, spline_order, (int)((dx[n_data_pts - 1] - dx[0]) / n_spline_coeffs));

	/* number of data points to fit */
	//#define N        200

	/* number of fit coefficients */
	//#define NCOEFFS  12

	/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
	//#define NBREAK   (NCOEFFS - 2)

	//const size_t n = n_data_pts;
	//const size_t ncoeffs = n_spline_coeffs;

	// Select breakpoints w/o the end points.
	vector<double>* internal_brkpts = NULL;
	if (breakpoint_type == DERIVATIVE_NU_BREAKPOINTS)
	{
		internal_brkpts = brkpt_select_per_derivative(dx, dy, n_data_pts, max_n_internal_breakpoints);
	}
	else if (breakpoint_type == VICINITY_DERIVATIVE_NU_BREAKPOINTS)
	{
		internal_brkpts = brkpt_select_per_vicinity_derivative(dx, dy, n_data_pts, max_n_internal_breakpoints);
	}
	else if (breakpoint_type == UNIFORM_BREAKPOINTS)
	{
		internal_brkpts = brkpt_select_uniform(dx, dy, n_data_pts, max_n_internal_breakpoints);
	}
	else if (breakpoint_type == RANDOM_NU_BREAKPOINTS)
	{
		internal_brkpts = brkpt_select_random(rng, dx, dy, n_data_pts, max_n_internal_breakpoints);
	}	
	//else if (breakpoint_type == PER_WINDOW_MAXIMA_BREAKPOINTS)
	//{
	//	internal_brkpts = brkpt_select_per_window_extrema(dx, dy, n_data_pts, max_n_internal_breakpoints);
	//}
	//else if (breakpoint_type == PER_WINDOW_DERIVATIVE_BREAKPOINTS)
	//{
	//	internal_brkpts = brkpt_select_per_window_derivative(dx, dy, n_data_pts, max_n_internal_breakpoints);
	//}
	else if (breakpoint_type == HILL_DERIVATIVE_NU_BREAKPOINTS)
	{
		internal_brkpts = brkpt_select_per_hill_derivative(dx, dy, n_data_pts, max_n_internal_breakpoints);
	}

	int n_breakpoints = (int)internal_brkpts->size() + 2; // Add the end points.

	// Set the number of coefficients: Note that the number of basis functions change based on the selected breakpoints.
	size_t ncoeffs = n_breakpoints + spline_order - 2;

	if (ncoeffs > (size_t)n_data_pts)
	{
		fprintf(stderr, "Too many coefficients to set non-uniformly.\n");
		exit(0);
	}

	// This is the number of knots.
	//const size_t nbreak = ncoeffs - (spline_order - 2);
	size_t nbreak = n_breakpoints;
	size_t i, j;
	gsl_bspline_workspace *bw;
	gsl_vector *B;
	//double dy;
	gsl_vector *c;
	//gsl_vector *w;
	gsl_vector *x, *y;
	gsl_matrix *X, *cov;
	gsl_multifit_linear_workspace *mw;
	double chisq;
	//double Rsq, dof, tss;

	//gsl_rng_env_setup();
	//r = gsl_rng_alloc(gsl_rng_default);

	/* allocate a cubic bspline workspace (k = 4) */
	bw = gsl_bspline_alloc(spline_order, nbreak);
	B = gsl_vector_alloc(ncoeffs);

	x = gsl_vector_alloc(n_data_pts);
	y = gsl_vector_alloc(n_data_pts);
	X = gsl_matrix_alloc(n_data_pts, ncoeffs);
	gsl_vector* breakpoints = gsl_vector_alloc(n_breakpoints);
	c = gsl_vector_alloc(ncoeffs);
	//w = gsl_vector_alloc(n_data_pts);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	mw = gsl_multifit_linear_alloc(n_data_pts, ncoeffs);

	//printf("#m=0,S=0\n");
	/* this is the data to be fitted */
	FILE* f_data = NULL;
	if (__DUMP_POLYFIT_MESSAGES__)
	{
		f_data = open_f("data.txt", "w");
	}

	for (i = 0; i < (unsigned int)n_data_pts; ++i)
	{
		gsl_vector_set(x, i, dx[i]);
		gsl_vector_set(y, i, dy[i]);

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(f_data, "%.2f\t%.2f\n", dx[i], dy[i]);
			fprintf(stderr, "%.2f\t%.2f\n", dx[i], dy[i]);
		}
	} // i loop.

	  // Close the data file if is it open.
	if (f_data != NULL)
	{
		fclose(f_data);
	}

	// Set the breakpoints vector, make sure to add the end points.
	// Add the first end point.
	gsl_vector_set(breakpoints, 0, dx[0]);

	// Add the breakpoints selected from derivative.
	for (int i_brkpt = 0; i_brkpt < (int)internal_brkpts->size(); i_brkpt++)
	{
		gsl_vector_set(breakpoints, i_brkpt+1, internal_brkpts->at(i_brkpt));
	} // i loop.

	// Add the last end point.
	gsl_vector_set(breakpoints, (int)internal_brkpts->size()+1, dx[n_data_pts - 1] + 0.01);
	delete internal_brkpts;

	  /* Use uniform breakpoints on [0, 15] */
	  //gsl_bspline_knots_uniform(dx[0], dx[n_data_pts - 1], bw);

	if (__DUMP_POLYFIT_MESSAGES__)
	{
		fprintf(stderr, "Setting the knots.\n");
	}

	gsl_bspline_knots(breakpoints, bw);

	// https://lists.gnu.org/archive/html/help-gsl/2012-01/msg00011.html
	if (__DUMP_POLYFIT_MESSAGES__)
	{
		gsl_vector_fprintf(stdout, bw->knots, "%f ");
	}

	/* construct the fit matrix X */
	FILE* f_splines = NULL;
	if (__DUMP_POLYFIT_MESSAGES__ ||
		dump_bsplines)
	{
		f_splines = fopen("basis_splines.txt", "a");
	}

	for (i = 0; i < (unsigned int)n_data_pts; ++i)
	{
		double xi = gsl_vector_get(x, i);

		/* compute B_j(xi) for all j */
		gsl_bspline_eval(xi, B, bw);

		/* fill in row i of X */
		for (j = 0; j < ncoeffs; ++j)
		{
			double Bj = gsl_vector_get(B, j);
			gsl_matrix_set(X, i, j, Bj);

			if (__DUMP_POLYFIT_MESSAGES__)
			{
				fprintf(f_splines, "%lf\t", Bj);
			}
		}

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(f_splines, "\n");
		}
	} // i loop.

	if (__DUMP_POLYFIT_MESSAGES__)
	{
		fclose(f_splines);
	}

	/* do the fit */
	gsl_multifit_linear(X, y, c, cov, &chisq, mw);

	//dof = n_data_pts - ncoeffs;
	//tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
	//Rsq = 1.0 - chisq / tss;

	//fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", 
	//                 chisq / dof, Rsq);

	// Reset the covariance matrix of the coefficients.
	for (unsigned int i = 0; i < ncoeffs; i++)
	{
		for (unsigned int j = 0; j < ncoeffs; j++)
		{
			gsl_matrix_set(cov, i, j, 0);
		} // j loop.
	} // i loop.

	/* output the smoothed curve */
	{
		double xi, yi, yerr;

		// Compute over the whole range.
		for (int i = 0; i < n_data_pts; i++)
		{
			xi = dx[i];
			gsl_bspline_eval(xi, B, bw);
			gsl_multifit_linear_est(B, c, cov, &yi, &yerr);

			if (gsl_finite(yi))
			{
				reconst_y[i] = yi;
			}
			else
			{
				reconst_y[i] = 0;
			}
			//printf("%f %f\n", xi, yi);
		}
	}

	//gsl_rng_free(r);
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	//gsl_vector_free(w);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);

	return true;
}

// This is an example of non-uniform bspline fitting.
bool bsplinefit_nonuniform_example(int n_data_pts,
	//int n_spline_coeffs,
	int n_breakpoints,
	int spline_order,
	double *dx, double *dy,
	double* reconst_y,
	double *store)
{
	//if (__DUMP_POLYFIT_MESSAGES__)
	//	fprintf(stderr, "Fitting %d points with B-spline with %d coeffs of order %d (%d points per fit).\n", n_data_pts, n_spline_coeffs, spline_order, (int)((dx[n_data_pts - 1] - dx[0]) / n_spline_coeffs));

	/* number of data points to fit */
	//#define N        200

	/* number of fit coefficients */
	//#define NCOEFFS  12

	/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
	//#define NBREAK   (NCOEFFS - 2)

	//const size_t n = n_data_pts;
	//const size_t ncoeffs = n_spline_coeffs;
	size_t ncoeffs = n_breakpoints + spline_order - 2;

	if (ncoeffs > (size_t)n_data_pts)
	{
		fprintf(stderr, "Too many coefficients to set non-uniformly.\n");
		exit(0);
	}

	// Use the data points x axis as the breakpoints.

	// This is the number of knots.
	//const size_t nbreak = ncoeffs - (spline_order - 2);
	size_t nbreak = n_breakpoints;
	size_t i, j;
	gsl_bspline_workspace *bw;
	gsl_vector *B;
	//double dy;
	gsl_vector *c;
	//gsl_vector *w;
	gsl_vector *x, *y;
	gsl_matrix *X, *cov;
	gsl_multifit_linear_workspace *mw;
	double chisq;
	//double Rsq, dof, tss;

	//gsl_rng_env_setup();
	//r = gsl_rng_alloc(gsl_rng_default);

	/* allocate a cubic bspline workspace (k = 4) */
	bw = gsl_bspline_alloc(spline_order, nbreak);
	B = gsl_vector_alloc(ncoeffs);

	x = gsl_vector_alloc(n_data_pts);
	y = gsl_vector_alloc(n_data_pts);
	X = gsl_matrix_alloc(n_data_pts, ncoeffs);
	gsl_vector* breakpoints = gsl_vector_alloc(n_breakpoints);
	c = gsl_vector_alloc(ncoeffs);
	//w = gsl_vector_alloc(n_data_pts);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	mw = gsl_multifit_linear_alloc(n_data_pts, ncoeffs);

	//printf("#m=0,S=0\n");
	/* this is the data to be fitted */
	FILE* f_data = NULL;
	if (__DUMP_POLYFIT_MESSAGES__)
	{
		f_data = open_f("data.txt", "w");
	}

	for (i = 0; i < (unsigned int)n_data_pts; ++i)
	{
		gsl_vector_set(x, i, dx[i]);
		gsl_vector_set(y, i, dy[i]);

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(f_data, "%.2f\t%.2f\n", dx[i], dy[i]);
		}
	} // i loop.

	// Close the data file if is it open.
	if (f_data != NULL)
	{
		fclose(f_data);
	}

	// This is the most vital portion, the selection of the breakpoints.
	// Set uniform breakpoints.
	// Add a small delta to make sure it the division is stable.
	double delta = (dx[n_data_pts - 1] + 0.01 - dx[0]) / (n_breakpoints - 1);
	double cur_brkpt = dx[0];
	for (i = 0; i < (unsigned int)n_breakpoints; i++)
	{
		gsl_vector_set(breakpoints, i, cur_brkpt);
		cur_brkpt += delta;
	} // i loop.

	/* Use uniform breakpoints on [0, 15] */
	//gsl_bspline_knots_uniform(dx[0], dx[n_data_pts - 1], bw);
	gsl_bspline_knots(breakpoints, bw);

	// https://lists.gnu.org/archive/html/help-gsl/2012-01/msg00011.html
	if (__DUMP_POLYFIT_MESSAGES__)
	{
		gsl_vector_fprintf(stdout, bw->knots, "%f ");
	}

	/* construct the fit matrix X */
	FILE* f_splines = NULL;
	if (__DUMP_POLYFIT_MESSAGES__)
	{
		f_splines = fopen("splines.txt", "a");
	}

	for (i = 0; i < (unsigned int)n_data_pts; ++i)
	{
		double xi = gsl_vector_get(x, i);

		/* compute B_j(xi) for all j */
		gsl_bspline_eval(xi, B, bw);

		/* fill in row i of X */
		for (j = 0; j < ncoeffs; ++j)
		{
			double Bj = gsl_vector_get(B, j);
			gsl_matrix_set(X, i, j, Bj);

			if (__DUMP_POLYFIT_MESSAGES__)
			{
				fprintf(f_splines, "%lf\t", Bj);
			}
		}

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(f_splines, "\n");
		}
	} // i loop.

	if (__DUMP_POLYFIT_MESSAGES__)
	{
		fclose(f_splines);
	}

	/* do the fit */
	gsl_multifit_linear(X, y, c, cov, &chisq, mw);

	//dof = n_data_pts - ncoeffs;
	//tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
	//Rsq = 1.0 - chisq / tss;

	//fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", 
	//                 chisq / dof, Rsq);

	// Reset the covariance matrix of the coefficients.
	for (unsigned int i = 0; i < ncoeffs; i++)
	{
		for (unsigned int j = 0; j < ncoeffs; j++)
		{
			gsl_matrix_set(cov, i, j, 0);
		}
	}

	/* output the smoothed curve */
	{
		double xi, yi, yerr;

		// Compute over the whole range.
		for (int i = 0; i < n_data_pts; i++)
		{
			xi = dx[i];
			gsl_bspline_eval(xi, B, bw);
			gsl_multifit_linear_est(B, c, cov, &yi, &yerr);

			if (gsl_finite(yi))
			{
				reconst_y[i] = yi;
			}
			else
			{
				reconst_y[i] = 0;
			}
			//printf("%f %f\n", xi, yi);
		}
	}

	//gsl_rng_free(r);
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	//gsl_vector_free(w);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);

	return true;

}

bool bsplinefit(int n_data_pts, 
				int n_spline_coeffs, 
				int spline_order,
				double *dx, double *dy, 
				double* reconst_y, 
				double *store)
{
	if(__DUMP_POLYFIT_MESSAGES__)
		fprintf(stderr, "Fitting %d points with B-spline with %d coeffs of order %d (%d points per fit).\n", n_data_pts, n_spline_coeffs, spline_order, (int)((dx[n_data_pts-1]-dx[0]) / n_spline_coeffs));

	/* number of data points to fit */
	//#define N        200

	/* number of fit coefficients */
	//#define NCOEFFS  12

	/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
	//#define NBREAK   (NCOEFFS - 2)

  //const size_t n = n_data_pts;
  const size_t ncoeffs = n_spline_coeffs;

  // This is the number of knots.
  const size_t nbreak = ncoeffs - (spline_order - 2);
  size_t i, j;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  //double dy;
  gsl_vector *c;
  //gsl_vector *w;
  gsl_vector *x, *y;
  gsl_matrix *X, *cov;
  gsl_multifit_linear_workspace *mw;
  double chisq;
  //double Rsq, dof, tss;

  //gsl_rng_env_setup();
  //r = gsl_rng_alloc(gsl_rng_default);

  /* allocate a cubic bspline workspace (k = 4) */
  bw = gsl_bspline_alloc(spline_order, nbreak);
  B = gsl_vector_alloc(ncoeffs);

  x = gsl_vector_alloc(n_data_pts);
  y = gsl_vector_alloc(n_data_pts);
  X = gsl_matrix_alloc(n_data_pts, ncoeffs);
  c = gsl_vector_alloc(ncoeffs);
  //w = gsl_vector_alloc(n_data_pts);
  cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
  mw = gsl_multifit_linear_alloc(n_data_pts, ncoeffs);

  //printf("#m=0,S=0\n");
  /* this is the data to be fitted */
  for (i = 0; i < (unsigned int)n_data_pts; ++i)
    {
      gsl_vector_set(x, i, dx[i]);
      gsl_vector_set(y, i, dy[i]);
    }

  /* use uniform breakpoints on [0, 15] */
  gsl_bspline_knots_uniform(dx[0], dx[n_data_pts-1], bw);

  // https://lists.gnu.org/archive/html/help-gsl/2012-01/msg00011.html
  if (__DUMP_POLYFIT_MESSAGES__)
  {
	  gsl_vector_fprintf(stdout, bw->knots, "%f ");
  }

  /* construct the fit matrix X */
	FILE* f_splines = NULL;
	if (__DUMP_POLYFIT_MESSAGES__)
	{
		f_splines = fopen("splines.txt", "a");
	}

	for (i = 0; i < (unsigned int)n_data_pts; ++i)
	{
		double xi = gsl_vector_get(x, i);

		/* compute B_j(xi) for all j */
		gsl_bspline_eval(xi, B, bw);

		/* fill in row i of X */
		for (j = 0; j < ncoeffs; ++j)
		{
			double Bj = gsl_vector_get(B, j);
			gsl_matrix_set(X, i, j, Bj);

			if (__DUMP_POLYFIT_MESSAGES__)
			{
				fprintf(f_splines, "%lf\t", Bj);
			}
		}

		if (__DUMP_POLYFIT_MESSAGES__)
		{
			fprintf(f_splines, "\n");
		}
	} // i loop.

	if (__DUMP_POLYFIT_MESSAGES__)
	{
		fclose(f_splines);
	}

  /* do the fit */
  gsl_multifit_linear(X,  y, c, cov, &chisq, mw);

  //dof = n_data_pts - ncoeffs;
  //tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
  //Rsq = 1.0 - chisq / tss;

  //fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", 
  //                 chisq / dof, Rsq);

  // Reset the covariance matrix of the coefficients.
	for(unsigned int i = 0; i < ncoeffs; i++)
	{
		   for(unsigned int j = 0; j < ncoeffs; j++)
		   {
				   gsl_matrix_set(cov, i, j, 0);
		   }
	}

  /* output the smoothed curve */
  {
    double xi, yi, yerr;

	// Compute over the whole range.
	for(int i = 0; i < n_data_pts; i++)
	{
		xi = dx[i];
        gsl_bspline_eval(xi, B, bw);
        gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
		
		if (gsl_finite(yi))
		{
			reconst_y[i] = yi;
		}
		else
		{
			reconst_y[i] = 0;
		}
        //printf("%f %f\n", xi, yi);
	}
  }

  //gsl_rng_free(r);
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  //gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);

  return true;
	
}


