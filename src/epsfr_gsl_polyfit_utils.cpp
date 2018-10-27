#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include "epsfr_gsl_polyfit_utils.h"

bool __DUMP_POLYFIT_MESSAGES__ = false;

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

  /* construct the fit matrix X */
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
        }
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
		reconst_y[i] = yi;
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

// Select changing points to use in the polyfit.
void select_points_of_interest_per_RD_signal_profile(double* signal_profile, 
													int start_i, int end_i, 
													int max_dist_between_cons_pts, 
													vector<double>* x_vec, vector<double>* y_vec)
{
	//fprintf(stderr, "Selecting POI for signal track in [%d-%d] with maximum distance of points %d\n", start_i, end_i, max_dist_between_cons_pts);

	// Update the points at each change in signal.
	double cur_y = signal_profile[start_i];
	double cur_x = start_i;
	x_vec->push_back(start_i);
	y_vec->push_back(cur_y);
	for(int i = start_i; i <= end_i; i++)
	{
		if(cur_y != signal_profile[i] ||
			cur_x + max_dist_between_cons_pts <= i)
		{
			// Update the current point.
			cur_y = signal_profile[i];
			cur_x = i;

			x_vec->push_back(i);
			y_vec->push_back(signal_profile[i]);
		}
		else
		{
			// Do nothing: Do not update the data vectors.
		}
	} // i loop.
}