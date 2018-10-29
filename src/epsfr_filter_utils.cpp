#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "epsfr_filter_utils.h"
#include "epsfr_utils.h"
#include "epsfr_ansi_string.h"
#include "epsfr_genomics_coords.h"
#include "epsfr_annot_region_tools.h"
#include "epsfr_signal_track_tools.h"
#include "epsfr_min_max_utils.h"
#include <algorithm>
#include "string.h"

using namespace std;

bool __DUMP_FILTER_MSGS__ = false;

#ifdef __unix__
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_eigen.h>
	#include <gsl/gsl_sort.h>
	//#include <gsl/gsl_wavelet.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_errno.h>
	#include <gsl/gsl_fft_complex.h>

	#define REAL(z,i) ((z)[2*(i)])
	#define IMAG(z,i) ((z)[2*(i)+1])
#endif

#ifdef __unix__
gsl_vector* get_gsl_vector_per_array(double* data, int l_data)
{
	// Allocate the block.
	gsl_block* vec_block = new gsl_block();
	vec_block->size = l_data;
	vec_block->data = data;

	// Allocate the gsl_vector, then copy the data.
	gsl_vector* vec = new gsl_vector();
	vec->size = l_data;
	vec->data = data;
	vec->block = vec_block;
	vec->owner = 0;
	vec->stride = 1;

	return(vec);
}
#endif


double* median_filter_data(double* track_data,
	int l_track_data,
	int l_averaging_win,
	int skip_value) // This is the value to be skipped from the window values.
{
	// Do copying from start to end.
	double* cur_filtered_track = new double[l_track_data + 2];

	int n_signal_wins = 0;
	int half_l_averaging = l_averaging_win / 2;

	// Set the maximum value for the histogram.
	int MAX_VAL = -1000 * 1000;
	int MIN_VAL = 1000 * 1000;
	for (int i_sig = 1; i_sig <= l_track_data; i_sig++)
	{
		if (MAX_VAL < (int)(track_data[i_sig]))
		{
			MAX_VAL = (int)(track_data[i_sig]);
		}

		if (MIN_VAL >(int)(track_data[i_sig]))
		{
			MIN_VAL = (int)(track_data[i_sig]);
		}
	} // i_sig loop.

	  // Make sure that this is the maximum.
	MAX_VAL += 1000;
	MIN_VAL -= 1000;

	// Initialize the current pdf.
	int* cur_win_pdf = NULL;
	int cur_win_max = 0;
	int cur_win_min = 1000 * 1000;
	double* cur_win_vals = new double[l_averaging_win + 2];
	int prev_avg_start = 0;
	int prev_avg_end = 0;

	// Go over all the positions as the middle of the filtering window.
	for (int cur_win_mid = 1; cur_win_mid <= l_track_data; cur_win_mid++)
	{
		int cur_avg_start = (cur_win_mid > half_l_averaging) ? (cur_win_mid - half_l_averaging) : (1);
		int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data) ? (cur_win_mid + half_l_averaging) : (l_track_data);
		if (cur_win_pdf == NULL)
		{
			cur_win_pdf = new int[MAX_VAL - MIN_VAL + 2];
			memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
			//t_string::set_byte_buffer(cur_win_pdf, sizeof(int) * (MAX_VAL - MIN_VAL + 1), 0);
			cur_win_pdf -= MIN_VAL;

			// Generate the pdf, get the minimum and maximum in the current window.
			for (int i = cur_avg_start; i <= cur_avg_end; i++)
			{
				if (track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]++;

					if (cur_win_max < track_data[i])
					{
						cur_win_max = track_data[i];
					}

					if (cur_win_min > track_data[i])
					{
						cur_win_min = track_data[i];
					}
				} // skip_Value check.
			} // i loop.
		} // cur_win_pdf is NULL check.
		else
		{
			// Remove the old values from the pdf, add the new values.
			for (int i = prev_avg_start; i < cur_avg_start; i++)
			{
				if (track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]--;
				}
			} // i loop.

			  // Update the min and max only for the values that are new in this window.
			for (int i = prev_avg_end + 1; i <= cur_avg_end; i++)
			{
				if (track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]++;

					if (cur_win_max < track_data[i])
					{
						cur_win_max = (int)track_data[i];
					}

					if (cur_win_min > track_data[i])
					{
						cur_win_min = (int)track_data[i];
					}
				} // skip_value check.
			}

			// Sanity Check: The total # of points must be equal to the window length.
			int n_total_pts = 0;
			for (int i = cur_win_min; i <= cur_win_max; i++)
			{
				if (i != skip_value)
				{
					n_total_pts += cur_win_pdf[i];
				}
			} // i loop.
		} // cur_win_pdf is NULL check.

		  // Count the total number of points without the skip value.
		int n_total = 0;
		for (int i = cur_win_min; i <= cur_win_max; i++)
		{
			if (i != skip_value)
			{
				n_total += cur_win_pdf[i];
			}
		} // i loop.

		  // Generate the window median.
		int cur_win_median = 0;
		int n_cur_total = 0;
		for (int i = cur_win_min; i <= cur_win_max; i++)
		{
			if (i != skip_value)
			{
				n_cur_total += cur_win_pdf[i];
				if (n_cur_total > n_total / 2)
				{
					// We found the median, can break out of the loop.
					cur_win_median = i;
					break;
				}
			}
		} // i loop.

		  // Track the minimum and maximum from the histogram to update them for the current window.
		int updated_win_min = 0;
		for (int i = cur_win_min; i <= cur_win_max; i++)
		{
			if (i != skip_value)
			{
				if (cur_win_pdf[i] > 0)
				{
					updated_win_min = i;
					break;
				}
			} // skip_value check.
		} // i loop.

		int updated_win_max = 0;
		for (int i = cur_win_max; i >= cur_win_min; i--)
		{
			if (i != skip_value)
			{
				if (cur_win_pdf[i] > 0)
				{
					updated_win_max = i;
					break;
				}
			} // skip_value check.
		} // i loop.

		  // Set the median.
		  //int median_per_pdf = cur_win_median;
		cur_filtered_track[cur_win_mid] = cur_win_median;

		// Update the previous averaging window limits.
		prev_avg_start = cur_avg_start;
		prev_avg_end = cur_avg_end;
		cur_win_min = updated_win_min;
		cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
		// Get the median via qsort and compare as a sanity check.
		int n_valid_vals = 0;
		for (int i = cur_avg_start; i <= cur_avg_end; i++)
		{
			if (track_data[i] != skip_value)
			{
				cur_win_vals[i - cur_avg_start] = track_data[i];
				n_valid_vals++;
			}
		} // i loop.
		qsort(cur_win_vals, n_valid_vals, sizeof(double), sort_doubles_descending);
		//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

		if (cur_win_vals[n_valid_vals / 2] != median_per_pdf)
		{
			fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[n_valid_vals / 2], median_per_pdf);
			for (int i = cur_avg_start; i <= cur_avg_end; i++)
			{
				fprintf(stderr, "%lf ", track_data[i]);
			} // i loop.
			fprintf(stderr, "\n");
			getc(stdin);
		}
#endif // _QSORT_CHECK_

		n_signal_wins++;
	} // main signal filtering loop.

	  // Free pdf memory.
	delete[] cur_win_vals;
	delete[](cur_win_pdf + MIN_VAL);

	return(cur_filtered_track);
}