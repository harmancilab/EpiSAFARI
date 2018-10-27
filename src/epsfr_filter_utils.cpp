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
