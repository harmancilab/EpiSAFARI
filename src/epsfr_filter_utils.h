#ifndef __GSL_FFT_FILTER_UTILS__
#define __GSL_FFT_FILTER_UTILS__

#include <vector>
using namespace std;

double* median_filter_data(double* track_data,
	int l_track_data,
	int l_averaging_win,
	int skip_value);

#endif // __GSL_FFT_FILTER_UTILS__