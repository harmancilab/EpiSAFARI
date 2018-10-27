#ifndef __SIGNAL_TRACK_FILE_INTERFACE__
#define __SIGNAL_TRACK_FILE_INTERFACE__

#include <vector>

using namespace std;

struct t_annot_region;

double* load_per_nucleotide_binary_profile(char* binary_per_nucleotide_profile_fp, int& l_profile);
void dump_per_nucleotide_binary_profile_per_bedgraph(const char* bgr_fp, bool dump_binary, const char* op_dir);
void dump_per_nucleotide_binary_profile(double* profile, int l_profile, const char* op_fp);
void dump_bedGraph_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, const char* chrom, const char* op_fp);

double* load_per_nucleotide_BGR_track(char* bgr_fp, int& l_profile);

unsigned char* load_per_nucleotide_binary_uchar_profile(char* binary_per_nucleotide_profile_fp, int& l_profile);
void dump_per_nucleotide_uchar_binary_profile(unsigned char* signal_profile, int l_profile, char* op_fp);

void floorize_profile(double* signal_profile, int l_profile);
double* get_zero_profile(int l_profile);

#endif // __SIGNAL_TRACK_FILE_INTERFACE__

