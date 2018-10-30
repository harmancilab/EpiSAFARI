#ifndef __EPISFR_UTILS__
#define __EPISFR_UTILS__

#include <vector>
#include <algorithm>
using namespace std;

struct t_extrema_statistic_defition;
struct t_annot_region;

void get_summit_dip_posns_per_ER(double* signal_profile, int l_profile,
	double* multi_mapp_signal, int l_multi_map_signal, double max_multi_mapp_val,
	int ER_start, int ER_end,
	int l_trough_win,
	double minimum_summit_height_2_max_summit_ratio,
	double maximum_dip_height_2_max_summit_ratio,
	vector<int>* selected_minima_posns,
	vector<int>* selected_maxima_posns);

struct t_extrema_statistic_defition
{
	double max_signal_at_trough;
	double min_summit2trough_ratio_per_trough;
	double max_mapp_signal_at_trough;

	double min_signal_at_summit;	

	int max_summit2trough_dist_in_bp;

	double max_multimapp_signal_at_trough;
};

struct t_episfr_annot_info
{
	char* element_type;
	char* element_name;
};

char* load_binary_sequence_file(char* bin_seq_fp, int& l_seq);

void binarize_fasta_file(char* fasta_fp, char* bin_dir);

vector<t_annot_region*>* load_annotation(char* gff_fp, int l_half_prom);

vector<t_annot_region*>* get_significant_extrema_per_signal_profile(const char* op_dir, 
	const char* chr_id, 
	double* signal_profile, int l_profile,
	double* multimapp_signal_profile, int l_multimapp_profile,
	char* chrom_seq,
	t_extrema_statistic_defition* extrema_statistic_defn);

void annotate_features(char* signal_directory, char* gff_fp, int l_half_prom, char* op_fp);


void select_points_of_interest_per_RD_signal_profile(double* signal_profile,
	int start_i, int end_i,
	int max_dist_between_cons_pts,
	vector<double>* x_vec, vector<double>* y_vec,
	bool sparse_signal);

void bspline_encode_mapped_read_profile(char* signal_dir,
	char* chr_id,
	int l_frag,
	int n_spline_coeff,
	int bspline_order,
	int min_n_pts_2_encode,
	int max_dist_between_cons_pts,
	double max_max_err,
	double max_avg_err,
	int l_win,
	bool sparse_profile,
	double top_err_perc_frac,
	int l_med_filt_win);

#endif // __EPISFR_UTILS__