#ifndef __EPISFR_UTILS__
#define __EPISFR_UTILS__

#include <vector>
#include <algorithm>
using namespace std;

struct t_extrema_statistic_defition;
struct t_annot_region;

enum {
	VALLEY_SIGNIFICANCE_BINOMIAL_INTERSECTED_NULLS,
	VALLEY_SIGNIFICANCE_BINOMIAL_UNION_NULLS,
	VALLEY_SIGNIFICANCE_MULTINOMIAL,
};

enum {
	HEIGHT_BASED_HILL_SCORE,
	DIST_BASED_HILL_SCORE
};

struct t_valley_significance_info
{
	double log_p_val;
	double log_q_val;
};

double atof_null(char* str);
int atoi_null(char* str);

double* load_signal_covg_per_directory_chr_id(const char* dat_dir,
	char* chr_id,
	int l_fragment,
	int& l_loaded_covg,
	bool& reads_loaded);

double* load_signal_covg_per_signal_file(const char* cur_dat_fp,
	int l_fragment,
	int& l_loaded_covg,
	bool& reads_loaded);

void get_2_sample_differential_valleys(char* sample1_valleys_bed_fp, char* sample1_dir, char* sample2_valleys_bed_fp, char* sample2_dir, t_extrema_statistic_defition* extrema_statistic_defn);

double* buffer_log_factorials(int n);

double get_valley_significance_per_binomial_test(int p_val_type, double* signal_profile, int l_profile, int left_max_posn, int right_max_posn, int min_posn, int l_normalizer, int scaling_factor, double* _log_factorials);

double get_valley_significance_per_multinomial_test(double* signal_profile, int l_profile, int left_max_posn, int right_max_posn, int min_posn, int l_normalizer, int scaling_factor, double* _log_factorials);
double get_multinomial_pval_per_extrema_vals(double norm_min_vic_sig, double norm_left_max_vic_sig, double norm_right_max_vic_sig, double* log_factorials);
double get_binomial_pval_per_extrema_vals(int p_val_type, double norm_min_vic_sig, double norm_left_max_vic_sig, double norm_right_max_vic_sig, double* log_factorials);

void append_signal_2_regions(char* signal_directory, int l_frag, char* bed_fp, char* op_fp);

void uniformize_dip_relative_window_positions(char* valleys_BED_fp, int l_win, int n_bins_per_win, char* op_fp);

void get_benjamini_hochberg_corrected_p_values_per_valleys(vector<t_annot_region*>* valleys);

double pearson_correlation(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins);
double slope(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins);
double PeakSeq_slope(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins);
void get_per_window_scaling_factors_for_profile1_per_profile2(double* signal_profile1, int l_prof1, double* signal_profile2, int l_prof2, int l_win,
	double& per_win_2DOF_lls_scaling_factor,
	double& per_win_1DOF_lls_scaling_factor,
	double& total_sig_scaling_factor);

void dump_valleys(vector<t_annot_region*>* significant_valleys, char* valleys_bed_fp);

vector<t_annot_region*>* merge_overlapping_valleys_per_pval_minimization(char* valleys_BED_fp, int l_minima_vicinity);

void refine_valleys_per_pval_min(double* signal_profile, int l_profile, vector<t_annot_region*>* valley_regs, int l_vic_expand);
void refine_valleys_per_height_trim(double* signal_profile, int l_profile, vector<t_annot_region*>* valley_regs, double height_frac);

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
	double min_signal_at_summit;
	
	double min_summit2trough_ratio_per_trough;
	
	int min_summit2trough_dist_in_bp;
	int max_summit2trough_dist_in_bp;
	
	double max_multimapp_signal_at_trough;

	double log_q_val_threshold;

	int p_val_estimate_extrema_vic_window_length;
	double p_val_estimate_signal_scaling_factor;

	bool sparse_profile;

	char p_val_type;

	int l_minima_vicinity_per_merging;

	int hill_score_type;
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
	unsigned char* multimapp_signal_profile, int l_multimapp_profile,
	char* chrom_seq,
	t_extrema_statistic_defition* extrema_statistic_defn);

void annotate_features(char* valleys_bed_fp, char* gff_fp, int l_half_prom, char* op_fp);

void select_points_of_interest_per_RD_signal_profile(double* signal_profile, int l_profile,
	int start_i, int end_i,
	int max_dist_between_cons_pts,
	vector<double>* x_vec, vector<double>* y_vec,
	bool sparse_signal);

void bspline_encode_mapped_read_profile(char* signal_dir,
	char* chr_id,
	int l_frag,
	int n_spline_coeff,
	int bspline_order,
	char brkpts_type,
	int min_n_pts_2_encode,
	int max_dist_between_cons_pts,
	double max_top_err,
	double max_avg_err,
	int l_win,
	bool sparse_profile,
	double top_err_perc_frac,
	int l_step_win,
	int l_med_filt_win);

void bspline_encode_mapped_read_profile_expanded_windows(char* signal_dir,
	char* chr_id,
	int l_frag,
	int n_spline_coeff,
	int bspline_order,
	char brkpts_type,
	int min_n_pts_2_encode,
	int max_dist_between_cons_pts,
	double max_top_err,
	double max_avg_err,
	int l_win,
	bool sparse_profile,
	double top_err_perc_frac,
	int l_med_filt_win);

#endif // __EPISFR_UTILS__