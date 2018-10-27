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

vector<t_annot_region*>* load_annotation(char* gff_fp, int l_half_prom);

vector<t_annot_region*>* get_significant_extrema_per_signal_profile(const char* op_dir, const char* chr_id, double* signal_profile, int l_profile,
	double* multimapp_signal_profile, int l_multimapp_profile,
	t_extrema_statistic_defition* extrema_statistic_defn);

void annotate_features(char* signal_directory, char* gff_fp, int l_half_prom, char* op_fp);

#endif // __EPISFR_UTILS__