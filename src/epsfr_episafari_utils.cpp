#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "epsfr_episafari_utils.h"
#include "epsfr_xlog_math.h"
#include "epsfr_rng.h"
#include "epsfr_seed_manager.h"
#include "epsfr_nucleotide.h"
#include "epsfr_genomics_coords.h"
#include "epsfr_nomenclature.h"
#include "epsfr_ansi_string.h"
#include "epsfr_min_max_utils.h"
#include "epsfr_filter_utils.h"
#include "epsfr_nomenclature.h"
#include "epsfr_utils.h"
#include "epsfr_annot_region_tools.h"
#include "epsfr_mapped_read_tools.h"
#include "epsfr_signal_track_tools.h"
#include "epsfr_gsl_polyfit_utils.h"
#include "ctype.h"

#ifdef __unix__
	#include <unistd.h> 
#endif

bool __DUMP_EPISAFARI_UTILS_MESSAGES__ = false;

int atoi_null(char* str)
{
	if (str == NULL)
	{
		return(-1);
	}
	else
	{
		return(atoi(str));
	}
}

double atof_null(char* str)
{
	if (str == NULL)
	{
		return(-1);
	}
	else
	{
		return(atof(str));
	}
}

void uniformize_dip_relative_window_positions(char* valleys_BED_fp, int l_win, int n_bins_per_win, char* op_fp)
{
	fprintf(stderr, "Uniformizing valleys in %s over %d long windows using %d bins on the histogram.\n", valleys_BED_fp, l_win, n_bins_per_win);

	vector<t_annot_region*>* valley_regs = load_BED_with_line_information(valleys_BED_fp);
	fprintf(stderr, "Loaded %d valley regions.\n", (int)valley_regs->size());

	FILE* f_bed = open_f(valleys_BED_fp, "r");
	char* header_line = getline(f_bed);
	if (header_line[0] != '#')
	{
		fprintf(stderr, "Could not read header line.\n");
		exit(0);
	}
	close_f(f_bed, valleys_BED_fp);

	t_string_tokens* header_line_toks = t_string::tokenize_by_chars(header_line, "\t");
	vector<char*>* header_line_cols = t_string::copy_tokens_2_strs(header_line_toks);	

	fprintf(stderr, "%d tokens in the header: %s\n", (int)header_line_cols->size(), header_line);

	//for (int i_col = 0; i_col < header_line_cols->size(); i_col++)
	//{
	//	fprintf(stderr, "%d: %s (%s)\n", i_col, header_line_cols->at(i_col), header_line_toks->at(i_col)->str());
	//} // i_col loop.

	t_string::clean_tokens(header_line_toks);

	int MIN_POSN_col_i = t_string::get_i_str_ci(header_line_cols, "MIN_POSN");
	fprintf(stderr, "Found minimum position column @ %d\n", MIN_POSN_col_i);

	int qVal_col_i = t_string::get_i_str_ci(header_line_cols, "qVal");
	fprintf(stderr, "Found q-value column @ %d\n", qVal_col_i);	

	fprintf(stderr, "Parsing minima and scores.\n");
	for (int i_reg = 0; i_reg < (int)valley_regs->size(); i_reg++)
	{
		char* cur_reg_line = (char*)(valley_regs->at(i_reg)->data);
		t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");
		double cur_valley_pval = atof(toks->at(qVal_col_i)->str());
		int min_posn = atoi(toks->at(MIN_POSN_col_i)->str());
		t_string::clean_tokens(toks);

		valley_regs->at(i_reg)->score = min_posn % l_win;
		valley_regs->at(i_reg)->dbl_score = cur_valley_pval;
	} // i_reg loop.

	fprintf(stderr, "Generating bins.\n");
	vector<t_annot_region*>** per_bin_regions = new vector<t_annot_region*>*[n_bins_per_win + 1];
	int l_bin = l_win / n_bins_per_win;
	fprintf(stderr, "l_bin=%d\n", l_bin);
	for (int i_bin = 0; i_bin < n_bins_per_win; i_bin++)
	{
		per_bin_regions[i_bin] = new vector<t_annot_region*>();
	} // i_bin loop.

	fprintf(stderr, "Assigning %d valleys to bins.\n", (int)valley_regs->size());
	for (int i_reg = 0; i_reg < (int)valley_regs->size(); i_reg++)
	{
		int i_bin = (valley_regs->at(i_reg)->score / l_bin);
		per_bin_regions[i_bin]->push_back(valley_regs->at(i_reg));

		//fprintf(stderr, "%s: bin_i: %d, minima_rel_posn: %d (%.4f)\n", (char*)(valley_regs->at(i_reg)->data), i_bin, valley_regs->at(i_reg)->score, valley_regs->at(i_reg)->dbl_score);
	} // i_reg loop.

	// For each bin, sort the valleys and save.
	int n_valleys_per_bin_2_dump = (int)(valley_regs->size());
	for (int i_bin = 0; i_bin < n_bins_per_win; i_bin++)
	{
		n_valleys_per_bin_2_dump = MIN((int)(per_bin_regions[i_bin]->size()), n_valleys_per_bin_2_dump);
		sort(per_bin_regions[i_bin]->begin(), per_bin_regions[i_bin]->end(), sort_regions_per_dbl_score);
	} // i_bin loop.

	fprintf(stderr, "Saving top %d valleys per bin.\n", n_valleys_per_bin_2_dump);

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "%s\n", header_line);
	for (int i_bin = 0; i_bin < n_bins_per_win; i_bin++)
	{
		for (int i_valley = 0; i_valley < n_valleys_per_bin_2_dump; i_valley++)
		{
			char* cur_reg_line = (char*)(per_bin_regions[i_bin]->at(i_valley)->data);
			fprintf(f_op, "%s\n", cur_reg_line);
		} // i_valley loop
	} // i_bin loop.
	fclose(f_op);
}

void append_signal_2_regions(char* signal_directory, int l_frag, char* bed_fp, char* op_fp)
{
	fprintf(stderr, "Appending signal from %s (l_frag=%d) to regions in %s\n", signal_directory, l_frag, bed_fp);
	vector<t_annot_region*>* regs = load_BED_with_line_information(bed_fp);
	fprintf(stderr, "Loaded %d regions.\n", (int)regs->size());

	t_restr_annot_region_list* restr_regs = restructure_annot_regions(regs);

	double total_signal = 0;
	for (int i_chr = 0; i_chr < (int)restr_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s.\n", restr_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_regs = restr_regs->regions_per_chrom[i_chr];

		vector<t_annot_region*>* cur_chr_merged_regs = merge_annot_regions(cur_chr_regs, 1);

		bool reads_loaded = false;
		int l_loaded_sig = 0;
		double* signal_profile = load_signal_covg_per_directory_chr_id(signal_directory, restr_regs->chr_ids->at(i_chr), l_frag, l_loaded_sig, reads_loaded);

		//// Update the total signal.
		//for (int i = 0; i < l_loaded_sig; i++)
		//{
		//	total_signal += signal_profile[i];
		//} // i loop.
		for (int i_reg = 0; i_reg < (int)cur_chr_merged_regs->size(); i_reg++)
		{
			for (int i = cur_chr_merged_regs->at(i_reg)->start; i <= cur_chr_merged_regs->at(i_reg)->end; i++)
			{
				if (i < l_loaded_sig)
				{
					total_signal += signal_profile[i];
				}
			} // i loop.
		} // i_reg loop.

		delete_annot_regions(cur_chr_merged_regs);

		for (int i_reg = 0; i_reg < (int)cur_chr_regs->size(); i_reg++)
		{
			double cur_reg_sig = 0;
			for (int i = cur_chr_regs->at(i_reg)->start; i <= cur_chr_regs->at(i_reg)->end; i++)
			{
				if (i < l_loaded_sig)
				{
					cur_reg_sig += signal_profile[i];
				}
			} // i loop.

			cur_chr_regs->at(i_reg)->dbl_score = cur_reg_sig;
		} // i_reg loop.

		delete[] signal_profile;
	} // i_chr loop.

	fprintf(stderr, "Total signal: %.2f\n", total_signal);

	// Save.
	FILE* f_reg = open_f(bed_fp, "r");
	char* header_line = getline(f_reg);
	fclose(f_reg);

	// Set the column id and add it to the header.
	char* sig_col_id = t_string::copy_me_str(signal_directory);
	t_string::replace_avoid_list(sig_col_id, "/.", '_');	

	FILE* f_op = open_f(op_fp, "w");
	if (header_line[0] == '#')
	{
		fprintf(f_op, "%s\t%s\n", header_line, sig_col_id);
	}

	// Save the signal for each region.
	for (int i_reg = 0; i_reg < (int)regs->size(); i_reg++)
	{
		double cur_reg_norm_sig = 1000 * 1000 * regs->at(i_reg)->dbl_score / total_signal;
		fprintf(f_op, "%s\t%.3f\n", (char*)(regs->at(i_reg)->data), cur_reg_norm_sig);
	} // i_reg loop.
	fclose(f_op);
}

void get_per_window_scaling_factors_for_profile1_per_profile2(double* signal_profile1, int l_prof1, double* signal_profile2, int l_prof2, int l_win,
	double& per_win_2DOF_lls_scaling_factor,
	double& per_win_1DOF_lls_scaling_factor,
	double& total_sig_scaling_factor)
{
	/*vector<int>* chip_seq_frag_cnt = new vector<int>();
	vector<int>* control_frag_cnt = new vector<int>();*/

	int l_processable_prof = (l_prof1 < l_prof2) ? (l_prof1) : (l_prof2);
	int n_wins_2_process = l_processable_prof / l_win;

	double* prof1_win_total_sigs = new double[n_wins_2_process + 2];
	double* prof2_win_total_sigs = new double[n_wins_2_process + 2];

	double total_prof1 = 0;
	double total_prof2 = 0;
	for (int i = 1; i <= l_prof1; i++)
	{
		total_prof1 += signal_profile1[i];
	} // i loop

	for (int i = 1; i <= l_prof2; i++)
	{
		total_prof2 += signal_profile2[i];
	} // i loop

	total_sig_scaling_factor = total_prof2 / total_prof1;

	// Count the numbers per window. 
	//int n_scaling_win = n_meg_wins * MEG_BASE / win_size;

	//int chip_seq_background_cnt = 0;
	//int control_cnt = 0;

	// Exploit the fact that the fragments are sorted with respect to base indices in the lists:
	// Keep track of the current fragment index for both fragment lists.
	//int cur_chr_fragment_list_base_i = 0;
	//int cur_chr_control_fragment_list_base_i = 0;

	int n_processed_wins = 0;

	FILE* f_per_win_vals = NULL;
	if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
		f_per_win_vals = fopen("per_win_vals.txt", "w");
	for (int i_win = 0; i_win < n_wins_2_process; i_win++)
	{
		int cur_win_start = l_win * i_win + 1;
		int cur_win_end = (i_win + 1) * l_win;

		double cur_total_prof1 = 0;
		double cur_total_prof2 = 0;
		for (int i_nuc = cur_win_start; i_nuc < cur_win_end; i_nuc++)
		{
			if (i_nuc < l_prof1 && i_nuc < l_prof2)
			{
				cur_total_prof1 += signal_profile1[i_nuc];
				cur_total_prof2 += signal_profile2[i_nuc];
			}
		} // i_nuc loop.

		if (cur_total_prof1 > 0 && cur_total_prof2 > 0)
		{
			prof1_win_total_sigs[n_processed_wins] = cur_total_prof1;
			prof2_win_total_sigs[n_processed_wins] = cur_total_prof2;
			n_processed_wins++;
		}

		if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
			fprintf(f_per_win_vals, "%lf\t%lf\n", cur_total_prof1, cur_total_prof2);
	} // i_win loop.

	if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
		fclose(f_per_win_vals);

	// Compute the slopes.
	per_win_2DOF_lls_scaling_factor = PeakSeq_slope(prof1_win_total_sigs, prof2_win_total_sigs, n_processed_wins);
	per_win_1DOF_lls_scaling_factor = slope(prof1_win_total_sigs, prof2_win_total_sigs, n_processed_wins);

	double r = pearson_correlation(prof1_win_total_sigs, prof2_win_total_sigs, n_processed_wins);

	if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
		fprintf(stderr, "Pearson Coeff: %lf\nSignal Division: %lf\n1-DOF: %lf\n2-DOF: %lf\n", r, total_sig_scaling_factor, per_win_1DOF_lls_scaling_factor,
			per_win_2DOF_lls_scaling_factor);

	//// If the total profile1 versus profile2 is much larger, use that.
	//if(total_prof1 / total_prof2 > per_win_1DOF_lls_scaling_factor)
	//{
	//	per_win_1DOF_lls_scaling_factor = total_prof1 / total_prof2;
	//}

	delete[] prof1_win_total_sigs;
	delete[] prof2_win_total_sigs;
}

//double PeakSeq_slope(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt, int n_processed_wins)
double PeakSeq_slope(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins)
{
	// Compute the slope from 2 DOF (slope) LLS regression.
	double num = 0.0f;
	double den = 0.0f;

	double sx = 0.0f;
	double sy = 0.0f;
	double sxx = 0.0f;
	double sxy = 0.0f;

	for (int i = 0; i < n_processed_wins; i++)
	{
		sx += prof1_total_sigs_per_win[i];
		sy += prof2_total_sigs_per_win[i];
		sxx += prof1_total_sigs_per_win[i] * prof1_total_sigs_per_win[i];
		sxy += prof1_total_sigs_per_win[i] * prof2_total_sigs_per_win[i];

		//sx += control_frag_cnt->at(i);
		//sy += chip_seq_frag_cnt->at(i);
		//sxx += control_frag_cnt->at(i) * control_frag_cnt->at(i);
		//sxy += control_frag_cnt->at(i) * chip_seq_frag_cnt->at(i);
	} // i loop

	num = n_processed_wins * sxy - sx * sy;
	den = n_processed_wins * sxx - sx * sx;

	if (num == 0 || den == 0)
	{
		fprintf(stderr, "Normalization(1): Very low signal mapped.\n");
		return(1.0);
	}

	return(num / den);
}

double slope(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins)
{
	// Compute the slope from 1 DOF (slope) LLS regression.
	double num = 0.0f;
	double den = 0.0f;
	for (int i = 0; i < n_processed_wins; i++)
	{
		num += prof1_total_sigs_per_win[i] * prof2_total_sigs_per_win[i];
		den += prof1_total_sigs_per_win[i] * prof1_total_sigs_per_win[i];
		//num += control_frag_cnt->at(i) * chip_seq_frag_cnt->at(i);
		//den += control_frag_cnt->at(i) * control_frag_cnt->at(i);
	} // i loop

	if (num == 0 || den == 0)
	{
		fprintf(stderr, "Normalization(2): Very low signal mapped.\n");
		return(1.0);
	}

	return(num / den);
}


// Copied from http://www.vias.org/tmdatanaleng/cc_corr_coeff.html
// Also check http://www.statsdirect.com/help/regression_and_correlation/sreg.htm
//double pearson_correlation(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt)
double pearson_correlation(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins)
{
	double num = 0.0f;
	//double den = 0.0f;

	double prof1_mean = 0.0f;
	double prof2_mean = 0.0f;
	for (int i_win = 0; i_win < n_processed_wins; i_win++)
	{
		//control_mean += control_frag_cnt->at(i_frag);
		prof1_mean += prof1_total_sigs_per_win[i_win];
		prof2_mean += prof2_total_sigs_per_win[i_win];
	} // i_win loop.
	prof1_mean /= n_processed_wins;
	prof2_mean /= n_processed_wins;

	//for(int i_frag = 0; i_frag < chip_seq_frag_cnt->size(); i_frag++)
	for (int i_win = 0; i_win < n_processed_wins; i_win++)
	{
		num += (prof1_total_sigs_per_win[i_win] - prof1_mean) * (prof2_total_sigs_per_win[i_win] - prof2_mean);
	} // i_win loop.

	double prof1_var = 0.0f;
	for (int i_win = 0; i_win < n_processed_wins; i_win++)
	{
		prof1_var += (prof1_total_sigs_per_win[i_win] - prof1_mean) * (prof1_total_sigs_per_win[i_win] - prof1_mean);
	}
	prof1_var = pow(prof1_var, .5);

	double prof2_var = 0.0f;
	for (int i_win = 0; i_win < n_processed_wins; i_win++)
	{
		prof2_var += (prof2_total_sigs_per_win[i_win] - prof2_mean) * (prof2_total_sigs_per_win[i_win] - prof2_mean);
	}
	prof2_var = pow(prof2_var, .5);

	double r = num / (prof1_var * prof2_var);

	return(r);
}

//void generate_valley_length_statistics(char* signal_dir, )
//{
//
//}


void get_2_sample_differential_valleys(char* sample1_valleys_bed_fp, char* sample1_dir, char* sample2_valleys_bed_fp, char* sample2_dir, t_extrema_statistic_defition* extrema_statistic_defn)
{
	fprintf(stderr, "Computing differential valleys between %s (%s) and %s (%s).\n", 
			sample1_valleys_bed_fp, sample1_dir,
			sample2_valleys_bed_fp, sample2_dir);

	// Report the significance method.
	if (extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_BINOMIAL_INTERSECTED_NULLS)
	{
		fprintf(stderr, "Using binomial distribution with multiplication for assessing significance.\n");
	}
	else if (extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_BINOMIAL_UNION_NULLS)
	{
		fprintf(stderr, "Using binomial distribution with merging for assessing significance.\n");
	}
	else if (extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_MULTINOMIAL)
	{
		fprintf(stderr, "Using multinomial distribution for assessing significance.\n");
	}

	vector<t_annot_region*>* sample1_valley_regs = load_BED_with_line_information(sample1_valleys_bed_fp);
	vector<t_annot_region*>* sample2_valley_regs = load_BED_with_line_information(sample2_valleys_bed_fp);

	fprintf(stderr, "%s: %d valleys\n%s: %d valleys\n",
			sample1_valleys_bed_fp, (int)(sample1_valley_regs->size()),
			sample2_valleys_bed_fp, (int)(sample2_valley_regs->size()));

	// Read the header line.
	FILE* f_bed = open_f(sample1_valleys_bed_fp, "r");
	char* header_line = getline(f_bed);
	t_string_tokens* header_line_toks = t_string::tokenize_by_chars(header_line, "\t");
	vector<char*>* header_line_cols = t_string::copy_tokens_2_strs(header_line_toks);
	t_string::clean_tokens(header_line_toks);
	// Get the extrema positions, 
	/*
	fprintf(f_valleys, "#CHROM\t\
	START\t\
	END\t\
	MIN_POSN\t\
	HEIGHT_AT_LEFT\t\
	HEIGHT_AT_RIGHT\t\
	HEIGHT_AT_MIN\t\
	AVG_MMAP\t\
	MAX_MMAP\t\
	LEFT_HILL_SCORE\t\
	RIGHT_HILL_SCORE\t\
	nA\t\
	nC\t\
	nG\t\
	nT\t\
	nCpG\t\
	pVal\t\
	qVal\n");
	*/

	int MIN_POSN_col_i = t_string::get_i_str_ci(header_line_cols, "MIN_POSN");
	fprintf(stderr, "Found minimum position column @ %d\n", MIN_POSN_col_i);
	close_f(f_bed, sample1_valleys_bed_fp);
	
	// Pool the regoins.
	vector<t_annot_region*>* pooled_valleys = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)sample1_valley_regs->size(); i_reg++)
	{
		t_annot_region* cur_valley1 = duplicate_region(sample1_valley_regs->at(i_reg));
		cur_valley1->data = sample1_valley_regs->at(i_reg);
		pooled_valleys->push_back(cur_valley1);
	} // i_reg loop.

	for (int i_reg = 0; i_reg < (int)sample2_valley_regs->size(); i_reg++)
	{
		t_annot_region* cur_valley2 = duplicate_region(sample2_valley_regs->at(i_reg));
		cur_valley2->data = sample2_valley_regs->at(i_reg);
		pooled_valleys->push_back(cur_valley2);
	} // i_reg loop.

	// Re-format, then process.
	t_restr_annot_region_list* restr_pooled_valleys = restructure_annot_regions(pooled_valleys);

	double* log_factorials = buffer_log_factorials(10 * 1000 * 1000);

	int l_frag = 200;

	FILE* f_diff_stats = open_f("2_sample_differential_stats.txt", "w");

	// Write the header.
	fprintf(f_diff_stats, "#CHROM\tSTART\tEND\tDIP\t\
SAMPLE1_LEFT_HEIGHT\tSAMPLE1_RIGHT_HEIGHT\tSAMPLE1_DIP_HEIGHT\t\
SAMPLE2_LEFT_HEIGHT\tSAMPLE2_RIGHT_HEIGHT\tSAMPLE2_DIP_HEIGHT\t\
DIFF12_LEFT_HEIGHT\tDIFF12_RIGHT_HEIGHT\tDIFF12_DIP_HEIGHT\t\
DIFF21_LEFT_HEIGHT\tDIFF21_RIGHT_HEIGHT\tDIFF21_DIP_HEIGHT\t\
SAMPLE1vs2_pval\tSAMPLE2vs1_pval\n");

	for (int i_chr = 0; i_chr < (int)restr_pooled_valleys->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_pooled_valleys = restr_pooled_valleys->regions_per_chrom[i_chr];
		fprintf(stderr, "Processing %d regions on %s\n", (int)cur_chr_pooled_valleys->size(), restr_pooled_valleys->chr_ids->at(i_chr));

		int l_sample1_profile = 0;
		bool sample1_reads_loaded = false;
		double* sample1_profile = load_signal_covg_per_directory_chr_id(sample1_dir, restr_pooled_valleys->chr_ids->at(i_chr), l_frag, l_sample1_profile, sample1_reads_loaded);
		fprintf(stderr, "Loaded %d long profile for %s @ %s\n", l_sample1_profile, sample1_dir, restr_pooled_valleys->chr_ids->at(i_chr));

		int l_sample2_profile = 0;
		bool sample2_reads_loaded = false;
		double* sample2_profile = load_signal_covg_per_directory_chr_id(sample2_dir, restr_pooled_valleys->chr_ids->at(i_chr), l_frag, l_sample2_profile, sample2_reads_loaded);
		fprintf(stderr, "Loaded %d long profile for %s @ %s\n", l_sample2_profile, sample2_dir, restr_pooled_valleys->chr_ids->at(i_chr));

		if (sample1_profile == NULL ||
			sample2_profile == NULL)
		{
			if (sample1_profile != NULL)
			{
				delete[] sample1_profile;
			}

			if (sample2_profile != NULL)
			{
				delete[] sample2_profile;
			}

			fprintf(stderr, "Could not load profile on %s, skipping.\n", restr_pooled_valleys->chr_ids->at(i_chr));
			continue;
		}

		int l_scaling_win = 10000;
		double per_win_2DOF_lls_scaling_factor = 1.0;
		double per_win_1DOF_lls_scaling_factor = 1.0;
		double total_sig_scaling_factor = 1.0;

		get_per_window_scaling_factors_for_profile1_per_profile2(sample1_profile, l_sample1_profile, sample2_profile, l_sample2_profile, l_scaling_win,
			per_win_2DOF_lls_scaling_factor,
			per_win_1DOF_lls_scaling_factor,
			total_sig_scaling_factor);

		// Scale profile: Following ensures we scale the smaller of the profiles, even if they are exchanged, we use the same signal.
		if(total_sig_scaling_factor > 1)
		{
			for (int i = 1; i <= l_sample1_profile; i++)
			{
				sample1_profile[i] *= total_sig_scaling_factor;
			} // i loop.
		}
		else
		{
			for (int i = 1; i <= l_sample2_profile; i++)
			{
				sample2_profile[i] /= total_sig_scaling_factor;
			} // i loop.
		}

		fprintf(stderr, "Scaled 1st profile with %.2f (%lf)\n", per_win_1DOF_lls_scaling_factor, total_sig_scaling_factor);
		 
		// Go over all the pooled valleys and assign differential stats.
		for (int i_pooled_reg = 0; i_pooled_reg < (int)cur_chr_pooled_valleys->size(); i_pooled_reg++)
		{
			t_annot_region* valley_reg = (t_annot_region*)(cur_chr_pooled_valleys->at(i_pooled_reg)->data);

			if (valley_reg == NULL)
			{
				fprintf(stderr, "Sanity check failed, both regions are non-existing for %s:%d-%d\n", 
						cur_chr_pooled_valleys->at(i_pooled_reg)->chrom,
						cur_chr_pooled_valleys->at(i_pooled_reg)->start,
						cur_chr_pooled_valleys->at(i_pooled_reg)->end);

				exit(0);
			}

			// If this region is beyond end of either of the profiles, do not process.
			if (valley_reg->end >= l_sample1_profile ||
				valley_reg->end >= l_sample2_profile)
			{
				continue;
			}

			// Compute the p-value for this valley.
			int l_vic = extrema_statistic_defn->p_val_estimate_extrema_vic_window_length;
			double scaling_factor = extrema_statistic_defn->p_val_estimate_signal_scaling_factor;
			t_valley_significance_info* sig_info = new t_valley_significance_info();
			sig_info->log_q_val = 0;

			char* valley_info_line = (char*)(valley_reg->data);
			t_string_tokens* valley_info_line_toks = t_string::tokenize_by_chars(valley_info_line, "\t");

			int valley_min_posn = atoi(valley_info_line_toks->at(MIN_POSN_col_i)->str());

			t_string::clean_tokens(valley_info_line_toks);

			// Divide the valley into 3 regions and compute enrichment over min position for each of these.
			double valley1_min_vic_sig = 0;
			double valley1_left_max_vic_sig = 0;
			double valley1_right_max_vic_sig = 0;

			// Get valley2 signals.
			double valley2_min_vic_sig = 0;
			double valley2_left_max_vic_sig = 0;
			double valley2_right_max_vic_sig = 0;

			for (int pos = valley_reg->start - l_vic; pos <= valley_reg->start + l_vic; pos++)
			{
				valley1_left_max_vic_sig += sample1_profile[pos];
				valley2_left_max_vic_sig += sample2_profile[pos];
			} // pos loop.

			for (int pos = valley_reg->end - l_vic; pos <= valley_reg->end + l_vic; pos++)
			{
				valley1_right_max_vic_sig += sample1_profile[pos];
				valley2_right_max_vic_sig += sample2_profile[pos];
			} // pos loop.

			for (int pos = valley_min_posn - l_vic; pos <= valley_min_posn + l_vic; pos++)
			{
				valley1_min_vic_sig += sample1_profile[pos];
				valley2_min_vic_sig += sample2_profile[pos];
			} // pos loop.

			// Normalize signals.
			double norm_valley1_min_vic_sig = valley1_min_vic_sig * scaling_factor / (2 * l_vic);
			double norm_valley1_left_max_vic_sig = valley1_left_max_vic_sig * scaling_factor / (2 * l_vic);
			double norm_valley1_right_max_vic_sig = valley1_right_max_vic_sig * scaling_factor / (2 * l_vic);

			double norm_valley2_min_vic_sig = valley2_min_vic_sig * scaling_factor / (2 * l_vic);
			double norm_valley2_left_max_vic_sig = valley2_left_max_vic_sig * scaling_factor / (2 * l_vic);
			double norm_valley2_right_max_vic_sig = valley2_right_max_vic_sig * scaling_factor / (2 * l_vic);

			// Get normalized difference signals.
			double diff21_min_height = MAX(1, norm_valley2_min_vic_sig - norm_valley1_min_vic_sig);
			double diff21_left_max_height = MAX(1, norm_valley2_left_max_vic_sig - norm_valley1_left_max_vic_sig);
			double diff21_right_max_height = MAX(1, norm_valley2_right_max_vic_sig - norm_valley1_right_max_vic_sig);			

			double diff12_min_height = MAX(1, norm_valley1_min_vic_sig - norm_valley2_min_vic_sig);
			double diff12_left_max_height = MAX(1, norm_valley1_left_max_vic_sig - norm_valley2_left_max_vic_sig);
			double diff12_right_max_height = MAX(1, norm_valley1_right_max_vic_sig - norm_valley2_right_max_vic_sig);

			double diff21_valley_pval = xlog(1.0);
			double diff12_valley_pval = xlog(1.0);

			if (extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_BINOMIAL_INTERSECTED_NULLS ||
				extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_BINOMIAL_UNION_NULLS)
			{
				diff21_valley_pval = get_binomial_pval_per_extrema_vals(extrema_statistic_defn->p_val_type, diff21_min_height, diff21_left_max_height, diff21_right_max_height, log_factorials);
				diff12_valley_pval = get_binomial_pval_per_extrema_vals(extrema_statistic_defn->p_val_type, diff12_min_height, diff12_left_max_height, diff12_right_max_height, log_factorials);
			}
			else if (extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_MULTINOMIAL)
			{
				diff21_valley_pval = get_multinomial_pval_per_extrema_vals(diff21_min_height, diff21_left_max_height, diff21_right_max_height, log_factorials);
				diff12_valley_pval = get_multinomial_pval_per_extrema_vals(diff12_min_height, diff12_left_max_height, diff12_right_max_height, log_factorials);
			}

			if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
			{
				fprintf(stderr, "%s:%d-%d-%d: %lf, %lf, %lf: %.5f; %lf, %lf, %lf: %.5f\n",
					cur_chr_pooled_valleys->at(i_pooled_reg)->chrom, valley_reg->start, valley_reg->end, valley_min_posn,
					diff21_min_height, diff21_left_max_height, diff21_right_max_height, diff21_valley_pval,
					diff12_min_height, diff12_left_max_height, diff12_right_max_height, diff12_valley_pval);
			}

			// Save stats.
			double vall1_left_max_height = sample1_profile[valley_reg->start];
			double vall1_right_max_height = sample1_profile[valley_reg->end];
			double vall1_min_height = sample1_profile[valley_min_posn];

			double vall2_left_max_height = sample2_profile[valley_reg->start];
			double vall2_right_max_height = sample2_profile[valley_reg->end];
			double vall2_min_height = sample2_profile[valley_min_posn];

			fprintf(f_diff_stats, "%s\t%d\t%d\t%d\t\
%.1f\t%.1f\t%.1f\t\
%.1f\t%.1f\t%.1f\t\
%.1f\t%.1f\t%.1f\t\
%.1f\t%.1f\t%.1f\t\
%.4f\t%.4f\n",
					cur_chr_pooled_valleys->at(i_pooled_reg)->chrom, cur_chr_pooled_valleys->at(i_pooled_reg)->start, cur_chr_pooled_valleys->at(i_pooled_reg)->end, valley_min_posn,
					vall1_left_max_height, vall1_right_max_height, vall1_min_height,
					vall2_left_max_height, vall2_right_max_height, vall2_min_height,
					diff12_left_max_height, diff12_right_max_height, diff12_min_height,
					diff21_left_max_height, diff21_right_max_height, diff21_min_height,					
					diff12_valley_pval, diff21_valley_pval);
		} // i_pooled_reg loop.

		delete[] sample1_profile;
		delete[] sample2_profile;
	} // i_chr loop.

	fclose(f_diff_stats);
}

bool sort_regions_per_increasing_p_value(t_annot_region* reg1, t_annot_region* reg2)
{
	t_valley_significance_info* reg1_sig_info = (reg1->significance_info);
	t_valley_significance_info* reg2_sig_info = (reg2->significance_info);

	return(reg1_sig_info->log_p_val < reg2_sig_info->log_p_val);
}

bool sort_regions_per_increasing_q_value(t_annot_region* reg1, t_annot_region* reg2)
{
	t_valley_significance_info* reg1_sig_info = (reg1->significance_info);
	t_valley_significance_info* reg2_sig_info = (reg2->significance_info);

	return(reg1_sig_info->log_q_val < reg2_sig_info->log_q_val);
}

void get_benjamini_hochberg_corrected_p_values_per_valleys(vector<t_annot_region*>* regions)
{
	sort(regions->begin(), regions->end(), sort_regions_per_increasing_p_value);

	// Do BH corrections.
	for (int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		double cur_region_rank = (double)i_reg + 1;

		t_valley_significance_info* cur_reg_sig_info = regions->at(i_reg)->significance_info;
		cur_reg_sig_info->log_q_val = (cur_reg_sig_info->log_p_val + log((double)regions->size() / (double)cur_region_rank));
	} // i_reg loop.

	  // Sort with respect to q-values.
	sort(regions->begin(), regions->end(), sort_regions_per_increasing_q_value);
}

double* load_signal_covg_per_signal_file(const char* cur_dat_fp,
	int l_fragment,
	int& l_loaded_covg,
	bool& reads_loaded)
{
	double* covg_signal = NULL;

	reads_loaded = false;

	// Search for mapped reads.
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

		return(covg_signal);
	}

	// Search for BGR.
	if (check_file(cur_dat_fp) &&
		t_string::ends_with(cur_dat_fp, ".bgr.gz"))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	if (check_file(cur_dat_fp) &&
		t_string::ends_with(cur_dat_fp, ".bgr"))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	l_loaded_covg = -1;
	return(NULL);
}

double* load_signal_covg_per_directory_chr_id(const char* dat_dir,
	char* chr_id,
	int l_fragment,
	int& l_loaded_covg,
	bool& reads_loaded)
{
	if (dat_dir == NULL)
	{
		return(NULL);
	}

	double* covg_signal = NULL;

	reads_loaded = false;

	// Search for mapped reads.
	char cur_dat_fp[1000];
	sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;

		double* buffer = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			buffer, NULL, NULL,
			l_buffer, l_loaded_covg);

		double* covg_signal = new double[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = buffer[i];
		} // i loop.
		delete[] buffer;

		return(covg_signal);
	}

	sprintf(cur_dat_fp, "%s/%s_mapped_reads.txt.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		reads_loaded = true;
		int l_buffer = 300 * 1000 * 1000;

		double* buffer = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			buffer, NULL, NULL,
			l_buffer, l_loaded_covg);

		double* covg_signal = new double[l_loaded_covg + 2];
		for (int i = 0; i <= l_loaded_covg; i++)
		{
			covg_signal[i] = buffer[i];
		} // i loop.
		delete[] buffer;

		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/signal_%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/%s.bgr.gz", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	// Search for BGR.
	sprintf(cur_dat_fp, "%s/%s.bgr", dat_dir, chr_id);
	if (check_file(cur_dat_fp))
	{
		covg_signal = load_per_nucleotide_BGR_track(cur_dat_fp, l_loaded_covg);
		return(covg_signal);
	}

	l_loaded_covg = -1;
	return(NULL);
}

void binarize_fasta_file(char* fasta_fp, char* bin_dir)
{
	printf("Saving %s.\n", fasta_fp);
	FILE* f_fasta = open_f(fasta_fp, "r");

	char* cur_line = NULL;
	char* cur_entry_buffer = new char[250 * 1000 * 1000];
	int cur_entry_i = 1; // This ensures that we are 1 based.
	char cur_entry_id[1000];

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", bin_dir);
	FILE* f_chr_ids = open_f(chr_ids_fp, "a");

	while (1)
	{
		// Process the current buffer.
		cur_line = getline(f_fasta);

		if (cur_line == NULL)
		{
			// File ended, dump the last entry if there are values in it.
			if (cur_entry_i > 1)
			{
				char cur_entry_bin_fp[1000];
				normalize_chr_id(cur_entry_id);
				sprintf(cur_entry_bin_fp, "%s/%s.bin", bin_dir, cur_entry_id);
				fprintf(f_chr_ids, "%s\n", cur_entry_id);
				FILE* f_bin = open_f(cur_entry_bin_fp, "wb");

				// Open the new binary file.
				fprintf(stderr, "Dumping %s (%d)\n", cur_entry_id, cur_entry_i);

				// Dump the sequence buffer.
				fwrite(&cur_entry_i, sizeof(int), 1, f_bin);
				fwrite(cur_entry_buffer, sizeof(char), cur_entry_i, f_bin);

				// Close current file.
				fclose(f_bin);
			}
			break;
		}
		else if (cur_line[0] == '>')
		{
			if (cur_entry_i > 1)
			{
				char cur_entry_bin_fp[1000];
				normalize_chr_id(cur_entry_id);
				sprintf(cur_entry_bin_fp, "%s/%s.bin", bin_dir, cur_entry_id);
				fprintf(f_chr_ids, "%s\n", cur_entry_id);
				FILE* f_bin = open_f(cur_entry_bin_fp, "wb");

				// Open the new binary file.
				fprintf(stderr, "Dumping %s (%d)\n", cur_entry_id, cur_entry_i);

				// Dump the sequence buffer.
				fwrite(&cur_entry_i, sizeof(int), 1, f_bin);
				fwrite(cur_entry_buffer, sizeof(char), cur_entry_i, f_bin);

				// Close current file.
				fclose(f_bin);
			}

			// Update the id, reset the counter.
			strcpy(cur_entry_id, &(cur_line[1]));
			cur_entry_i = 1;
		}
		else
		{
			// Concatenate the current sequence line.
			int l_cur_line = t_string::string_length(cur_line);
			for (int i = 0; i < l_cur_line; i++)
			{
				cur_entry_buffer[cur_entry_i] = cur_line[i];
				cur_entry_i++;
			} // i loop.
		}

		delete[] cur_line;
	} // file reading loop.

	fclose(f_fasta);
	fclose(f_chr_ids);

	delete[] cur_entry_buffer;
}

char* load_binary_sequence_file(char* bin_seq_fp, int& l_seq)
{
	FILE* f_bin_seq = open_f(bin_seq_fp, "rb");
	int l_bin_seq = 0;

	// Read the length.
	fread(&l_bin_seq, sizeof(int), 1, f_bin_seq);

	l_seq = l_bin_seq;

	// Read the sequence.
	char* bin_seq = new char[l_bin_seq + 2];
	fread(bin_seq, sizeof(char), l_bin_seq + 1, f_bin_seq);
	fclose(f_bin_seq);

	return(bin_seq);
}

double* buffer_log_factorials(int n)
{
	double* factorials = new double[n + 2];
	factorials[0] = xlog(1.0);
	for (int i = 1; i <= n; i++)
	{
		factorials[i] = xlog_mul(factorials[i - 1], xlog(i));
	} // i loop.

	return(factorials);
}

double get_multinomial_pval_per_extrema_vals(double norm_min_vic_sig, double norm_left_max_vic_sig, double norm_right_max_vic_sig, double* log_factorials)
{
	norm_min_vic_sig = round(norm_min_vic_sig);
	norm_left_max_vic_sig = round(norm_left_max_vic_sig);
	norm_right_max_vic_sig = round(norm_right_max_vic_sig);

	norm_min_vic_sig = MAX(norm_min_vic_sig, 1);
	norm_left_max_vic_sig = MAX(norm_left_max_vic_sig, 1);
	norm_right_max_vic_sig = MAX(norm_right_max_vic_sig, 1);

	// Add the probability of enrichment on the left side.
	int grand_total_int = (int)(norm_min_vic_sig + norm_left_max_vic_sig + norm_right_max_vic_sig);
	double log_multinom_p_val = xlog(0);
	double log_flip = -1 * log(3);

	// Go over all the signal levels at the minimum and distribute.
	for (int cur_min_sig = norm_min_vic_sig; cur_min_sig >= 0; cur_min_sig--)
	{
		int n_reads_2_distribute = norm_min_vic_sig - cur_min_sig;

		// Update the current read configuration.
		for (int n_new_left_reads = 0; n_new_left_reads <= n_reads_2_distribute; n_new_left_reads++)
		{
			int n_new_right_reads = n_reads_2_distribute - n_new_left_reads;

			int n_left_reads = n_new_left_reads + norm_left_max_vic_sig;
			int n_right_reads = n_new_right_reads + norm_right_max_vic_sig;

			if (cur_min_sig + n_left_reads + n_right_reads != grand_total_int)
			{
				fprintf(stderr, "Sanity check failed in multinomial test: %d, %d, %d reads vs %d reads; (%.1f, %.1f, %.1f)\n",
					cur_min_sig, n_left_reads, n_right_reads, grand_total_int,
					norm_min_vic_sig, norm_left_max_vic_sig, norm_right_max_vic_sig);

				exit(0);
			}

			// Compute the probability.
			double log_flip_prob = log_flip * grand_total_int;
			double log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[n_left_reads], xlog_mul(log_factorials[n_right_reads], log_factorials[cur_min_sig])));

			log_multinom_p_val = xlog_sum(log_multinom_p_val, xlog_mul(log_cur_perm, log_flip_prob));
		} // left_i loop.
	} // cur_min_sig loop.

	return(log_multinom_p_val);
}

double get_valley_significance_per_multinomial_test(double* signal_profile, int l_profile, int left_max_posn, int right_max_posn, int min_posn, int l_normalizer, int scaling_factor, double* _log_factorials)
{
	// Divide the valley into 3 regions and compute enrichment over min position for each of these.
	double min_vic_sig = 0;
	double left_max_vic_sig = 0;
	double right_max_vic_sig = 0;

	for (int pos = left_max_posn - l_normalizer; pos <= left_max_posn + l_normalizer; pos++)
	{
		left_max_vic_sig += signal_profile[pos];
	} // pos loop.

	for (int pos = right_max_posn - l_normalizer; pos <= right_max_posn + l_normalizer; pos++)
	{
		right_max_vic_sig += signal_profile[pos];
	} // pos loop.

	for (int pos = min_posn - l_normalizer; pos <= min_posn + l_normalizer; pos++)
	{
		min_vic_sig += signal_profile[pos];
	} // pos loop.

	double norm_min_vic_sig = min_vic_sig * scaling_factor / l_normalizer;
	double norm_left_max_vic_sig = left_max_vic_sig * scaling_factor / l_normalizer;
	double norm_right_max_vic_sig = right_max_vic_sig * scaling_factor / l_normalizer;	

	// Setup the log factorials.
	double* log_factorials = NULL;
	if (_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(100 * 1000 + 3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	double log_multinom_p_val = get_multinomial_pval_per_extrema_vals(norm_min_vic_sig, norm_left_max_vic_sig, norm_right_max_vic_sig, log_factorials);

	// Free memory if it is allocated in the function.
	if (_log_factorials == NULL)
	{
		delete[] log_factorials;
	}

	return(log_multinom_p_val);

	//norm_min_vic_sig = round(norm_min_vic_sig);
	//norm_left_max_vic_sig = round(norm_left_max_vic_sig);
	//norm_right_max_vic_sig = round(norm_right_max_vic_sig);

	//norm_min_vic_sig = MAX(norm_min_vic_sig, 1);
	//norm_left_max_vic_sig = MAX(norm_left_max_vic_sig, 1);
	//norm_right_max_vic_sig = MAX(norm_right_max_vic_sig, 1);

	//// Add the probability of enrichment on the left side.
	//int grand_total_int = (int)(norm_min_vic_sig + norm_left_max_vic_sig + round(norm_right_max_vic_sig));
	//double log_multinom_p_val = xlog(0);
	//double log_flip = -1 * log(3);

	//// Go over all the signal levels at the minimum and distribute.
	//for (int cur_min_sig = norm_min_vic_sig; cur_min_sig >= 0; cur_min_sig--)
	//{
	//	int n_reads_2_distribute = norm_min_vic_sig - cur_min_sig;

	//	// Update the current read configuration.
	//	for (int n_new_left_reads = 0; n_new_left_reads <= n_reads_2_distribute; n_new_left_reads++)
	//	{
	//		int n_new_right_reads = n_reads_2_distribute - n_new_left_reads;

	//		int n_left_reads = n_new_left_reads + norm_left_max_vic_sig;
	//		int n_right_reads = n_new_right_reads + norm_right_max_vic_sig;

	//		if (cur_min_sig + n_left_reads + n_right_reads != grand_total_int)
	//		{
	//			fprintf(stderr, "Sanity check failed in multinomial test: %d, %d, %d reads vs %d reads; (%.1f, %.1f, %.1f)\n", 
	//				cur_min_sig, n_left_reads, n_right_reads, grand_total_int,
	//				norm_min_vic_sig, norm_left_max_vic_sig, norm_right_max_vic_sig);

	//			exit(0);
	//		}

	//		// Compute the probability.
	//		double log_flip_prob = log_flip * grand_total_int;
	//		double log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[n_left_reads], xlog_mul(log_factorials[n_right_reads], log_factorials[cur_min_sig])));

	//		log_multinom_p_val = xlog_sum(log_multinom_p_val, xlog_mul(log_cur_perm, log_flip_prob));
	//	} // left_i loop.
	//} // cur_min_sig loop.

	//if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Valley: %d-%d-%d: %lf, %lf: %lf (flip_prob: %.3f)\n", left_max_posn, min_posn, right_max_posn, norm_min_vic_sig, norm_left_max_vic_sig, norm_right_max_vic_sig, exp(log_flip));
	//}

	//// Free memory if it is allocated in the function.
	//if (_log_factorials == NULL)
	//{
	//	delete[] log_factorials;
	//}

	//return(log_multinom_p_val);
}

double get_binomial_pval_per_extrema_vals(int p_val_type, double norm_min_vic_sig, double norm_left_max_vic_sig, double norm_right_max_vic_sig, double* log_factorials)
{
	norm_min_vic_sig = round(norm_min_vic_sig);
	norm_left_max_vic_sig = round(norm_left_max_vic_sig);
	norm_right_max_vic_sig = round(norm_right_max_vic_sig);

	norm_min_vic_sig = MAX(norm_min_vic_sig, 1);
	norm_left_max_vic_sig = MAX(norm_left_max_vic_sig, 1);
	norm_right_max_vic_sig = MAX(norm_right_max_vic_sig, 1);

	// Add the probability of enrichment on the left side.
	int grand_total_int = (int)(norm_min_vic_sig + norm_left_max_vic_sig);
	double left_region_log_p_val = xlog(0);
	double flip_prob = 0.5;
	double log_flip = log(flip_prob);
	for (int i = 0; i <= norm_min_vic_sig; i++)
	{
		double log_cur_half_pow = log_flip * grand_total_int;
		double log_cur_perm = 0.0; // = xlog(1.0).

								   // Compute the current permutation.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int - i]));
		left_region_log_p_val = xlog_sum(left_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	  // Add the probability of enrichment on the right side.
	double right_region_log_p_val = xlog(0);
	grand_total_int = (int)(norm_min_vic_sig + norm_right_max_vic_sig);
	for (int i = 0; i <= norm_min_vic_sig; i++)
	{
		double log_cur_half_pow = log_flip * grand_total_int;
		double log_cur_perm = 0.0; // = xlog(1.0).

								   // Compute the current permutation.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int - i]));
		right_region_log_p_val = xlog_sum(right_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
	{
		//fprintf(stderr, "Valley: %d-%d-%d: %lf, %lf: %lf (flip_prob: %.3f)\n", left_max_posn, min_posn, right_max_posn, norm_min_vic_sig, norm_left_max_vic_sig, norm_right_max_vic_sig, exp(log_flip));
		fprintf(stderr, "Valley: %lf, %lf: %lf (flip_prob: %.3f)\n", norm_min_vic_sig, norm_left_max_vic_sig, norm_right_max_vic_sig, exp(log_flip));
	}

	//// Free memory if it is allocated in the function.
	//if (_log_factorials == NULL)
	//{
	//	delete[] log_factorials;
	//}

	//double merged_p_val = xlog(0);

	if (p_val_type == VALLEY_SIGNIFICANCE_BINOMIAL_INTERSECTED_NULLS)
	{
		return(right_region_log_p_val + left_region_log_p_val);
	}
	else if (p_val_type == VALLEY_SIGNIFICANCE_BINOMIAL_UNION_NULLS)
	{
		double left_right = left_region_log_p_val + right_region_log_p_val;

		// 1st way.
		//double one_minus_left_region_log_p_val = xlog_sub(xlog(1.0), left_region_log_p_val);
		//double one_minus_right_region_log_p_val = xlog_sub(xlog(1.0), right_region_log_p_val);

		//double right_min_left = one_minus_left_region_log_p_val + right_region_log_p_val;
		//double left_min_right = one_minus_right_region_log_p_val + left_region_log_p_val;
		//return(xlog_sum(left_right, xlog_sum(left_min_right, right_min_left)));

		// 2nd way.		
		double merged_p_val = xlog_sum(left_region_log_p_val, right_region_log_p_val);
		merged_p_val = xlog_sub(merged_p_val, left_right);

		return(merged_p_val);
	}
	else
	{
		fprintf(stderr, "No p-value type selected @ %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
}

// Compute binomial p-value of the signal imbalance on the valley.
double get_valley_significance_per_binomial_test(int p_val_type, double* signal_profile, int l_profile, int left_max_posn, int right_max_posn, int min_posn, int l_normalizer, int scaling_factor, double* _log_factorials)
{
	// Divide the valley into 3 regions and compute enrichment over min position for each of these.
	double min_vic_sig = 0;
	double left_max_vic_sig = 0;
	double right_max_vic_sig = 0;

	for (int pos = left_max_posn - l_normalizer; pos <= left_max_posn + l_normalizer; pos++)
	{
		left_max_vic_sig += signal_profile[pos];
	} // pos loop.

	for (int pos = right_max_posn - l_normalizer; pos <= right_max_posn + l_normalizer; pos++)
	{
		right_max_vic_sig += signal_profile[pos];
	} // pos loop.

	for (int pos = min_posn - l_normalizer; pos <= min_posn + l_normalizer; pos++)
	{
		min_vic_sig += signal_profile[pos];
	} // pos loop.

	double norm_min_vic_sig = min_vic_sig* scaling_factor / l_normalizer;
	double norm_left_max_vic_sig = left_max_vic_sig * scaling_factor / l_normalizer;
	double norm_right_max_vic_sig = right_max_vic_sig * scaling_factor / l_normalizer;

	// Setup the log factorials.
	double* log_factorials = NULL;
	if (_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(100*1000 + 3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	double merged_p_val = get_binomial_pval_per_extrema_vals(p_val_type, norm_min_vic_sig, norm_left_max_vic_sig, norm_right_max_vic_sig, log_factorials);

	// Free memory if it is allocated in the function.
	if (_log_factorials == NULL)
	{
		delete[] log_factorials;
	}

	return(merged_p_val);
}

double get_binomial_p_val_per_enriched_background_values(int enriched_value, int background_value, double* _log_factorials)
{
	int grand_total_int = (int)(enriched_value + background_value);

	// Setup the log factorials.
	double* log_factorials = NULL;
	if (_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(grand_total_int + 3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	// Compute the p-values.
	double log_flip = log(.5);
	double cur_region_log_p_val = xlog(0.0);
	for (int i = 0; i <= background_value; i++)
	{
		double log_cur_half_pow = log_flip * grand_total_int;
		double log_cur_perm = 0.0; // = xlog(1.0).

								   // Compute the current permutation, using shortcuts, bypassing the function calls from xlog_math library.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int - i]));
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	  // Free memory if it is allocated in the function.
	if (_log_factorials == NULL)
	{
		delete[] log_factorials;
	}
	else
	{
	}

	return(cur_region_log_p_val);
}

void refine_valleys_per_pval_min(double* signal_profile, int l_profile, vector<t_annot_region*>* valley_regs, int l_vic_expand)
{
	fprintf(stderr, "Refining %d valleys.\n", (int)(valley_regs->size()));
	double* log_factorials = buffer_log_factorials(1000 * 1000);
	for (int i_reg = 0; i_reg < (int)(valley_regs->size()); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Refining %d. valley.         \r", i_reg);
		}

		void** valley_info = (void**)(valley_regs->at(i_reg)->data);

		t_extrema_node** valley_points = (t_extrema_node**)(valley_info[0]);

		int exp_start = MAX(1, valley_regs->at(i_reg)->start - l_vic_expand);
		int exp_end = MIN(l_profile, valley_regs->at(i_reg)->end + l_vic_expand);

		int opt_start = valley_regs->at(i_reg)->start;
		int opt_end = valley_regs->at(i_reg)->end;
		double opt_enrichment = -1234;

		double max_signal = -1 * 1000 * 1000;
		double min_signal = 1000 * 1000;
		for (int i = exp_start; i < exp_end; i++)
		{
			max_signal = MAX(max_signal, signal_profile[i]);
			min_signal = MIN(min_signal, signal_profile[i]);
		} // i loop.

		int min_pos = valley_points[2]->extrema_posn;	

		// Thresh loop.
		for(int thresh = (int)min_signal; thresh < (int)max_signal; thresh++)
		{
			int cur_start = min_pos - 1;
			while (cur_start > exp_start &&
				signal_profile[cur_start] < thresh)
			{
				cur_start--;
			}

			int cur_end = min_pos + 1;
			while (cur_end < exp_end &&
				signal_profile[cur_end] < thresh)
			{
				cur_end++;
			}

			double cur_tot_signal = 0;
			for (int i = cur_start; i < cur_end; i++)
			{
				cur_tot_signal += max_signal - signal_profile[i];
			} // i loop.

			double cur_enrichment = get_binomial_p_val_per_enriched_background_values((int)(cur_tot_signal / 200), 5, log_factorials);

			if (opt_enrichment == -1234 ||
				cur_enrichment < opt_enrichment)
			{
				opt_start = cur_start;
				opt_end = cur_end;
				opt_enrichment = cur_enrichment;
			}
		} // thresh loop.

		if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
		{
			fprintf(stderr, "%d-%d -> %d-%d\n", valley_regs->at(i_reg)->start, valley_regs->at(i_reg)->end, opt_start, opt_end);
		}

		valley_regs->at(i_reg)->start = opt_start;
		valley_regs->at(i_reg)->end = opt_end;
	} // i_reg loop.

	delete[] log_factorials;
}

void refine_valleys_per_height_trim(double* signal_profile, int l_profile, vector<t_annot_region*>* valley_regs, double height_frac)
{
	fprintf(stderr, "Refining %d valleys.\n", (int)(valley_regs->size()));

	for (int i_reg = 0; i_reg < (int)(valley_regs->size()); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Refining %d. valley.         \r", i_reg);
		}

		void** valley_info = (void**)(valley_regs->at(i_reg)->data);

		t_extrema_node** valley_points = (t_extrema_node**)(valley_info[0]);

		int exp_start = MAX(1, valley_regs->at(i_reg)->start);
		int exp_end = MIN(l_profile, valley_regs->at(i_reg)->end);

		int opt_start = valley_regs->at(i_reg)->start;
		int opt_end = valley_regs->at(i_reg)->end;

		double max_signal = -1 * 1000 * 1000;
		double min_signal = 1000 * 1000;
		for (int i = exp_start; i < exp_end; i++)
		{
			max_signal = MAX(max_signal, signal_profile[i]);
			min_signal = MIN(min_signal, signal_profile[i]);
		} // i loop.

		int min_pos = valley_points[2]->extrema_posn;

		int cur_start = min_pos - 1;
		while (cur_start > exp_start &&
			signal_profile[cur_start] < signal_profile[exp_start] * height_frac)
		{
			cur_start--;
		}

		int cur_end = min_pos + 1;
		while (cur_end < exp_end &&
			signal_profile[cur_end] < signal_profile[exp_end] * height_frac)
		{
			cur_end++;
		}

		if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
		{
			fprintf(stderr, "%d-%d -> %d-%d\n", valley_regs->at(i_reg)->start, valley_regs->at(i_reg)->end, opt_start, opt_end);
		}

		valley_regs->at(i_reg)->start = cur_start;
		valley_regs->at(i_reg)->end = cur_end;
	} // i_reg loop.
}

double* bin_profile(double* profile, int l_profile, int l_bin, int& n_bins)
{
	n_bins = 0;
	int max_n_bins = l_profile / l_bin + 100;
	int cur_bin_start = 1;
	double* binned_profile = new double[max_n_bins];
	memset(binned_profile, 0, sizeof(double) * max_n_bins);
	while(cur_bin_start < l_profile)
	{
		double cur_bin_sig = 0;
		for (int i = cur_bin_start; MAX(l_profile, cur_bin_start + l_bin); i++)
		{
			cur_bin_sig += profile[i];
		} // i loop.

		// Set the binned profile value.
		binned_profile[n_bins] = cur_bin_sig / l_bin;
		n_bins++;

		// Update bin start.
		cur_bin_start += l_bin;
	} // binning loop.

	return(binned_profile);
}

vector<t_annot_region*>* get_significant_extrema_per_signal_profile(const char* op_dir,
																	const char* chr_id,
																	double* signal_profile, int l_profile,
																	unsigned char* multimapp_signal_profile, int l_multimapp_profile,
																	char* chrom_seq,
																	t_extrema_statistic_defition* extrema_statistic_defn)
{
	// Remove the negative values.
	for (int i = 0; i <= l_profile; i++)
	{
		signal_profile[i] = MAX(0, signal_profile[i]);
	} // i loop.

	if (extrema_statistic_defn->hill_score_type == DIST_BASED_HILL_SCORE)
	{
		fprintf(stderr, "Assigning distance based hill scores.\n");
	}
	else if (extrema_statistic_defn->hill_score_type == HEIGHT_BASED_HILL_SCORE)
	{
		fprintf(stderr, "Assigning height based hill scores.\n");
	}

	// This is necessary to make computes faster for very long valleys.
	//// Bin the signal for very long regions?
	//int l_bin = 1;
	//if (extrema_statistic_defn->max_summit2trough_dist_in_bp > 10000)
	//{
	//	// Create 10 bp bins.
	//	l_bin = 10;
	//}
	//else if(extrema_statistic_defn->max_summit2trough_dist_in_bp > 100000)
	//{
	//	// Create 100 base pair bins
	//	l_bin = 100;
	//}
	//else if (extrema_statistic_defn->max_summit2trough_dist_in_bp > 1000000)
	//{
	//	l_bin = 1000;
	//}

	//// Normalize signal?
	//double* binned_profile = bin_profile(signal_profile, l_profile, l_bin);
	//double* orig_signal_profile = signal_profile;
	//signal_profile = binned_profile;

	fprintf(stderr, "Computing extremas.\n");
	vector<t_extrema_node*>* maxima_nodes = new vector<t_extrema_node*>();
	vector<t_extrema_node*>* minima_nodes = new vector<t_extrema_node*>();
	int* derivative_map = new int[l_profile + 2];
	memset(derivative_map, 0, sizeof(int) * (l_profile+1));
	get_extrema_per_plateaus(signal_profile, l_profile, maxima_nodes, minima_nodes, derivative_map, 0);

	fprintf(stderr, "Identified %d minima and %d maxima.\n", (int)minima_nodes->size(), (int)maxima_nodes->size());

	vector<t_extrema_node*>* all_extrema_regs = new vector<t_extrema_node*>();
	all_extrema_regs->insert(all_extrema_regs->end(), minima_nodes->begin(), minima_nodes->end());
	all_extrema_regs->insert(all_extrema_regs->end(), maxima_nodes->begin(), maxima_nodes->end());

	// Sort the extrema with respect to position.
	sort(all_extrema_regs->begin(), all_extrema_regs->end(), sort_extremas_per_posn);

	fprintf(stderr, "Pooled %d extrema, processing the valleys with trough2summit distance in %d-%d bps.\n", 
			(int)all_extrema_regs->size(), 
			extrema_statistic_defn->min_summit2trough_dist_in_bp,
			extrema_statistic_defn->max_summit2trough_dist_in_bp);

	double* log_factorials = buffer_log_factorials(1000 * 1000);

	// Dump the extremas: Dump the minima and maxima depending on their "peakiness": Find maxima around minima that satisfies the peakiness.
	vector<t_annot_region*>* significant_valleys = new vector<t_annot_region*>();
	for (int i_ext = 0; i_ext < (int)all_extrema_regs->size(); i_ext++)
	{
		if (i_ext % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. extrema        \r", i_ext);
		}

		// Process minima: For every minima, check if it is a significant minima.
		if (all_extrema_regs->at(i_ext)->extrema_type == EXTREMA_MIN)
		{
			// Check the signal at this minima here to make it faster to go through the data.
			if (all_extrema_regs->at(i_ext)->height_at_extrema > extrema_statistic_defn->max_signal_at_trough)
			{
				continue;
			}

			vector<t_extrema_node*>** left_right_maxima = new vector<t_extrema_node*>*[2];
			left_right_maxima[0] = new vector<t_extrema_node*>();
			left_right_maxima[1] = new vector<t_extrema_node*>();

			// Find the maxima around the minima.
			for (int dir = 0; dir < 2; dir++)
			{
				int index_update = (dir == 0) ? (-1) : (1);

				for (int j_ext = i_ext;
						j_ext > 0 && 
						j_ext < (int)all_extrema_regs->size() &&
						fabs((double)(all_extrema_regs->at(i_ext)->extrema_posn - all_extrema_regs->at(j_ext)->extrema_posn)) < extrema_statistic_defn->max_summit2trough_dist_in_bp;
					j_ext += index_update)
				{
					if (all_extrema_regs->at(j_ext)->extrema_type == EXTREMA_MAX)
					{
						if (all_extrema_regs->at(j_ext)->height_at_extrema > extrema_statistic_defn->min_signal_at_summit &&
							all_extrema_regs->at(j_ext)->height_at_extrema / all_extrema_regs->at(i_ext)->height_at_extrema > extrema_statistic_defn->min_summit2trough_ratio_per_trough &&
							fabs((double)(all_extrema_regs->at(i_ext)->extrema_posn - all_extrema_regs->at(j_ext)->extrema_posn)) > extrema_statistic_defn->min_summit2trough_dist_in_bp)
						{
							left_right_maxima[dir]->push_back(all_extrema_regs->at(j_ext));
						} // extrema is a summit check.
					} // min2max distance check.
				} // j_ext loop.
			} // dir loop.

			// Find the peakiest valley among the list of summits around this minima.
			for (int left_max_i = 0; left_max_i < (int)(left_right_maxima[0]->size()); left_max_i++)
			{
				for (int right_max_i = 0; right_max_i < (int)(left_right_maxima[1]->size()); right_max_i++)
				{
					t_extrema_node* cur_left_maxima = left_right_maxima[0]->at(left_max_i);
					t_extrema_node* cur_right_maxima = left_right_maxima[1]->at(right_max_i);
					t_extrema_node* cur_minima = all_extrema_regs->at(i_ext);

					// Compute the total derivative map to right and left.
					double* left_right_pos_der_frac = new double[2];

					for (int dir = 0; dir < 2; dir++)
					{
						double cur_dir_total_pos_deriv_signs = 0;
						double cur_dir_total_neg_deriv_signs = 0;

						double cur_dir_total_pos_signal_increase = 0;
						double cur_dir_total_pos_signal_decrease = 0;

						int update = (dir == 0) ? (-1) : (1);

						for (int i = cur_minima->extrema_posn; 
							i > cur_left_maxima->extrema_posn && i < cur_right_maxima->extrema_posn; 
							i += update)
						{
							double cur_sig_diff = signal_profile[i] - signal_profile[i - 1];
							double cur_deriv = derivative_map[i];

							// If moving to left, revert the signs of derivative and signal difference.
							if (dir == 0)
							{
								cur_deriv *= -1;
								cur_sig_diff *= -1;
							}

							if(cur_deriv >= 0)
							{
								cur_dir_total_pos_deriv_signs++;
							}

							if(cur_deriv <= 0)
							{
								cur_dir_total_neg_deriv_signs++;
							}

							// Following updates signal weighted hill score.
							if (cur_sig_diff > 0)
							{
								cur_dir_total_pos_signal_increase += cur_sig_diff;
							}

							if (cur_sig_diff < 0)
							{
								cur_dir_total_pos_signal_decrease += -1*cur_sig_diff;
							}
						} // i loop.

						// Increase/decrease count based hill score.
						if (extrema_statistic_defn->hill_score_type == DIST_BASED_HILL_SCORE)
						{
							left_right_pos_der_frac[dir] = cur_dir_total_pos_deriv_signs / (cur_dir_total_pos_deriv_signs + cur_dir_total_neg_deriv_signs);
						}

						// Increase/decrease signal based hill score.
						if (extrema_statistic_defn->hill_score_type == HEIGHT_BASED_HILL_SCORE)
						{
							left_right_pos_der_frac[dir] = cur_dir_total_pos_signal_increase / (cur_dir_total_pos_signal_increase + cur_dir_total_pos_signal_decrease);
						}
					} // dir loop.

					// Compute the average multi-map signal within the valley.
					double total_multimapp_signal = 0;
					double max_multimapp_signal = -1;

					if (multimapp_signal_profile != NULL)
					{
						for (int i = cur_left_maxima->extrema_posn; i < cur_right_maxima->extrema_posn; i++)
						{
							if (i < l_multimapp_profile)
							{
								double multimapp_val = ((double)(multimapp_signal_profile[i])) / 100.0;

								if (multimapp_val > max_multimapp_signal)
								{
									max_multimapp_signal = multimapp_val;
								}

								total_multimapp_signal += multimapp_val;
							}
						} // i loop.
					}
					else
					{
						max_multimapp_signal = 0;
					}

					// Count the nucleotides.
					int* nuc_counts = new int[10];
					memset(nuc_counts, 0, sizeof(int) * 10);
					int n_CpGs = 0;

					if (chrom_seq != NULL)
					{
						for (int i = cur_left_maxima->extrema_posn; i < cur_right_maxima->extrema_posn; i++)
						{
							nuc_counts[nuc_2_num(chrom_seq[i])]++;

							if (toupper(chrom_seq[i]) == 'C' &&
								toupper(chrom_seq[i + 1]) == 'G')
							{
								n_CpGs++;
							}
						} // i loop.
					}

					nuc_counts[5] = n_CpGs;

					total_multimapp_signal /= fabs((double)(cur_left_maxima->extrema_posn - cur_right_maxima->extrema_posn));
				
					t_annot_region* new_valley = get_empty_region();
					new_valley->chrom = t_string::copy_me_str(chr_id);
					new_valley->start = cur_left_maxima->extrema_posn;
					new_valley->end = cur_right_maxima->extrema_posn;
				
					t_extrema_node** valley_points = new t_extrema_node*[3];
					valley_points[0] = cur_left_maxima;
					valley_points[1] = cur_right_maxima;
					valley_points[2] = all_extrema_regs->at(i_ext);

					// Set the points.
					void** valley_info = new void*[10];
					valley_info[0] = valley_points;

					// Set the multimapp signal levels.
					double* cur_valley_multimapp_signals = new double[3];
					cur_valley_multimapp_signals[0] = total_multimapp_signal;
					cur_valley_multimapp_signals[1] = max_multimapp_signal;
					valley_info[1] = cur_valley_multimapp_signals;

					// Set the positive derivative fractions.
					valley_info[2] = left_right_pos_der_frac;

					// Set the nucleotide counts as information.
					valley_info[3] = nuc_counts;

					// Compute the p-value for this valley.
					int l_vic = extrema_statistic_defn->p_val_estimate_extrema_vic_window_length;
					double scaling_factor = extrema_statistic_defn->p_val_estimate_signal_scaling_factor;
					t_valley_significance_info* sig_info = new t_valley_significance_info();
					sig_info->log_q_val = 0;

					if (extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_BINOMIAL_INTERSECTED_NULLS ||
						extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_BINOMIAL_UNION_NULLS)
					{
						sig_info->log_p_val = get_valley_significance_per_binomial_test(extrema_statistic_defn->p_val_type, signal_profile, l_profile, cur_left_maxima->extrema_posn, cur_right_maxima->extrema_posn, cur_minima->extrema_posn, l_vic, scaling_factor, log_factorials);
					}
					else if (extrema_statistic_defn->p_val_type == VALLEY_SIGNIFICANCE_MULTINOMIAL)
					{
						sig_info->log_p_val = get_valley_significance_per_multinomial_test(signal_profile, l_profile, cur_left_maxima->extrema_posn, cur_right_maxima->extrema_posn, cur_minima->extrema_posn, l_vic, scaling_factor, log_factorials);
					}

					//double left_2_ get_binomial_p_val_per_enriched_background_values(signal_profile[cur_left_maxima->extrema_posn], int background_value, double* _log_factorials)
					
					new_valley->significance_info = sig_info;

					new_valley->data = valley_info;

					significant_valleys->push_back(new_valley);
				} // right_max_i loop.
			} // left_max_i loop.
		} // minima check.
	} // i_ext loop.

	// Refine the valleys.
	if (extrema_statistic_defn->sparse_profile)
	{
		refine_valleys_per_height_trim(signal_profile, l_profile, significant_valleys, 0.99);
	}
	else
	{
		refine_valleys_per_pval_min(signal_profile, l_profile, significant_valleys, 0);
	}

	// Map the binned signal back.

	// Dump the valleys file.
	char valleys_bed_fp[1000];
	sprintf(valleys_bed_fp, "%s/significant_valleys_%s.bed.gz", op_dir, chr_id);
	dump_valleys(significant_valleys, valleys_bed_fp);

	delete[] log_factorials;

	delete[] derivative_map;
	delete all_extrema_regs;
	delete minima_nodes;
	delete maxima_nodes;
	return(significant_valleys);
}

void dump_valleys(vector<t_annot_region*>* significant_valleys, char* valleys_bed_fp)
{
	// Dump the valleys file.
	//char valleys_bed_fp[1000];
	//sprintf(valleys_bed_fp, "%s/significant_valleys_%s.bed", op_dir, chr_id);
	FILE* f_valleys = open_f(valleys_bed_fp, "w");

	// Write a header line?
	fprintf(f_valleys, "#CHROM\t\
START\t\
END\t\
MIN_POSN\t\
HEIGHT_AT_LEFT\t\
HEIGHT_AT_RIGHT\t\
HEIGHT_AT_MIN\t\
AVG_MMAP\t\
MAX_MMAP\t\
LEFT_HILL_SCORE\t\
RIGHT_HILL_SCORE\t\
nA\t\
nC\t\
nG\t\
nT\t\
nCpG\t\
pVal\t\
qVal\n");

	for (int i_val = 0; i_val < (int)(significant_valleys->size()); i_val++)
	{
		void** valley_info = (void**)(significant_valleys->at(i_val)->data);

		t_extrema_node** valley_points = (t_extrema_node**)(valley_info[0]);
		double* cur_valley_multimapp_signals = (double*)(valley_info[1]);
		double* left_right_pos_der_frac = (double*)(valley_info[2]);
		int* nuc_counts = (int*)(valley_info[3]);

		t_valley_significance_info* significance_info = significant_valleys->at(i_val)->significance_info;

		fprintf(f_valleys, "%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\n",
			significant_valleys->at(i_val)->chrom,
			translate_coord(significant_valleys->at(i_val)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(significant_valleys->at(i_val)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			valley_points[2]->extrema_posn,
			valley_points[0]->height_at_extrema,
			valley_points[1]->height_at_extrema,
			valley_points[2]->height_at_extrema,
			cur_valley_multimapp_signals[0],
			cur_valley_multimapp_signals[1],
			left_right_pos_der_frac[0],
			left_right_pos_der_frac[1],
			nuc_counts[0],
			nuc_counts[1],
			nuc_counts[2],
			nuc_counts[3],
			nuc_counts[5],
			significance_info->log_p_val,
			significance_info->log_q_val);
	} // i_val loop.

	close_f(f_valleys, valleys_bed_fp);
}

// Sort valleys per their minima position.
bool sort_valleys_per_min_posn(t_annot_region* reg1, t_annot_region* reg2)
{
	void** reg1_info = (void**)(reg1->data);
	void** reg2_info = (void**)(reg2->data);
	t_extrema_node** reg1_points = (t_extrema_node**)(reg1_info[0]);
	t_extrema_node** reg2_points = (t_extrema_node**)(reg2_info[0]);

	return(reg1_points[2]->extrema_posn < reg2_points[2]->extrema_posn);
}

vector<t_annot_region*>* merge_overlapping_valleys_per_pval_minimization(char* valleys_BED_fp, int l_minima_vicinity)
{
	// Read the header.
	FILE* f_valleys = open_f(valleys_BED_fp, "r");
	char* header_line = getline(f_valleys);
	if (header_line[0] != '#')
	{
		fprintf(stderr, "Could not load header.\n");
		exit(0);
	}
	vector<char*>* header_cols = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(header_line, "\t"));
	fclose(f_valleys);

	fprintf(stderr, "Read %d columns.\n", (int)header_cols->size());

	int min_posn_col_i = t_string::get_i_str(header_cols, "MIN_POSN");
	int pval_col_i = t_string::get_i_str(header_cols, "pVal");
	fprintf(stderr, "Found min posn col @ %d, p-val col @ %d\n", min_posn_col_i, pval_col_i);

	vector<t_annot_region*>* merged_significant_valleys = new vector<t_annot_region*>();

	vector<t_annot_region*>* all_valleys = load_BED_with_line_information(valleys_BED_fp);

	// Set the minimum position as the score.
	for (int i_reg = 0; i_reg < (int)all_valleys->size(); i_reg++)
	{
		char* i_reg_line = (char*)(all_valleys->at(i_reg)->data);
		t_string_tokens* toks = t_string::tokenize_by_chars(i_reg_line, "\t");
		int i_reg_min_posn = atoi(toks->at(min_posn_col_i)->str());
		all_valleys->at(i_reg)->score = i_reg_min_posn;
		t_string::clean_tokens(toks);
	} // i_regl oop.

	fprintf(stderr, "Merging %d valleys whose minima overlap within %d base pairs per p-value minimization.\n", (int)all_valleys->size(), l_minima_vicinity);

	// Divide wrt chromosomes.
	t_restr_annot_region_list* restr_valleys = restructure_annot_regions(all_valleys);

	for (int i_chr = 0; i_chr < (int)restr_valleys->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* significant_valleys = restr_valleys->regions_per_chrom[i_chr];

		fprintf(stderr, "%s: Merging %d valleys.\n", restr_valleys->chr_ids->at(i_chr), (int)significant_valleys->size());

		sort(significant_valleys->begin(), significant_valleys->end(), sort_regions_per_score);

		int i_reg = 0;
		while (i_reg < (int)significant_valleys->size())
		{
			char* i_reg_line = (char*)(significant_valleys->at(i_reg)->data);
			t_string_tokens* toks = t_string::tokenize_by_chars(i_reg_line, "\t");
			double i_reg_p_val = atof(toks->at(pval_col_i)->str());
			int i_reg_min_posn = atoi(toks->at(min_posn_col_i)->str());
			t_string::clean_tokens(toks);

			int j_reg = i_reg;

			fprintf(stderr, "@%d-%d-%d: %lf           \r", significant_valleys->at(i_reg)->start, i_reg_min_posn, significant_valleys->at(i_reg)->end, i_reg_p_val);

			t_annot_region* most_significant_valley = NULL;
			double most_significant_valley_p_val = 1000;

			// Find the valleys to merge with this valley.
			while (j_reg < (int)significant_valleys->size())
			{
				char* j_reg_line = (char*)(significant_valleys->at(j_reg)->data);
				t_string_tokens* toks = t_string::tokenize_by_chars(j_reg_line, "\t");
				double j_reg_p_val = atof(toks->at(pval_col_i)->str());
				int j_reg_min_posn = atoi(toks->at(min_posn_col_i)->str());
				t_string::clean_tokens(toks);

				if (fabs(i_reg_min_posn - j_reg_min_posn) < l_minima_vicinity)
				{
					if (i_reg_min_posn != j_reg_min_posn)
					{
						fprintf(stderr, "Nearby Valley @ %d-%d-%d: %lf                     \r",
							significant_valleys->at(j_reg)->start, j_reg_min_posn, significant_valleys->at(j_reg)->end, j_reg_p_val);
					}

					// Check if this valley is more significant than previous overlapping valleys.
					if (most_significant_valley == NULL ||
						most_significant_valley_p_val > j_reg_p_val)
					{
						most_significant_valley = significant_valleys->at(j_reg);
						most_significant_valley_p_val = j_reg_p_val;
					}
				}
				else
				{
					// Add the most significant valley.
					merged_significant_valleys->push_back(most_significant_valley);

					break;
				}

				j_reg++;
			} // j_reg loop.

			// Update i_reg.
			i_reg = j_reg;
		} // i_reg loop.

		fprintf(stderr, "\n");
	} // i_chr loop.

	fprintf(stderr, "Merged into %d valleys.\n", (int)merged_significant_valleys->size());
	return(merged_significant_valleys);
}

vector<t_annot_region*>* load_annotation(char* gff_fp, int l_half_prom)
{
	fprintf(stderr, "Loading annotations from %s\n", gff_fp);

	vector<t_annot_region*>* annotation_regions = new vector<t_annot_region*>();
	FILE* f_gff = open_f(gff_fp, "r");
	if (f_gff == NULL)
	{
		fprintf(stderr, "Could not open gff file from %s\n", gff_fp);
		exit(0);
	}

	while (1)
	{
		char* cur_line = getline(f_gff);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");

		if (t_string::compare_strings(toks->at(2)->str(), "gene") ||
			t_string::compare_strings(toks->at(2)->str(), "transcript") ||
			t_string::compare_strings(toks->at(2)->str(), "exon") ||
			t_string::compare_strings_ci(toks->at(2)->str(), "ENCODE2_TF_Peak") ||
			t_string::compare_strings_ci(toks->at(2)->str(), "tf_peak"))
		{
		}
		else
		{
			delete[] cur_line;
			t_string::clean_tokens(toks);
			continue;
		}

		t_string_tokens* info_toks = toks->at(8)->tokenize_by_chars("; =");
		char* gene_name = NULL;
		//char* transcript_name = NULL;
		//char* element_type = NULL;
		for (int i_t = 0; i_t < (int)info_toks->size(); i_t++)
		{
			if (t_string::compare_strings(info_toks->at(i_t)->str(), "gene_name") ||
				t_string::compare_strings(info_toks->at(i_t)->str(), "TF_name"))
			{
				gene_name = t_string::copy_me_str(info_toks->at(i_t + 1)->str());
			}

			//if (t_string::compare_strings(info_toks->at(i_t)->str(), "transcript_id"))
			//{
			//	transcript_name = t_string::copy_me_str(info_toks->at(i_t + 1)->str());
			//}			

			//if (t_string::compare_strings(info_toks->at(i_t)->str(), "gene_type"))
			//{
			//	element_type = t_string::copy_me_str(info_toks->at(i_t + 1)->str());
			//}
		} // i_t loop.

		t_string::clean_tokens(info_toks);

		t_episfr_annot_info* cur_annot_info = new t_episfr_annot_info();
		cur_annot_info->element_type = t_string::copy_me_str(toks->at(2)->str());
		cur_annot_info->element_name = t_string::copy_me_str(gene_name);

		char* cur_chrom = t_string::copy_me_str(toks->at(0)->str());
		normalize_chr_id(cur_chrom);

		// chr1    HAVANA  gene    11869   14412   .       +       .       ID=ENSG00000223972.4;gene_id=ENSG00000223972.4;tra
		t_annot_region* cur_element = get_empty_region();
		cur_element->chrom = t_string::copy_me_str(cur_chrom);
		cur_element->start = atoi(toks->at(3)->str());
		cur_element->end = atoi(toks->at(4)->str());
		cur_element->strand = toks->at(6)->str()[0];
		cur_element->name = gene_name;
		cur_element->data = cur_annot_info;

		if (false)
		{
			fprintf(stderr, "%s (%s): %s:%d-%d (%c)\n", gene_name, cur_annot_info->element_type, cur_chrom, cur_element->start, cur_element->end, cur_element->strand);
		}

		annotation_regions->push_back(cur_element);

		// Add the transcript promoters.
		if (t_string::compare_strings(cur_annot_info->element_type, "transcript"))
		{
			t_annot_region* prom_element = get_empty_region();
			prom_element->chrom = t_string::copy_me_str(cur_chrom);
			prom_element->start = atoi(toks->at(3)->str());
			prom_element->end = atoi(toks->at(4)->str());
			prom_element->strand = toks->at(6)->str()[0];
			prom_element->name = t_string::copy_me_str(gene_name);

			t_episfr_annot_info* prom_annot_info = new t_episfr_annot_info();
			prom_annot_info->element_type = t_string::copy_me_str("promoter");
			prom_annot_info->element_name = t_string::copy_me_str(gene_name);

			int prom_start = prom_element->start - l_half_prom;
			int prom_end = prom_element->start + l_half_prom;

			if (prom_element->strand == '-')
			{
				prom_start = prom_element->end - l_half_prom;
				prom_end = prom_element->end + l_half_prom;
			}

			prom_element->start = prom_start;
			prom_element->end = prom_end;
			prom_element->data = prom_annot_info;

			if (false)
			{
				fprintf(stderr, "%s (%s): %s:%d-%d (%c)\n", gene_name, cur_annot_info->element_type, cur_chrom, prom_element->start, prom_element->end, cur_element->strand);
			}

			annotation_regions->push_back(prom_element);
		}

		if ((int)annotation_regions->size() % (100 * 1000) == 0)
		{
			fprintf(stderr, "@ %d. element..\r", (int)annotation_regions->size());
		}

		t_string::clean_tokens(toks);
		delete[] cur_line;
	} // gff reading loop.	

	// Close the bedgraph file.
	close_f(f_gff, gff_fp);

	return(annotation_regions);
}

void annotate_features(char* valleys_bed_fp, char* gff_fp, int l_half_prom, char* op_fp)
{
	// Load the gff file.
	vector<t_annot_region*>* annotations = load_annotation(gff_fp, l_half_prom);
	fprintf(stderr, "Loaded %d annotation elements.\n", (int)annotations->size());

	// Load the valleys.
	vector<t_annot_region*>* valley_regs = load_BED_with_line_information(valleys_bed_fp);
	fprintf(stderr, "Loaded %d valleys.\n", (int)valley_regs->size());

	// Replace info's.
	for (int i_val = 0; i_val < (int)valley_regs->size(); i_val++)
	{
		void* cur_info = valley_regs->at(i_val)->data;
		void** new_info = new void*[2];
		new_info[0] = cur_info;
		new_info[1] = new vector<t_episfr_annot_info*>();

		valley_regs->at(i_val)->data = new_info;
	} // i_val loop.

	// Annotate.
	vector<t_annot_region*>* intersects = intersect_annot_regions(valley_regs, annotations, true);
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);

		t_annot_region* int_valley_reg = int_info->src_reg;
		t_annot_region* int_annot_reg = int_info->dest_reg;

		// Write the annotation for the current valley.
		t_episfr_annot_info* cur_annot_info = (t_episfr_annot_info*)(int_annot_reg->data);


		void** cur_valley_info = (void**)(int_valley_reg->data);
		vector<t_episfr_annot_info*>* annots = (vector<t_episfr_annot_info*>*)(cur_valley_info[1]);
		annots->push_back(cur_annot_info);

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

	// Read the header.
	FILE* f_valleys = open_f(valleys_bed_fp, "r");
	char* header_line = getline(f_valleys);
	fclose(f_valleys);

	// Save header then write annotated regions.
	FILE* f_op = open_f(op_fp, "w");

	// Add the header line.
	if (header_line[0] == '#')
	{
		char* gff_fn = get_file_name(gff_fp);
		fprintf(f_op, "%s\t%s_Annotation\n", header_line, gff_fn);
	}

	for (int i_val = 0; i_val < (int)valley_regs->size(); i_val++)
	{
		void** cur_valley_info = (void**)(valley_regs->at(i_val)->data);
		char* cur_valley_line = (char*)(cur_valley_info[0]);
		vector<t_episfr_annot_info*>* annots = (vector<t_episfr_annot_info*>*)(cur_valley_info[1]);

		fprintf(f_op, "%s\t", cur_valley_line);
		if (annots->size() == 0)
		{
			fprintf(f_op, ".");
		}
		else
		{
			// Process all the annotations from this file.
			for (int i_ann = 0; i_ann < (int)annots->size(); i_ann++)
			{
				fprintf(f_op, "%s:%s;", annots->at(i_ann)->element_type, annots->at(i_ann)->element_name);
			} // i_ann loop.
		}

		fprintf(f_op, "\n");
	} // i_val loop.

	fclose(f_op);
}

// Select changing points to use in the polyfit.
void select_points_of_interest_per_RD_signal_profile(double* signal_profile, int l_profile,
	int start_i, int end_i,
	int min_dist_between_cons_pts,
	vector<double>* x_vec, vector<double>* y_vec,
	bool sparse_signal)
{
	//fprintf(stderr, "Selecting POI for signal track in [%d-%d] with maximum distance of points %d\n", start_i, end_i, max_dist_between_cons_pts);

	//// Update the points at each change in signal.
	//double cur_y = signal_profile[start_i];
	//double cur_x = start_i;

	//if (!sparse_signal)
	//{
	//	x_vec->push_back(start_i);
	//}

	if (!sparse_signal)
	{
#undef __USE_PLATEAU_MID_POIS__
#ifdef __USE_PLATEAU_MID_POIS__
		int l_win = end_i - start_i;
		int first_plateau_start_pos = MAX(1, start_i - l_win);
		int last_plateau_end_pos = MIN(l_profile, end_i + l_win);

		// Add mid plateaus.
		int i = first_plateau_start_pos;
		while (i < last_plateau_end_pos)
		{
			// Find the mid point.			
			int plateau_start = i;
			int plateau_end = i;
			while (plateau_end < last_plateau_end_pos &&
				signal_profile[plateau_end] == signal_profile[i])
			{
				plateau_end++;
			}
			i = plateau_end;

			int plateau_mid_pt = (plateau_start + (plateau_end - 1)) / 2;

			if (plateau_mid_pt >= start_i && plateau_mid_pt < end_i)
			{
				x_vec->push_back(plateau_mid_pt);
			}
		} // i loop.
#endif
		// Add the periodic POIs: Make sure to add first point and not the last.
		for(int i = start_i; i < end_i; i+= min_dist_between_cons_pts)
		{
			bool added_per_height_change = false;
			for (int i_v = 0; i_v < (int)x_vec->size(); i_v++)
			{
				if (x_vec->at(i_v) == i)
				{
					added_per_height_change = true;
					break;
				}
			} // i_v loop.

			if (!added_per_height_change)
			{
				x_vec->push_back(i);
			}

			//fprintf(stderr, "%d-%d: Adding %d (cur_x=%.3f)\n", start_i, end_i, i, cur_x);
			//getc(stdin);
		} // i loop.

		sort(x_vec->begin(), x_vec->end());

		// Add heights.
		for (int i_v = 0; i_v < (int)x_vec->size(); i_v++)
		{
			// Make sure there are no duplicates.
			if (i_v > 0)
			{
				if (x_vec->at(i_v) == x_vec->at(i_v - 1))
				{
					fprintf(stderr, "Duplicated POIs.\n");
					exit(0);
				}
			}

			y_vec->push_back(signal_profile[(int)(x_vec->at(i_v))]);
		} // i_v loop.
	}
	else
	{
		for (int i = start_i; i <= end_i; i++)
		{
			// If sparse signal is requested, make sure the current signal is not 0, i.e. sparsity. We do not want to use missing regions for spling-fitting.
			if (signal_profile[i] != 0)
			{
				// Update the current point.
				x_vec->push_back(i);
				y_vec->push_back(signal_profile[i]);
			}
			else
			{
				// Do nothing: Do not update the data vectors.
			}
		} // i loop.
	}

	//for (int i = start_i; i <= end_i; i++)
	//{
	//	// If sparse signal is requested, make sure the current signal is not 0, i.e. sparsity. We do not want to use missing regions for spling-fitting.
	//	bool sparsity_check = (!sparse_signal || (sparse_signal && signal_profile[i] != 0));

	//	// If sparsity check holds, 
	//	if (sparsity_check &&
	//		(sparse_signal ||
	//			cur_y != signal_profile[i] ||
	//			cur_x + min_dist_between_cons_pts <= i))
	//	{
	//		// Update the current point.
	//		cur_y = signal_profile[i];
	//		cur_x = i;

	//		x_vec->push_back(i);
	//		y_vec->push_back(signal_profile[i]);
	//	}
	//	else
	//	{
	//		// Do nothing: Do not update the data vectors.
	//	}
	//} // i loop.

	//if (!sparse_signal)
	//{
	//	x_vec->push_back(end_i);
	//	y_vec->push_back(signal_profile[end_i]);
	//}

}


// Select changing points to use in the polyfit.
void select_points_of_interest_per_RD_signal_profile_OBSOLETE(double* signal_profile,
	int start_i, int end_i,
	int min_dist_between_cons_pts,
	vector<double>* x_vec, vector<double>* y_vec,
	bool sparse_signal)
{
	//fprintf(stderr, "Selecting POI for signal track in [%d-%d] with maximum distance of points %d\n", start_i, end_i, max_dist_between_cons_pts);

	// Update the points at each change in signal.
	double cur_y = signal_profile[start_i];
	double cur_x = start_i;

	// Add the first position in the region.
	if (!sparse_signal)
	{
		x_vec->push_back(start_i);
		y_vec->push_back(cur_y);
	}

	if (!sparse_signal)
	{
		for (int i = start_i + 1; i < end_i; i++)
		{
			bool added_per_height_change = false;

#undef __USE_HEIGHT_CHANGES__
#ifdef __USE_HEIGHT_CHANGES__
			if (cur_y != signal_profile[i])
			{
				added_per_height_change = true;
				x_vec->push_back(i);
				y_vec->push_back(signal_profile[i]);
			}
#endif // __USE_HEIGHT_CHANGES__

			if (cur_x + min_dist_between_cons_pts == i)
			{
				// Update the current point.
				cur_x = i;

				if (!added_per_height_change)
				{
					x_vec->push_back(i);
					y_vec->push_back(signal_profile[i]);
				}

				//fprintf(stderr, "%d-%d: Adding %d (cur_x=%.3f)\n", start_i, end_i, i, cur_x);
				//getc(stdin);
			}
			else
			{
				// Do nothing: Do not update the data vectors.
			}

			cur_y = signal_profile[i];
		} // i loop.
	}
	else
	{
		for (int i = start_i; i <= end_i; i++)
		{
			// If sparse signal is requested, make sure the current signal is not 0, i.e. sparsity. We do not want to use missing regions for spling-fitting.
			bool sparsity_check = (!sparse_signal || (sparse_signal && signal_profile[i] != 0));

			// If sparsity check holds, 
			if (sparsity_check &&
				(sparse_signal ||
					cur_y != signal_profile[i] ||
					cur_x + min_dist_between_cons_pts <= i))
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

	//for (int i = start_i; i <= end_i; i++)
	//{
	//	// If sparse signal is requested, make sure the current signal is not 0, i.e. sparsity. We do not want to use missing regions for spling-fitting.
	//	bool sparsity_check = (!sparse_signal || (sparse_signal && signal_profile[i] != 0));

	//	// If sparsity check holds, 
	//	if (sparsity_check &&
	//		(sparse_signal ||
	//			cur_y != signal_profile[i] ||
	//			cur_x + min_dist_between_cons_pts <= i))
	//	{
	//		// Update the current point.
	//		cur_y = signal_profile[i];
	//		cur_x = i;

	//		x_vec->push_back(i);
	//		y_vec->push_back(signal_profile[i]);
	//	}
	//	else
	//	{
	//		// Do nothing: Do not update the data vectors.
	//	}
	//} // i loop.

	//if (!sparse_signal)
	//{
	//	x_vec->push_back(end_i);
	//	y_vec->push_back(signal_profile[end_i]);
	//}

}

//#define __COMBINE_PER_SUM__
#define __COMBINE_PER_MAX__
//#define __COMBINE_PER_MIN__
//#define __COMBINE_PER_SELF__

#define COMBINE_SUM(x, y) ((x)+(y))
#define COMBINE_MAX(x, y) (((x)>(y))?(x):(y))

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
										int l_med_filt_win)
{
	double* coeff = new double[n_spline_coeff + 2];

	bool reads_loaded = false;
	int l_track = 0;
	double* signal_profile = load_signal_covg_per_directory_chr_id(signal_dir, chr_id, l_frag, l_track, reads_loaded);
	if (signal_profile == NULL)
	{
		fprintf(stderr, "Could not load the coverage from %s/%s\n", signal_dir, chr_id);
		exit(0);
	}

	if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Saving the original bedgraph file.\n");
		char orig_bgr_fp[1000];
		sprintf(orig_bgr_fp, "%s/original_%s.bgr.gz", signal_dir, chr_id);
		dump_bedGraph_per_per_nucleotide_binary_profile(signal_profile, l_track, chr_id, orig_bgr_fp);
	}

	double* spline_fit_profile = new double[l_track + 2];
	memset(spline_fit_profile, 0, sizeof(double) * (l_track+2));

	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	double aggregate_total_err = 0;
	int start_i = 1;
	while (start_i < l_track)
	{
		double tot_sig = 0;
		for (int i = start_i; i < start_i + l_win; i++)
		{
			if (i < l_track)
			{
				tot_sig += signal_profile[i];
			}
		} // i loop.

		// If there is not signal in this window, do not process.
		if (tot_sig == 0)
		{
			start_i += l_step_win;
			continue;
		}

		int end_i = MIN(l_track-1, start_i + l_win);
		vector<double>* x_vec = new vector<double>();
		vector<double>* y_vec = new vector<double>();
		select_points_of_interest_per_RD_signal_profile(signal_profile, l_track,
			start_i, end_i,
			max_dist_between_cons_pts,
			x_vec, y_vec, sparse_profile);

		if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
		{
			FILE* f_POI = open_f("POIs.txt", "a");
			for (int i = 0; i < (int)x_vec->size(); i++)
			{
				fprintf(f_POI, "%.1f\n", x_vec->at(i));
			}
			fclose(f_POI);
		}

		// Allocate and generate the x and y data.
		//fprintf(stderr, "Generated %d data points.\n", x_vec->size());
		int n_data_points = (int)x_vec->size();

		// Make sure there are at least 5 points, otherwise there is not much change here to encode.
		if (n_data_points > min_n_pts_2_encode)
		{
			double* x = new double[n_data_points + 2];
			double* y = new double[n_data_points + 2];
			double* reconst_y = new double[n_data_points + 2];
			memset(reconst_y, 0, sizeof(double) * (n_data_points + 1));

			double top_perc_err = 0;

			double total_err = 0;
			double max_err = 0;
			int cur_n_spline_coeff = n_spline_coeff;
			int cur_bspline_order = bspline_order;
			double total_POI_signal = 0;

			// Copy vectors.
			for (int i_l = 0; i_l < (int)x_vec->size(); i_l++)
			{
				total_POI_signal += y_vec->at(i_l);

				x[i_l] = (x_vec->at(i_l) - start_i);
				y[i_l] = y_vec->at(i_l);
			} // i_l loop.

if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
{
			char data_win_fp[1000];
			sprintf(data_win_fp, "data_%d.txt", start_i);
			FILE* f_input = open_f(data_win_fp, "w");
			for (int i_l = 0; i_l < (int)x_vec->size(); i_l++)
			{
				fprintf(f_input, "%lf\t%lf\n", x[i_l], y[i_l]);
			} // i_l loop.
			fclose(f_input);
}

			// This is the error tracking loop.
			do
			{
				if (n_data_points < cur_n_spline_coeff)
				{
					break;
				}

				// Reset the errors.
				total_err = 0;
				max_err = 0;

				//// Do fitting.				
				//if (brkpts_type == UNIFORM_BREAKPOINTS)
				//{
				//	bsplinefit(n_data_points, cur_n_spline_coeff, cur_bspline_order, x, y, reconst_y, coeff);
				//}
				//else
				//{
				bsplinefit_nonuniform_per_breakpoint_type(n_data_points,
					cur_n_spline_coeff - cur_bspline_order, // This gives us the number of internal breakpoints.
					brkpts_type,
					cur_bspline_order,
					x, y,
					reconst_y,
					coeff,
					rng,
					false);
				//}
				//else if (HILL_DERIVATIVE_NU_BREAKPOINTS)
				//{
				//	bsplinefit_nonuniform_per_hill_derivative(n_data_points,
				//		cur_n_spline_coeff - cur_bspline_order,
				//		cur_bspline_order,
				//		x, y,
				//		reconst_y,
				//		coeff);
				//}

				// compute the errors only on the fit locations.
				vector<double>* errors = new vector<double>();
				for (int i_l = 0; i_l < n_data_points; i_l++)
				{
					double cur_reconst_y = reconst_y[i_l];

					if (cur_reconst_y < 0)
					{
						cur_reconst_y = 0;
					}

					double cur_err = fabs(cur_reconst_y - y[i_l]);
					total_err += cur_err;
					max_err = MAX(max_err, cur_err);
					errors->push_back(cur_err);
				} // i_l loop.

				sort(errors->begin(), errors->end());

				int top_err_perc_index = (int)(errors->size() * top_err_perc_frac);
				if (top_err_perc_index >= (int)errors->size())
				{
					top_err_perc_index = (int)errors->size() - 1;
				}

				top_perc_err = errors->at(top_err_perc_index);
				delete errors;

				  // Save reconstructured window.
if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
{
					char recons_fp[1000];
					sprintf(recons_fp, "recon_%d.txt", start_i);
					FILE* f_reconst = open_f(recons_fp, "w");
					for (int i_l = 0; i_l < n_data_points; i_l++)
					{
						double cur_reconst_y = reconst_y[i_l];
						fprintf(f_reconst, "%lf\t%lf\n", x[i_l], cur_reconst_y);
					} // i_l loop.
					fclose(f_reconst);
					//fprintf(stderr, "\n");
}

				// Add to the spline fit profile 
				for (int i_l = 0; i_l < n_data_points - 1; i_l++)
				{
					// copy this block.
					for (int i = (int)x[i_l]; i < (int)x[i_l + 1]; i++)
					{
#ifdef __COMBINE_PER_SUM__
						spline_fit_profile[i + start_i] = COMBINE_SUM(reconst_y[i_l], spline_fit_profile[i + start_i]);
#endif

#ifdef __COMBINE_PER_MAX__
						spline_fit_profile[i + start_i] = COMBINE_MAX(reconst_y[i_l], spline_fit_profile[i + start_i]);
#endif
					} // i loop.
				} // i_l loop.

				// To the left of first control point and to the right of the last control point.
				for (int i = 0; i < (int)x[0]; i++)
				{
#ifdef __COMBINE_PER_SUM__
					spline_fit_profile[i + start_i] = COMBINE_SUM(reconst_y[0], spline_fit_profile[i + start_i]);
#endif

#ifdef __COMBINE_PER_MAX__
					spline_fit_profile[i + start_i] = COMBINE_MAX(reconst_y[0], spline_fit_profile[i + start_i]);
#endif
					//spline_fit_profile[i + start_i] = reconst_y[0];
				} // i loop.

				for (int i = (int)x[n_data_points-1]; i < l_win; i++)
				{
					if (i + start_i < l_track)
					{
#ifdef __COMBINE_PER_SUM__
						spline_fit_profile[i + start_i] = COMBINE_SUM(reconst_y[n_data_points - 1], spline_fit_profile[i + start_i]);
#endif 

#ifdef __COMBINE_PER_MAX__
						spline_fit_profile[i + start_i] = COMBINE_MAX(reconst_y[n_data_points - 1], spline_fit_profile[i + start_i]);
#endif
				  		//spline_fit_profile[i + start_i] = reconst_y[n_data_points - 1];
					}
				} // i loop.

				fprintf(stderr, "Window %d-%d: %d coefficients, spline order %d, %d POIs (Total %.2f): Error: %lf, Avg Error: %lf, Top Error @ %d: %lf                 \r",
					start_i, end_i,
					cur_n_spline_coeff,
					cur_bspline_order,
					(int)x_vec->size(),
					total_POI_signal,
					total_err, 
					total_err / n_data_points, 
					top_err_perc_index, top_perc_err);

				// Increase the # of coefficients and the bpline order. This tunes the knot points that are uniformly distributed. This part can be potentially optimized.
				cur_n_spline_coeff += 2;

				// Do not change the order, this may create overfitting, which is not what we want.
				//cur_bspline_order += 2;
			} // spline fitting loop.
			while (total_err / n_data_points > max_avg_err ||
				top_perc_err > max_top_err);

			aggregate_total_err += total_err;

			delete[] reconst_y;
			delete[] x;
			delete[] y;
		}
		else
		{
if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
{
			fprintf(stderr, "No fit @ %s:%d-%d\n", chr_id, start_i, end_i);
}

			// Copy the signal profile in this window directly.
			for (int i = start_i; i < end_i; i++)
			{
				spline_fit_profile[i] = signal_profile[i];
			} // i loop.
		}

		delete x_vec;
		delete y_vec;

		start_i += l_step_win;
		//start_i += l_win;
	} // l_win loop.

#ifdef __COMBINE_PER_SUM__
	double norm_fact = (double)l_win / (double)l_step_win;
	fprintf(stderr, "\nNormalizing factor is %.2f\n", norm_fact);
	for (int i = 1; i <= l_track; i++)
	{
		spline_fit_profile[i] /= norm_fact;
	} // i loop.
#elif defined __COMBINE_PER_MAX__
	// Do nothing.
	//double norm_fact = 1.0;
#endif 

	double absolute_aggregate_total_err = 0;
	for (int i = 1; i < l_track; i++)
	{
		absolute_aggregate_total_err += fabs(spline_fit_profile[i] - signal_profile[i]);
	} // i loop.

	fprintf(stderr, "\nAggregate total error: %lf\nAbsolute error: %lf\n", aggregate_total_err, absolute_aggregate_total_err);

	// Write the error summary.
	FILE* f_errors = open_f("error_summary.txt", "a");
	fprintf(f_errors, "%s\t%lf\t%lf\t%d\n", chr_id, aggregate_total_err, absolute_aggregate_total_err, l_track);
	fclose(f_errors);

	// Median smooth the data.
	if (l_med_filt_win > 0)
	{
		fprintf(stderr, "Smoothing spline coded profile using window length of %d.\n", l_med_filt_win);

		//double min_val = 1000 * 1000 * 1000;
		//for (int i = 1; i < l_track; i++)
		//{
		//	if (spline_fit_profile[i] < min_val)
		//	{
		//		min_val = spline_fit_profile[i];
		//	}
		//} // i loop.

		for (int i = 1; i < l_track; i++)
		{
			//spline_fit_profile[i] -= min_val;
			spline_fit_profile[i] *= 100;
		}

		double* med_smoothed_profile = median_filter_data(spline_fit_profile,
														l_track,
														l_med_filt_win,
														-1000);

		for (int i = 1; i < l_track; i++)
		{
			med_smoothed_profile[i] /= 100;
			//med_smoothed_profile[i] += min_val;
		}

		delete[] spline_fit_profile;
		spline_fit_profile = med_smoothed_profile;
	}

	fprintf(stderr, "Finished, saving spline coded profile.\n");
	char encoded_profile_bin_fp[1000];
	sprintf(encoded_profile_bin_fp, "%s/spline_coded_%s.bin.gz", signal_dir, chr_id);
	dump_per_nucleotide_binary_profile(spline_fit_profile, l_track, encoded_profile_bin_fp);

	if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
	{
		char spline_bgr_fp[1000];
		sprintf(spline_bgr_fp, "%s/spline_coded_%s.bgr.gz", signal_dir, chr_id);
		dump_bedGraph_per_per_nucleotide_binary_profile(spline_fit_profile, l_track, chr_id, spline_bgr_fp);
	}

	delete[] coeff;
	delete[] signal_profile;
	delete[] spline_fit_profile;
}

