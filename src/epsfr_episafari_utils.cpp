#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "epsfr_episafari_utils.h"
#include "epsfr_xlog_math.h"
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

//#define MAX(x, y) (((x)>(y))?(x):(y))
//#define MIN(x, y) (((x)<(y))?(x):(y))

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

double* load_signal_covg_per_directory_chr_id(char* dat_dir,
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
		covg_signal = new double[l_buffer + 2];

		buffer_per_nucleotide_profile_no_buffer(cur_dat_fp, l_fragment,
			covg_signal, NULL, NULL,
			l_buffer, l_loaded_covg);

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

// Compute binomial p-value of the signal imbalance on the valley.
double get_valley_significance(double* signal_profile, int l_profile, int left_max_posn, int right_max_posn, int min_posn, int l_normalizer, int scaling_factor, double* _log_factorials)
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
		fprintf(stderr, "Valley: %d-%d-%d: %lf, %lf: %lf (flip_prob: %.3f)\n", left_max_posn, min_posn, right_max_posn, norm_min_vic_sig, norm_left_max_vic_sig, norm_right_max_vic_sig, exp(log_flip));
	}

	// Free memory if it is allocated in the function.
	if (_log_factorials == NULL)
	{
		delete[] log_factorials;
	}

	return(right_region_log_p_val + left_region_log_p_val);
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

// Compute binomial p-value of the signal imbalance on the valley.
double get_valley_significance_OBS(double* signal_profile, int l_profile, int left_max_posn, int right_max_posn, int min_posn, int l_frag, double* _log_factorials)
{
	// Get the total signal and left and right signal levels, too.
	double total_sig = 0;
	double left_sig = 0;
	double right_sig = 0;
	for (int pos = left_max_posn - l_frag; pos <= right_max_posn + l_frag; pos++)
	{
		total_sig += signal_profile[pos];

		if (pos <= min_posn)
		{
			left_sig += signal_profile[pos];
		}

		if (pos >= min_posn)
		{
			right_sig += signal_profile[pos];
		}
	} // pos loop.

	int signal_at_min = (int)(signal_profile[min_posn]);

	double total_norm_signal = total_sig / l_frag;
	double left_norm_signal = left_sig / l_frag;
	double right_norm_signal = right_sig / l_frag;

	double right_flip_prob = (double)(right_max_posn - min_posn) / (double)(right_max_posn - left_max_posn + 1);
	double left_flip_prob = (double)(min_posn - left_max_posn) / (double)(right_max_posn - left_max_posn + 1);
	right_flip_prob = MAX(0.0001, right_flip_prob);
	left_flip_prob = MAX(0.0001, left_flip_prob);

	int grand_total_int = (int)(total_norm_signal);

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
	double log_right_flip = log(right_flip_prob);
	double log_left_flip = log(left_flip_prob);
	double cur_region_log_p_val = xlog(0.0);

	// Move the signal at minimum to left side.
	for (int i = 0; i <= signal_at_min; i++)
	{
		int cur_left_total_sig = (int)left_norm_signal + i;
		int cur_right_total_sig = grand_total_int - cur_left_total_sig;

		double log_cur_half_pow = log_left_flip * cur_left_total_sig + cur_right_total_sig * log_right_flip;

		double log_cur_perm = 0.0; // = xlog(1.0).

		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[cur_left_total_sig], log_factorials[cur_right_total_sig]));

		// Compute the current permutation.
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	// Move the signal at minimum to right side; this has to be done twice bc of non-symmetry of flips.
	for (int i = 0; i <= signal_at_min; i++)
	{
		int cur_left_total_sig = (int)left_norm_signal - i;
		int cur_right_total_sig = grand_total_int - cur_left_total_sig;

		double log_cur_half_pow = log_left_flip * cur_left_total_sig + cur_right_total_sig * log_right_flip;

		double log_cur_perm = 0.0; // = xlog(1.0).

		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[cur_left_total_sig], log_factorials[cur_right_total_sig]));

		// Compute the current permutation.
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
	{
		fprintf(stderr, "Valley: %d-%d-%d: %lf, %lf: %lf (l/r flip_prob: %.3f, %.3f)\n", left_max_posn, min_posn, right_max_posn, left_norm_signal, right_norm_signal, cur_region_log_p_val, exp(log_left_flip), exp(log_right_flip));
	}

	//// Add the probability of enrichment on the left side.
	//flip_prob = 1 - flip_prob;
	//flip_prob = MAX(0.0001, flip_prob);
	//log_flip = log(flip_prob);
	//log_one_min_flip = log(1 - flip_prob);
	//for (int i = (int)left_norm_signal; i <= grand_total_int; i++)
	//{
	//	double log_cur_half_pow = log_flip * i + (grand_total_int - i) * log_one_min_flip;
	//	double log_cur_perm = 0.0; // = xlog(1.0).

	//							   // Compute the current permutation.
	//	log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int - i]));
	//	cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	//} // i loop.

	//if (__DUMP_EPISAFARI_UTILS_MESSAGES__)
	//{
	//	fprintf(stderr, "Valley: %d-%d-%d: %lf, %lf: %lf (flip_prob: %.3f)\n", left_max_posn, min_posn, right_max_posn, left_norm_signal, right_norm_signal, cur_region_log_p_val, exp(log_flip));
	//}

	  // Free memory if it is allocated in the function.
	if (_log_factorials == NULL)
	{
		delete[] log_factorials;
	}

	return(cur_region_log_p_val);
}

//vector<t_extrema_node*>* prune_nodes_per_vicinity_signal(vector<t_extrema_node*>* minima_nodes, vector<t_extrema_node*>* maxima_nodes, t_extrema_statistic_defition* extrema_statistic_defn)
//{
//	vector<t_extrema_node*>* all_nodes = new vector<t_extrema_node*>();
//	all_nodes->insert(all_nodes->end(), minima_nodes->begin(), minima_nodes->end());
//	all_nodes->insert(all_nodes->end(), maxima_nodes->begin(), maxima_nodes->end());
//
//	sort(all_nodes->begin(), all_nodes->end(), sort_extremas_per_posn);
//
//	for (int i_ext = 0; i_ext < all_nodes->size(); i_ext++)
//	{
//
//	} // i_ext loop.
//}

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
																	double* multimapp_signal_profile, int l_multimapp_profile,
																	char* chrom_seq,
																	t_extrema_statistic_defition* extrema_statistic_defn)
{
	// Remove the negative values.
	for (int i = 0; i <= l_profile; i++)
	{
		signal_profile[i] = MAX(0, signal_profile[i]);
	} // i loop.

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

							if (cur_sig_diff > 0)
							{
								cur_dir_total_pos_signal_increase += cur_sig_diff;
							}

							if (cur_sig_diff < 0)
							{
								cur_dir_total_pos_signal_decrease += -1*cur_sig_diff;
							}
						} // i loop.

						left_right_pos_der_frac[dir] = cur_dir_total_pos_deriv_signs / (cur_dir_total_pos_deriv_signs + cur_dir_total_neg_deriv_signs);
						//left_right_pos_der_frac[dir] = cur_dir_total_pos_signal_increase / (cur_dir_total_pos_signal_increase + cur_dir_total_pos_signal_decrease);
					} // dir loop.

					// Compute the average multi-map signal within the valley.
					double total_multimapp_signal = 0;
					double max_multimapp_signal = -1;
					for (int i = cur_left_maxima->extrema_posn; i < cur_right_maxima->extrema_posn; i++)
					{
						if (i < l_multimapp_profile)
						{
							if (multimapp_signal_profile[i] > max_multimapp_signal)
							{
								max_multimapp_signal = multimapp_signal_profile[i];
							}

							total_multimapp_signal += multimapp_signal_profile[i];
						}
					} // i loop.

					// Count the nucleotides.
					int* nuc_counts = new int[10];
					memset(nuc_counts, 0, sizeof(int) * 10);
					int n_CpGs = 0;
					for (int i = cur_left_maxima->extrema_posn; i < cur_right_maxima->extrema_posn; i++)
					{
						nuc_counts[nuc_2_num(chrom_seq[i])]++;

						if (toupper(chrom_seq[i]) == 'C' &&
							toupper(chrom_seq[i + 1]) == 'G')
						{
							n_CpGs++;
						}
					} // i loop.

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
					sig_info->log_p_val = get_valley_significance(signal_profile, l_profile, cur_left_maxima->extrema_posn, cur_right_maxima->extrema_posn, cur_minima->extrema_posn, l_vic, scaling_factor, log_factorials);
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
	sprintf(valleys_bed_fp, "%s/significant_valleys_%s.bed", op_dir, chr_id);
	dump_valleys(significant_valleys, valleys_bed_fp);

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
	//fclose(f_valleys);
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
						fprintf(stderr, "\nNearby Valley @ %d-%d-%d: %lf            \n",
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
	FILE* f_gff = NULL;

	if (t_string::ends_with(gff_fp, ".gz") ||
		t_string::ends_with(gff_fp, ".gzip"))
	{
		char ungzip_cmd[1000];
		sprintf(ungzip_cmd, "gzip -cd %s", gff_fp);
#ifdef _WIN32
		f_gff = _popen(ungzip_cmd, "r");
#else 
		f_gff = popen(ungzip_cmd, "r");
#endif		
	}
	else
	{
		f_gff = open_f(gff_fp, "r");
	}

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
	if (t_string::compare_strings(gff_fp, "stdin"))
	{
	}
	else if (t_string::ends_with(gff_fp, ".gff.gz") ||
		t_string::ends_with(gff_fp, ".gff.gzip"))
	{
#ifdef _WIN32
		_pclose(f_gff);
#else 
		pclose(f_gff);
#endif
	}
	else
	{
		fclose(f_gff);
	}

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
void select_points_of_interest_per_RD_signal_profile(double* signal_profile,
	int start_i, int end_i,
	int min_dist_between_cons_pts,
	vector<double>* x_vec, vector<double>* y_vec,
	bool sparse_signal)
{
	//fprintf(stderr, "Selecting POI for signal track in [%d-%d] with maximum distance of points %d\n", start_i, end_i, max_dist_between_cons_pts);

	// Update the points at each change in signal.
	double cur_y = signal_profile[start_i];
	double cur_x = start_i;
	//x_vec->push_back(start_i);
	//y_vec->push_back(cur_y);
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

void bspline_encode_mapped_read_profile(char* signal_dir,
										char* chr_id,
										int l_frag,
										int n_spline_coeff,
										int bspline_order,
										int min_n_pts_2_encode,
										int max_dist_between_cons_pts,
										double max_top_err, 
										double max_avg_err,
										int l_win,
										bool sparse_profile,
										double top_err_perc_frac,
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

	fprintf(stderr, "Saving the original bedgraph file.\n");
	char orig_bgr_fp[1000];
	sprintf(orig_bgr_fp, "%s/original_%s.bgr", signal_dir, chr_id);
	dump_bedGraph_per_per_nucleotide_binary_profile(signal_profile, l_track, chr_id, orig_bgr_fp);

	char comp_orig_bgr_fp[1000];
	sprintf(comp_orig_bgr_fp, "%s/original_%s.bgr.gz", signal_dir, chr_id);
	compressFile(orig_bgr_fp, comp_orig_bgr_fp);

	// Delete the uncompressed file.
#ifdef _WIN32
	_unlink(orig_bgr_fp);
#endif

#ifdef __unix__
	unlink(orig_bgr_fp);
#endif

	//if (!check_file(mapped_reads_fp))
	//{
	//	fprintf(stderr, "Could not find data file @ %s\n", mapped_reads_fp);
	//	exit(0);
	//}

	//int l_frag = 200;
	//int l_buff = 250 * 1000 * 1000;
	//int l_track = 0;
	//double* signal_profile = new double[l_buff];
	//buffer_per_nucleotide_profile_no_buffer(mapped_reads_fp, l_frag, signal_profile, NULL, NULL, l_buff, l_track);

	double* spline_fit_profile = new double[l_track + 2];
	memset(spline_fit_profile, 0, sizeof(double) * (l_track+2));

	//int max_dist_between_cons_pts = 50;

	//double top_err_perc_frac = 0.95;

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
			start_i += l_win;
			continue;
		}

		int end_i = MIN(l_track-1, start_i + l_win);
		vector<double>* x_vec = new vector<double>();
		vector<double>* y_vec = new vector<double>();
		select_points_of_interest_per_RD_signal_profile(signal_profile,
			start_i, end_i,
			max_dist_between_cons_pts,
			x_vec, y_vec, sparse_profile);

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

			// Copy vectors.
			for (int i_l = 0; i_l < (int)x_vec->size(); i_l++)
			{
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

			while (total_err == 0 ||
				total_err / n_data_points > max_avg_err ||
				top_perc_err > max_top_err)
			{
				if (n_data_points < cur_n_spline_coeff)
				{
					break;
				}

				// Reset the errors.
				total_err = 0;
				max_err = 0;

				// Do fitting.
				bsplinefit(n_data_points, cur_n_spline_coeff, cur_bspline_order, x, y, reconst_y, coeff);

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
					for (int i = (int)x[i_l]; i <= (int)x[i_l + 1]; i++)
					{
						spline_fit_profile[i + start_i] = reconst_y[i_l];
					} // i loop.
				} // i_l loop.

				// To the left of first control point and to the right of the last control point.
				for (int i = 0; i <= (int)x[0]; i++)
				{
					spline_fit_profile[i + start_i] = reconst_y[0];
				} // i loop.

				for (int i = (int)x[n_data_points-1]; i < l_win; i++)
				{
					if (i + start_i < l_track)
					{
						spline_fit_profile[i + start_i] = reconst_y[n_data_points - 1];
					}
				} // i loop.

				fprintf(stderr, "Window %d-%d: %d coefficients, spline order %d, %d POIs: Error: %lf, Avg Error: %lf, Top Error @ %d: %lf                 \r",
					start_i, end_i,
					cur_n_spline_coeff,
					cur_bspline_order,
					(int)x_vec->size(),
					total_err, 
					total_err / n_data_points, 
					top_err_perc_index, top_perc_err);

				// Increase the # of coefficients and the bpline order. This tunes the knot points that are uniformly distributed. This part can be potentially optimized.
				cur_n_spline_coeff += 10;

				// Do not change the order, this may create overfitting, which is not what we want.
				//cur_bspline_order += 2;
			} // spline fitting loop.

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

		start_i += l_win;
	} // l_win loop.

	// Median smooth the data.
	if (l_med_filt_win > 0)
	{
		fprintf(stderr, "\nSmoothing spline coded profile using window length of %d.\n", l_med_filt_win);

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

	fprintf(stderr, "\nFinished, saving spline coded profile.\n");
	char encoded_profile_bin_fp[1000];
	sprintf(encoded_profile_bin_fp, "%s/spline_coded_%s.bin", signal_dir, chr_id);
	dump_per_nucleotide_binary_profile(spline_fit_profile, l_track, encoded_profile_bin_fp);

	fprintf(stderr, "Compressing the profile.\n");
	char comp_encoded_profile_bin_fp[1000];
	sprintf(comp_encoded_profile_bin_fp, "%s/spline_coded_%s.bin.gz", signal_dir, chr_id);
	compressFile(encoded_profile_bin_fp, comp_encoded_profile_bin_fp);

	// Dump the bedgraph.
	//for (int i = 1; i < l_track; i++)
	//{
	//	spline_fit_profile[i] *= 10;
	//}
	//floorize_profile(spline_fit_profile, l_track);
	char spline_bgr_fp[1000];
	sprintf(spline_bgr_fp, "%s/spline_coded_%s.bgr", signal_dir, chr_id);
	dump_bedGraph_per_per_nucleotide_binary_profile(spline_fit_profile, l_track, chr_id, spline_bgr_fp);
	char comp_spline_bgr_fp[1000];
	sprintf(comp_spline_bgr_fp, "%s/spline_coded_%s.bgr.gz", signal_dir, chr_id);
	compressFile(spline_bgr_fp, comp_spline_bgr_fp);

	// Delete the uncompressed file.
#ifdef _WIN32
	_unlink(encoded_profile_bin_fp);
	_unlink(spline_bgr_fp);
#endif

#ifdef __unix__
	unlink(encoded_profile_bin_fp);
	unlink(spline_bgr_fp);
#endif

	delete[] coeff;
	delete[] signal_profile;
	delete[] spline_fit_profile;
}
