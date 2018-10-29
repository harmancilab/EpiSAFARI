#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "epsfr_episafari_utils.h"
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

#ifdef __unix__
	#include <unistd.h> 
#endif

bool __DUMP_EPISAFARI_UTILS_MESSAGES__ = false;

#define MAX(x, y) (((x)>(y))?(x):(y))
#define MIN(x, y) (((x)<(y))?(x):(y))

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

vector<t_annot_region*>* get_significant_extrema_per_signal_profile(const char* op_dir, const char* chr_id, double* signal_profile, int l_profile,
																	double* multimapp_signal_profile, int l_multimapp_profile,
																	t_extrema_statistic_defition* extrema_statistic_defn)
{
	// Normalize signal?

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

	fprintf(stderr, "Pooled %d extrema, processing the valleys with max trough2summit distance of %d bps.\n", (int)all_extrema_regs->size(), extrema_statistic_defn->max_summit2trough_dist_in_bp);

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
						if (all_extrema_regs->at(i_ext)->height_at_extrema <= extrema_statistic_defn->max_signal_at_trough &&
							all_extrema_regs->at(j_ext)->height_at_extrema > extrema_statistic_defn->min_signal_at_summit &&
							all_extrema_regs->at(j_ext)->height_at_extrema / all_extrema_regs->at(i_ext)->height_at_extrema > extrema_statistic_defn->min_summit2trough_ratio_per_trough)
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

						int update = (dir == 0) ? (-1) : (1);

						for (int i = cur_minima->extrema_posn; 
							i > cur_left_maxima->extrema_posn && i < cur_right_maxima->extrema_posn; 
							i += update)
						{
							double cur_diff = signal_profile[i] - signal_profile[i - 1];
							if (dir == 0)
							{
								cur_diff *= -1;
							}

							if(cur_diff >= 0)
							{
								cur_dir_total_pos_deriv_signs++;
							}

							if(cur_diff < 0)
							{
								cur_dir_total_neg_deriv_signs++;
							}
						} // i loop.

						left_right_pos_der_frac[dir] = cur_dir_total_pos_deriv_signs / (cur_dir_total_pos_deriv_signs + cur_dir_total_neg_deriv_signs);
					} // dir loop.

					// Compute the average multi-map signal within the valley.
					double total_multimapp_signal = 0;
					double max_multimapp_signal = 10000;
					for (int i = cur_left_maxima->extrema_posn; i < cur_right_maxima->extrema_posn; i++)
					{
						if (i < l_multimapp_profile &&
							max_multimapp_signal > multimapp_signal_profile[i])
						{
							max_multimapp_signal = multimapp_signal_profile[i];

							total_multimapp_signal += multimapp_signal_profile[i];
						}
					} // i loop.

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
					void** valley_info = new void*[3];
					valley_info[0] = valley_points;

					// Set the multimapp signal levels.
					double* cur_valley_multimapp_signals = new double[5];
					cur_valley_multimapp_signals[0] = total_multimapp_signal;
					cur_valley_multimapp_signals[1] = max_multimapp_signal;
					valley_info[1] = cur_valley_multimapp_signals;

					// Set the positive derivative fractions.
					valley_info[2] = left_right_pos_der_frac;

					new_valley->data = valley_info;

					significant_valleys->push_back(new_valley);
				} // right_max_i loop.
			} // left_max_i loop.
		} // minima check.
	} // i_ext loop.

	// Dump the valleys file.
	char valleys_bed_fp[1000];
	sprintf(valleys_bed_fp, "%s/significant_valleys_%s.bed", op_dir, chr_id);
	FILE* f_valleys = open_f(valleys_bed_fp, "w");
	for (int i_val = 0; i_val < (int)(significant_valleys->size()); i_val++)
	{
		void** valley_info = (void**)(significant_valleys->at(i_val)->data);

		t_extrema_node** valley_points = (t_extrema_node**)(valley_info[0]);
		double* cur_valley_multimapp_signals = (double*)(valley_info[1]);
		double* left_right_pos_der_frac = (double*)(valley_info[2]);

		fprintf(f_valleys, "%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", 
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
				left_right_pos_der_frac[1]);
	} // i_val loop.
	fclose(f_valleys);

	// Free minima maxima memory.

	delete[] derivative_map;
	delete all_extrema_regs;
	delete minima_nodes;
	delete maxima_nodes;
	return(significant_valleys);
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

void annotate_features(char* signal_directory, char* gff_fp, int l_half_prom, char* op_fp)
{
	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", signal_directory);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	if (chr_ids == NULL)
	{
		fprintf(stderr, "Could not load chromosome id's from %s.\n", chr_ids_fp);
		exit(0);
	}

	// Load the gff file.
	vector<t_annot_region*>* annotations = load_annotation(gff_fp, l_half_prom);
	fprintf(stderr, "Loaded %d annotation elements.\n", (int)annotations->size());

	FILE* f_op = open_f(op_fp, "w");
	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Annotating features on %s\n", chr_ids->at(i_chr));

		// Load the feats.
		char valleys_bed_fp[1000];
		sprintf(valleys_bed_fp, "%s/significant_valleys_%s.bed", signal_directory, chr_ids->at(i_chr));

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

		// Dump annotations.				
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
				// Process all the annotations.
				for (int i_ann = 0; i_ann < (int)annots->size(); i_ann++)
				{
					fprintf(f_op, "%s:%s;", annots->at(i_ann)->element_type, annots->at(i_ann)->element_name);
				} // i_ann loop.
			}

			fprintf(f_op, "\n");
		} // i_val loop.
	} // i_chr loop.
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
				total_err = 0;
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
			fprintf(stderr, "No fit @ %s:%d-%d\n", chr_id, start_i, end_i);

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
