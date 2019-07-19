#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "epsfr_ansi_string.h"
#include "epsfr_nomenclature.h"
#include "epsfr_annot_region_tools.h"
#include "epsfr_signal_track_tools.h"
#include "epsfr_utils.h"
#include "epsfr_mapped_read_tools.h"
#include "epsfr_gsl_polyfit_utils.h"
#include "epsfr_episafari_utils.h"
#include "epsfr_min_max_utils.h"
#include "epsfr_rng.h"
#include "epsfr_seed_manager.h"
#include "epsfr_genomics_coords.h"
#include "epsfr_ansi_cli.h"

#include <time.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>

void print_usage(char* argv[])
{
	fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Options:\n\
Input Data Processing:\n\
	Mapped Read Processing:\n\
		-preprocess_reads [File format (\"SAM\"/\"ELAND\"/\"bowtie\"/\"tagAlign\"/\"BED5\"/\"BED6\")] [Mapped reads file path (\"stdin\" for piped input)] [Output directory]\n\
		-sort_reads [Preprocessed read directory] [Sorted reads output directory]\n\
		-remove_duplicates [Sorted preprocessed read directory] [Max # of duplicates per position] [Output directory]\n\
	BedGraph/Sequence Processing:\n\
		-separate_bedGraph_2_chromosomes [bedGraph file path (\"stdin\" for piped input)] [Output directory]\n\
		-preprocess_FASTA [Directory with fasta files] [Sequence extension (e.g., fasta, fa)] [Output directory]\n\
	B-spline Processing:\n\
		-save_basis_splines [Number of breakpoints including ends] [Spline Degree]\n\
Valley Detection:\n\
	-bspline_encode [Options] [Values]\n\
	-get_valleys [Options] [Values]\n\
	-assign_valleys_2_regions [Regions BED file path] [Valleys BED file path]\n\
Differential Valley Analysis:\n\
	-get_2_sample_differential_valleys [Sample 1 valleys BED file path] [Sample 1 directory] [Sample 2 valleys BED file path] [Sample 2 directory] [Sparse data flag (0/1)] [P-value type: 0: Binomial (Mult.), 1: Binomial (Merge), 2: Multinomial]\n\
Valley Annotation:\n\
	-annotate_features [Valleys BED file path] [GFF file path] [Half promoter length] [Output file path]\n\
Signal Matrix Processing:\n\
	-append_signal_2_regions [Signals directory] [Fragment length] [BED file path] [Output file path]\n\
	-quantile_normalize_matrix [Signal matrix file path] [Output file path]\n\
Dip Normalization:\n\
	-uniformize_dip_relative_window_positions [Valleys BED file path] [Window length] [# bins] [Output file path]\n\
Profile Simulation:\n\
	-random_copy_profile_valleys [Signal file path] [Regions BED file path] [Random profile length] [Output profile file path]\n", argv[0]);
}

int main(int argc, char* argv[])
{
	clock_t start_c = clock();

	if (argc < 2)
	{
		print_usage(argv);
		exit(0);
	}

	if (t_string::compare_strings(argv[1], "-random_copy_profile_valleys"))
	{
		if (argc < 6)
		{
			fprintf(stderr, "USAGE: %s -random_copy_profile_valleys [Signal file path] [Regions BED file path] [Random profile length] [Output signal file path]\n", argv[0]);
			exit(0);
		}

		char* signal_fp = argv[2];
		char* regs_BED_fp = argv[3];
		int l_sim_profile = atoi(argv[4]);
		char* op_fp = argv[5];

		int l_profile = 0;
		int l_frag = 200;
		bool reads_loaded = false;
		double* signal_profile = load_signal_covg_per_signal_file(signal_fp,
																l_frag,
																l_profile,
																reads_loaded);
		fprintf(stderr, "Loaded %d nucleotide profile.\n", l_profile);

		// Load the valley regions.
		vector<t_annot_region*>* valley_regs = load_BED_with_line_information(regs_BED_fp);
		fprintf(stderr, "Assigning dips to valley regions.\n");
		for (int i_reg = 0; i_reg < (int)valley_regs->size(); i_reg++)
		{
			char* cur_reg_line = (char*)(valley_regs->at(i_reg)->data);
			int cur_valley_dip_posn = 0;
			if (sscanf(cur_reg_line, "%*s %*s %*s %d", &cur_valley_dip_posn) != 1)
			{
				fprintf(stderr, "Could not read the dip position from %s\n", cur_reg_line);
				exit(0);
			}

			valley_regs->at(i_reg)->score = cur_valley_dip_posn;
		} // i_reg loop.
		
		fprintf(stderr, "Generating simulated profile of length %d base pairs using %d valley regions.\n", l_sim_profile, (int)valley_regs->size());

		double* simulated_profile = new double[l_sim_profile + 10];
		memset(simulated_profile, 0, sizeof(double) * (l_sim_profile + 2));

		t_rng* rng = new t_rng(t_seed_manager::seed_me());
		fprintf(stderr, "Shuffling positions.\n");
		vector<int>* sim_profile_reg_posn_indices = rng->fast_permute_indices(1, l_sim_profile);
		fprintf(stderr, "Shuffled positions:\n");
		for (int i = 0; i < 10; i++)
		{
			fprintf(stderr, "%d ", sim_profile_reg_posn_indices->at(i));
		} // i loop.
		fprintf(stderr, "...\n");

		FILE* f_simulated_dips = open_f("simulated_dips.bed", "w");
		int simulated_posn_i = 0;
		for (int i_reg = 0; i_reg < (int)valley_regs->size(); i_reg++)
		{
			fprintf(stderr, "Adding %d. region.           \r", i_reg);

			// Copy the signal from profile to simulated profile.
			bool position_available = false;

			while (!position_available)
			{
				int cur_reg_sim_start = sim_profile_reg_posn_indices->at(simulated_posn_i);
				int cur_reg_sim_end = sim_profile_reg_posn_indices->at(simulated_posn_i) + valley_regs->at(i_reg)->end - valley_regs->at(i_reg)->start + 1;

				// Save the dip's location.
				int cur_valley_sim_dip = valley_regs->at(i_reg)->score - valley_regs->at(i_reg)->start + cur_reg_sim_start;

				// Check if this position is available.
				position_available = true;
				for (int i = cur_reg_sim_start; i < cur_reg_sim_end; i++)
				{
					if (simulated_profile[i] > 0)
					{
						position_available = false;
						break;
					}
				} // i loop.

				if (position_available)
				{
					fprintf(f_simulated_dips, "1\t%d\t%d\t%s\n", cur_valley_sim_dip - 100, cur_valley_sim_dip + 100, (char*)(valley_regs->at(i_reg)->data));

					for (int i = cur_reg_sim_start; i < cur_reg_sim_end; i++)
					{
						simulated_profile[i] = signal_profile[i - cur_reg_sim_start + valley_regs->at(i_reg)->start];
					} // i loop.
					
					break;
				}

				// Push the simulated position forward and move to next region/same region.
				simulated_posn_i++;
			} // position_available check.
		} // i loop.
		fclose(f_simulated_dips);

		// Save.
		dump_per_nucleotide_binary_profile(simulated_profile, l_sim_profile, op_fp);
	} // -random_copy_profile_regions
	else if (t_string::compare_strings(argv[1], "-uniformize_dip_relative_window_positions"))
	{
		if (argc < 5)
		{
			fprintf(stderr, "usage: %s -uniformize_dip_relative_window_positions [Valleys BED file path] [Window length] [# bins] [Output file path]\n", argv[0]);
		}

		char* valleys_BED_fp = argv[2];
		int l_win = atoi(argv[3]);
		int n_bins_per_win = atoi(argv[4]);
		char* op_fp = argv[5];

		uniformize_dip_relative_window_positions(valleys_BED_fp, l_win, n_bins_per_win, op_fp);
	} // -uniformize_dip_relative_window_positions option.

	if (strcmp(argv[1], "-help") == 0 ||
		strcmp(argv[1], "-version") == 0 ||
		strcmp(argv[1], "-h") == 0 ||
		strcmp(argv[1], "-v") == 0)
	{
		print_usage(argv);
		exit(0);
	}
	else if (strcmp(argv[1], "-quantile_normalize_matrix") == 0)
	{
		if (argc != 4)
		{
			fprintf(stderr, "%s -quantile_normalize_matrix [Signal matrix file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* signal_regions_BED_fp = argv[2];
		char* op_fp = argv[3];
		quantile_normalize_signal_matrix(signal_regions_BED_fp, op_fp);
	} // -quantile_normalize_matrix option.
	else if (strcmp(argv[1], "-get_2_sample_differential_valleys") == 0)
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s -get_2_sample_differential_valleys [Sample 1 valleys BED file path] [Sample 1 directory] [Sample 2 valleys BED file path] [Sample 2 directory] [Sparse data flag (0/1)] [P-value type: 0: Binomial (Mult.), 1: Binomial (Merge), 2: Multinomial]\n", argv[0]);
			exit(0);
		}

		char* sample1_valleys_bed_fp = argv[2];
		char* sample1_dir = argv[3];
		char* sample2_valleys_bed_fp = argv[4];
		char* sample2_dir = argv[5];
		bool sparse_valleys = (argv[6][0] == '1');
		char p_val_type = atoi(argv[7]);

		// set the extrema statistics to be used in computing significance.
		t_extrema_statistic_defition* extrema_statistic_defn = new t_extrema_statistic_defition();
		extrema_statistic_defn->max_signal_at_trough = -1;
		extrema_statistic_defn->min_signal_at_summit = -1;
		extrema_statistic_defn->min_summit2trough_ratio_per_trough = -1;
		extrema_statistic_defn->min_summit2trough_dist_in_bp = -1;
		extrema_statistic_defn->max_summit2trough_dist_in_bp = -1;
		extrema_statistic_defn->max_multimapp_signal_at_trough = -1;
		extrema_statistic_defn->log_q_val_threshold = -1;
		extrema_statistic_defn->p_val_estimate_extrema_vic_window_length = 50;
		extrema_statistic_defn->p_val_estimate_signal_scaling_factor = 1.0;
		extrema_statistic_defn->sparse_profile = sparse_valleys;
		extrema_statistic_defn->l_minima_vicinity_per_merging = 100;
		extrema_statistic_defn->p_val_type = p_val_type;

		//extrema_statistic_defn->hill_score_type = HEIGHT_BASED_HILL_SCORE;
		//extrema_statistic_defn->hill_score_type = DIST_BASED_HILL_SCORE;
		extrema_statistic_defn->hill_score_type = -1;

		if (sparse_valleys)
		{
			extrema_statistic_defn->p_val_estimate_signal_scaling_factor = 100;
		}

		get_2_sample_differential_valleys(sample1_valleys_bed_fp, sample1_dir, sample2_valleys_bed_fp, sample2_dir, extrema_statistic_defn);
	} // -get_2_sample_differential_valleys option.
	else if (strcmp(argv[1], "-save_basis_splines") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s -save_basis_splines [Number of breakpoints including ends] [Spline Degree] [Breakpoint type (0: Uniform / 1: Derivative / 2: Vicinity Derivative / 3: Random)]\n", argv[0]);
			exit(0);
		}

		int n_brkpts = atoi(argv[2]);
		int spline_degree = atoi(argv[3]);
		int brkpts_type_index = atoi(argv[4]);

		char brkpts_type = -1;
		if (brkpts_type_index == 0)
		{
			fprintf(stderr, "Uniform breakpoints.\n");
			brkpts_type = UNIFORM_BREAKPOINTS;
		}
		else if (brkpts_type_index == 1)
		{
			fprintf(stderr, "Cannot use Derivative-based breakpoints.\n");
			brkpts_type = DERIVATIVE_NU_BREAKPOINTS;
			exit(0);
		}
		else if (brkpts_type_index == 2)
		{
			fprintf(stderr, "Cannot use Vicinity derivative-based breakpoints.\n");
			brkpts_type = VICINITY_DERIVATIVE_NU_BREAKPOINTS;
			exit(0);
		}
		else if (brkpts_type_index == 3)
		{
			fprintf(stderr, "Random breakpoints.\n");
			brkpts_type = RANDOM_NU_BREAKPOINTS;
		}

		int n_spline_coeff = spline_degree + n_brkpts - 2;

		int max_n_internal_breakpoints = n_brkpts - 2;

		int n_data_pts = 1000;
		double* dx = new double[n_data_pts + 2];
		double* dy = new double[n_data_pts + 2];

		for (int i = 0; i < n_data_pts; i++)
		{
			dx[i] = i;
			dy[i] = i^2;
		} // i loop.

		double* reconst_y = new double[n_data_pts + 2];

		double* coeff = new double[n_spline_coeff + 2];

		fprintf(stderr, "Saving the basis splines for %d data points using %d breakpoints and %d-degree spline.\n", n_data_pts, n_brkpts, spline_degree);

		t_rng* rng = new t_rng(t_seed_manager::seed_me());

		// Save the basis splines using uniform breakpoints.
		bsplinefit_nonuniform_per_breakpoint_type(n_data_pts,
			max_n_internal_breakpoints,
			brkpts_type,
			spline_degree,
			dx, dy,
			reconst_y,
			coeff,
			rng,
			true);
	} // -save_basis_splines option.
	else if (strcmp(argv[1], "-assign_valleys_2_regions") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s -assign_valleys_2_regions [Regions BED file path] [Valleys BED file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* regs_BED_fp = argv[2];
		char* valleys_BED_fp = argv[3];
		char* op_fp = argv[4];

		vector<t_annot_region*>* regs = load_BED_with_line_information(regs_BED_fp);
		for (int i_reg = 0; i_reg < (int)regs->size(); i_reg++)
		{
			void* cur_dat = regs->at(i_reg)->data;
			void** new_dat = new void*[2];
			new_dat[0] = cur_dat;
			new_dat[1] = new vector<t_annot_region*>();
			regs->at(i_reg)->data = new_dat;
		} // i_reg loop.

		vector<t_annot_region*>* valleys_regs = load_BED_with_line_information(valleys_BED_fp);
		fprintf(stderr, "Assigning %d valleys to %d regions.\n", (int)valleys_regs->size(), (int)regs->size());

		vector<t_annot_region*>* intersects = intersect_annot_regions(regs, valleys_regs, false, true);
		for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
		{
			t_intersect_info* cur_int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* cur_reg = cur_int_info->src_reg;
			t_annot_region* cur_valley_reg = cur_int_info->dest_reg;

			vector<t_annot_region*>* cur_reg_valleys = (vector<t_annot_region*>*)(((void**)(cur_reg->data))[1]);
			cur_reg_valleys->push_back(cur_valley_reg);

			delete cur_int_info;
		} // i_int loop.
		delete_annot_regions(intersects);

		// Dump the regions with valleys.
		FILE* f_op = open_f(op_fp, "w");
		for (int i_reg = 0; i_reg < (int)regs->size(); i_reg++)
		{
			char* cur_reg_line = (char*)(((void**)(regs->at(i_reg)->data))[0]);
			vector<t_annot_region*>* cur_reg_valleys = (vector<t_annot_region*>*)(((void**)(regs->at(i_reg)->data))[1]);

			fprintf(f_op, "%s\t%d", cur_reg_line, (int)cur_reg_valleys->size());

			if ((int)cur_reg_valleys->size() > 0)
			{
				sort(cur_reg_valleys->begin(), cur_reg_valleys->end(), sort_regions);

				fprintf(f_op, "\t%d-%d", cur_reg_valleys->at(0)->start, cur_reg_valleys->back()->end);
				for (int i_val = 0; i_val < (int)cur_reg_valleys->size(); i_val++)
				{
					fprintf(f_op, ";%d-%d", cur_reg_valleys->at(i_val)->start, cur_reg_valleys->at(i_val)->end);
				} // i_val loop.
			}
			else
			{
				fprintf(f_op, "\t%d-%d", regs->at(i_reg)->start, regs->at(i_reg)->end);
			}

			fprintf(f_op, "\n");
		} // i_reg loop.
		fclose(f_op);
	} // -assign_valleys_2_regions option.
	else if (strcmp(argv[1], "-preprocess_FASTA") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s -preprocess_FASTA [Directory with fasta files] [Sequence extension (e.g., fasta, fa)] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* fasta_dir = argv[2];
		char* seq_ext = argv[3];
		char* op_dir = argv[4];

		vector<char*>* fasta_fns = load_directory_files(fasta_dir, seq_ext);
		printf("Preprocessing %d sequence FASTA files in %s\n", (int)fasta_fns->size(), fasta_dir);

		// Binarize each file separately.
		for (int i_fa = 0; i_fa < (int)fasta_fns->size(); i_fa++)
		{
			//binarize_genome_per_fasta(fasta_fps->at(i_fa), op_dir);
			char cur_chr_fa_fp[1000];
			sprintf(cur_chr_fa_fp, "%s/%s", fasta_dir, fasta_fns->at(i_fa));
			binarize_fasta_file(cur_chr_fa_fp, op_dir);
		} // i_fa loop.
	}
	else if (strcmp(argv[1], "-annotate_features") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s -annotate_features [Valleys BED file path] [GFF file path] [Half promoter length] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* valleys_bed_fp = argv[2];
		char* gff_fp = argv[3];
		int l_prom = atoi(argv[4]);
		char* op_fp = argv[5];

		annotate_features(valleys_bed_fp, gff_fp, l_prom/2, op_fp);
	} // -annotate_features option.
	else if (strcmp(argv[1], "-sort_reads") == 0)
	{
		if (argc != 4)
		{
			fprintf(stderr, "%s -sort_reads [Preprocessed read directory] [Sorted reads output directory]\n", argv[0]);
			exit(0);
		}

		char* preprocessed_reads_dir = argv[2];
		char* sorted_reads_op_dir = argv[3];

		int bucket_size = 50 * 1000 * 1000;

		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dir);

		// Process per chromosome.
		vector<char*>* chr_ids = buffer_file(chr_ids_fp);
		if (chr_ids == NULL)
		{
			fprintf(stderr, "Could not fine the chromosome id's list file @ %s\n", chr_ids_fp);
			exit(0);
		}

		char sorted_chr_ids_fp[1000];
		sprintf(sorted_chr_ids_fp, "%s/chr_ids.txt", sorted_reads_op_dir);
		FILE* f_sorted_chr_ids = open_f(sorted_chr_ids_fp, "w");
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(f_sorted_chr_ids, "%s\n", chr_ids->at(i_chr));
		} // i_chr loop.
		fclose(f_sorted_chr_ids);

		// After sorting, dump the reads into a new directory named xxxx_sorted/.
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Processing %s\n", chr_ids->at(i_chr));
			char cur_chr_reads_fp[1000];
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt.gz", preprocessed_reads_dir, chr_ids->at(i_chr));

			FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");

			// Load the entries, store them in a buffer.
			vector<FILE*>* bucket_f_list = new vector<FILE*>();
			vector<int>* bucket_starts = new vector<int>();
			vector<int>* bucket_sizes = new vector<int>();
			while (1)
			{
				// Get the read and determine which bucket it belongs to.
				char* cur_line = getline(f_cur_chr_reads);
				if (cur_line == NULL)
				{
					break;
				}

				char cur_cigar_str[1000];
				char cur_strand;
				int cur_start;
				sscanf(cur_line, "%s %c %d", cur_cigar_str, &cur_strand, &cur_start);
				int bucket_start = (int)(floor((double)cur_start / bucket_size) * bucket_size);

				bool bucket_found = false;
				for (int i_st = 0; i_st < (int)bucket_starts->size(); i_st++)
				{
					if (bucket_start == bucket_starts->at(i_st))
					{
						bucket_found = true;
						fprintf(bucket_f_list->at(i_st), "%s\n", cur_line);
						bucket_sizes->at(i_st)++;
					}
				} // i_st loop.

				  // Need a new bucket?
				if (bucket_found == false)
				{
					bucket_starts->push_back(bucket_start);
					char new_bucket_fp[1000];
					//sprintf(new_bucket_fp, "%s/bucket_%d.txt", preprocessed_reads_dir, bucket_start);
					sprintf(new_bucket_fp, "%s/bucket_%d.txt", sorted_reads_op_dir, bucket_start);

					//if (__DUMP_PEAK_MESSAGES__)
					//	fprintf(stderr, "Adding bucket %s\n", new_bucket_fp);

					FILE* f_bucket = open_f(new_bucket_fp, "w");
					bucket_f_list->push_back(f_bucket);
					bucket_sizes->push_back(0);

					bucket_found = false;
					for (int i_st = 0; i_st < (int)bucket_starts->size(); i_st++)
					{
						if (bucket_start == bucket_starts->at(i_st))
						{
							bucket_found = true;
							fprintf(bucket_f_list->at(i_st), "%s\n", cur_line);
							bucket_sizes->at(i_st)++;
						}
					} // i_st loop.				
				}

				delete[] cur_line;
			} // file reading loop.

			close_f(f_cur_chr_reads, cur_chr_reads_fp);

			for (int i_st = 0; i_st < (int)bucket_f_list->size(); i_st++)
			{
				fclose(bucket_f_list->at(i_st));
			} // i_st loop.
			delete(bucket_f_list);

			char cur_chr_sorted_reads_fp[1000];
			sprintf(cur_chr_sorted_reads_fp, "%s/%s_mapped_reads.txt.gz", sorted_reads_op_dir, chr_ids->at(i_chr));
			FILE* f_cur_chr_sorted_reads = open_f(cur_chr_sorted_reads_fp, "w");
			sort(bucket_starts->begin(), bucket_starts->end());
			for (int i_buck = 0; i_buck < (int)bucket_starts->size(); i_buck++)
			{
				// Load the reads in the current bucket, sort them, dump them.
				char cur_bucket_fp[1000];
				//sprintf(cur_bucket_fp, "%s/bucket_%d.txt", preprocessed_reads_dir, bucket_starts->at(i_buck));
				sprintf(cur_bucket_fp, "%s/bucket_%d.txt", sorted_reads_op_dir, bucket_starts->at(i_buck));

				//if (__DUMP_PEAK_MESSAGES__)
				//	fprintf(stderr, "Sorting reads in bucket %s.\n", cur_bucket_fp);

				vector<char*>* cur_sorted_bucket_read_lines = sort_bucket_read_lines(cur_bucket_fp);

				for (int i_l = 0; i_l < (int)cur_sorted_bucket_read_lines->size(); i_l++)
				{
					fprintf(f_cur_chr_sorted_reads, "%s\n", cur_sorted_bucket_read_lines->at(i_l));

					// Free memory.	
					delete[] cur_sorted_bucket_read_lines->at(i_l);
				} // i_l loop.

				delete cur_sorted_bucket_read_lines;
			} // i_buck loop.
			close_f(f_cur_chr_sorted_reads, cur_chr_sorted_reads_fp);

			delete bucket_starts;
			delete bucket_sizes;
		} // i_chr loop.
	} // -sort_reads
	else if (strcmp(argv[1], "-remove_duplicates") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -remove_duplicates [Sorted preprocessed read directory] [Max # of duplicates per position] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* sorted_preprocessed_reads_dir = argv[2];
		int max_n_amp_reads = atoi(argv[3]);
		char* pruned_reads_dir = argv[4];


		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", sorted_preprocessed_reads_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_fp);
		if (chr_ids == NULL)
		{
			fprintf(stderr, "Could not fine the chromosome id's list file @ %s\n", chr_ids_fp);
			exit(0);
		}

		char pruned_chr_ids_fp[1000];
		sprintf(pruned_chr_ids_fp, "%s/chr_ids.txt", pruned_reads_dir);
		FILE* f_pruned_chr_ids = open_f(pruned_chr_ids_fp, "w");
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(f_pruned_chr_ids, "%s\n", chr_ids->at(i_chr));
		} // i_chr loop.
		fclose(f_pruned_chr_ids);

		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Pruning %s\n", chr_ids->at(i_chr));
			char cur_chr_reads_fp[1000];
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt.gz", sorted_preprocessed_reads_dir, chr_ids->at(i_chr));
			FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");

			char cur_chr_pruned_reads_fp[1000];
			sprintf(cur_chr_pruned_reads_fp, "%s/%s_mapped_reads.txt.gz", pruned_reads_dir, chr_ids->at(i_chr));
			FILE* f_cur_chr_pruned_reads = open_f(cur_chr_pruned_reads_fp, "w");

			int n_processed_reads = 0;
			int n_pruned_reads = 0;
			int prev_read_start = 0;
			int n_amp_reads = 0;
			while (1)
			{
				char* cur_line = getline(f_cur_chr_reads);
				if (cur_line == NULL)
				{
					break;
				}

				//if (__DUMP_PEAK_MESSAGES__)
				//{
				//	if (n_processed_reads % 1000000 == 0)
				//	{
				//		fprintf(stderr, "Processing %d. read.              \r", n_processed_reads);
				//	}
				//}

				int cur_read_start = 0;
				if (sscanf(cur_line, "%*s %*s %d", &cur_read_start) != 1)
				{
					fprintf(stderr, "Could not parse: %s\n", cur_line);
					exit(0);
				}

				// Check if the read start updated.
				if (cur_read_start == prev_read_start)
				{
					n_amp_reads++;
				}
				else
				{
					prev_read_start = cur_read_start;
					n_amp_reads = 1;
				}

				if (n_amp_reads <= max_n_amp_reads)
				{
					n_pruned_reads++;
					fprintf(f_cur_chr_pruned_reads, "%s\n", cur_line);
				}

				n_processed_reads++;

				// Free line memory.
				delete[] cur_line;
			} // file reading loop.

			close_f(f_cur_chr_reads, cur_chr_reads_fp);
			close_f(f_cur_chr_pruned_reads, cur_chr_pruned_reads_fp);

			fprintf(stderr, "Processed %d reads, pruned to %d reads.\n", n_processed_reads, n_pruned_reads);
		} // i_chr loop.
	} // -remove_duplicates
	else if (strcmp(argv[1], "-separate_bedGraph_2_chromosomes") == 0)
	{
		if (argc != 4)
		{
			fprintf(stderr, "USAGE: %s -separate_bedGraph_2_chromosomes [bedGraph file path (\"stdin\" for piped input)] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* bedgraph_fp = argv[2];
		char* op_dir = argv[3];

		vector<char*>* chr_ids = new vector<char*>();
		vector<FILE*>* f_bgr_op_list = new vector<FILE*>();
		vector<char*>* bgr_op_list_fps = new vector<char*>();

		fprintf(stderr, "Separating %s with respect to chromosomes and saving to %s.\n", bedgraph_fp, op_dir);

		FILE* f_bgr = open_f(bedgraph_fp, "r");

		while (1)
		{
			char* cur_line = getline(f_bgr);
			if (cur_line == NULL)
			{
				break;
			}

			// Add the chromosome if it does not exist.
			char cur_chr_id[1000];
			sscanf(cur_line, "%s", cur_chr_id);
			char* copy_chr_id = t_string::copy_me_str(cur_chr_id);
			normalize_chr_id(copy_chr_id);

			int cur_chr_i = t_string::get_i_str(chr_ids, copy_chr_id);
			if (cur_chr_i == (int)chr_ids->size())
			{				
				chr_ids->push_back(t_string::copy_me_str(copy_chr_id));
				
				fprintf(stderr, "Adding %s\n", copy_chr_id);

				char bgr_op_fp[1000];
				sprintf(bgr_op_fp, "%s/%s.bgr.gz", op_dir, copy_chr_id);

				if (check_file(bgr_op_fp))
				{
					fprintf(stderr, "%s exists, pooling on it.\n", bgr_op_fp);
				}

				FILE* cur_f_bgr_op = open_f(bgr_op_fp, "a");
				f_bgr_op_list->push_back(cur_f_bgr_op);
				bgr_op_list_fps->push_back(bgr_op_fp);
			} // new chromosome check.

			// Write the current line.
			cur_chr_i = t_string::get_i_str(chr_ids, copy_chr_id);
			if (cur_chr_i == (int)chr_ids->size())
			{
				fprintf(stderr, "Sanity check failed, failed to add chromosome: %s\n", cur_line);
				exit(0);
			}

			fprintf(f_bgr_op_list->at(cur_chr_i), "%s\n", cur_line);

			delete[] copy_chr_id;
			delete[] cur_line;
		} // file reading loop.

		// Close file.
		close_f(f_bgr, bedgraph_fp);

		if (f_bgr_op_list->size() != chr_ids->size())
		{
			fprintf(stderr, "Sanity check failed, chromsome list is not same size as file list.\n");
			exit(0);
		}

		// Write the chromosome ids and close the bedgraph files.
		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.txt", op_dir);
		FILE* f_chr_ids_list = open_f(chr_ids_list_fp, "w");
		for (int i_chr = 0; i_chr < (int)f_bgr_op_list->size(); i_chr++)
		{
			fprintf(f_chr_ids_list, "%s\n", chr_ids->at(i_chr));
			close_f(f_bgr_op_list->at(i_chr), bgr_op_list_fps->at(i_chr));
		} // i_chr loop.
		fclose(f_chr_ids_list);
	} // -separate_bedGraph_2_chromosomes option.
	else if (strcmp(argv[1], "-preprocess_reads") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -preprocess [File format (\"SAM\"/\"ELAND\"/\"bowtie\"/\"tagAlign\"/\"BED5\"/\"BED6\")] [Mapped reads file path (\"stdin\" for piped input)] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* format_str = argv[2];
		char* mapped_reads_fp = argv[3];
		char* op_dir = argv[4];

		printf("Preprocessing:\nInput Format: %s\nChIP-Seq input file path: %s\nParsed Reads Output Directory: %s\n",
			format_str, mapped_reads_fp, op_dir);

		// Do not do validation while dumping the binary read files.
		if (strcmp(format_str, "ELAND") == 0)
		{
			// Read the ELAND file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_ELAND_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_ELAND_read_line, false);
		}
		else if (strcmp(format_str, "SAM") == 0)
		{
			// Read the SAM file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_SAM_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_SAM_read_line, false);
		}
		else if (strcmp(format_str, "tagAlign") == 0)
		{
			// Read the SAM file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_tagAlign_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_tagAlign_read_line, false);
		}
		else if (strcmp(format_str, "bowtie") == 0)
		{
		// Read the bowtie file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
		//parse_bowtie_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
		preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_bowtie_read_line, false);
		}
		else if (strcmp(format_str, "BED5") == 0)
		{
			// Read the bowtie file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_bowtie_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_BED5_read_line, false);
		}
		else if (strcmp(format_str, "BED6") == 0)
		{
			// Read the bowtie file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_bowtie_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_BED6_read_line, false);
		}
		else
		{
			printf("Unknown format for the mapped read file name, use SAM/ELAND/bowtie/tagAlign.\n");
			exit(0);
		}
	} // -preprocess option.
	if (t_string::compare_strings(argv[1], "-get_valleys"))
	{
		if (argc < 3)
		{
			fprintf(stderr, "USAGE: %s -get_valleys [Options] [Values]\n\
	-signal_dir [Encoded signals data directory path]\n\
	-mmap_dir [Multi-mappability profile directory]\n\
	-genome_dir [Genome sequence directory]\n\
	-max_signal_at_trough [Maximum signal at trough]\n\
	-min_signal_at_summit [Minimum signal at summit]\n\
	-f_min [Minimum summit2trough ratio per trough]\n\
	-l_min [Minimum summit2trough distance in bps]\n\
	-l_max [Maximum summit2trough distance in bps]\n\
	-max_mmap [Maximum multimapp signal at trough]\n\
	-max_qval [Q-val threshold]\n\
	-sparse_profile [Sparse data flag (0/1)]\n\
	-pval_type [P-value type: 0: Binomial (Mult.), 1: Binomial (Merge), 2: Multinomial]\n\
	-l_p [P-value estimation window length]\n", argv[0]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "-");

		bool success = false;

		char* signal_data_dir = cli->get_value_by_option("-signal_dir", success);if (!success) { fprintf(stderr, "Could not read signal data directory.\n");exit(0); }
		char* multi_mapp_signal_profile_dir = cli->get_value_by_option("-mmapp_dir", success);if (!success) { fprintf(stderr, "Could not read multimapp directory.\n");exit(0); }
		char* genome_seq_dir = cli->get_value_by_option("-genome_dir", success);if (!success) { fprintf(stderr, "Could not read genome sequence directory.\n");exit(0); }

		double max_signal_at_trough = atof_null(cli->get_value_by_option("-max_signal_at_trough", success));if (!success) { max_signal_at_trough = 10000; }
		double min_signal_at_summit = atof_null(cli->get_value_by_option("-min_signal_at_summit", success));if (!success) { min_signal_at_summit = 5; }
		double min_summit2trough_ratio_per_trough = atof_null(cli->get_value_by_option("-f_min", success));if (!success) { min_summit2trough_ratio_per_trough = 1.2; }
		int min_summit2trough_dist_in_bp = atoi_null(cli->get_value_by_option("-l_min", success));if (!success) { min_summit2trough_dist_in_bp = 0; }
		int max_summit2trough_dist_in_bp = atoi_null(cli->get_value_by_option("-l_max", success));if (!success) { max_summit2trough_dist_in_bp = 1000; }
		double max_multimapp_signal_at_trough = atof_null(cli->get_value_by_option("-max_mmap", success));if (!success) { max_multimapp_signal_at_trough = 1.2; }
		double log_q_val_threshold = log(atof_null(cli->get_value_by_option("-max_qval", success)));if (!success) { log_q_val_threshold = 0.01; }

		int l_p = atoi_null(cli->get_value_by_option("-l_p", success));if (!success) { l_p = 50; }
		
		bool sparse_valleys = false;
		if (!cli->get_value_by_option("-sparse_profile", success)) { sparse_valleys = false; }
		else { sparse_valleys = (cli->get_value_by_option("-sparse_profile", success)[0] == '1'); }

		char p_val_type = atoi(cli->get_value_by_option("-pval_type", success));if (!success) { p_val_type = 0; }

		fprintf(stderr, "Parameters:\n\
signal_data_dir=%s\n\
max_signal_at_trough=%.3f\n\
min_signal_at_summit=%.3f\n\
min_summit2trough_ratio_per_trough=%.3f\n\
min_summit2trough_dist_in_bp=%d\n\
max_summit2trough_dist_in_bp=%d\n\
multi_mapp_signal_profile_dir=%s\n\
max_multimapp_signal_at_trough=%.3f\n\
genome_seq_dir=%s\n\
log_q_val_threshold=%.3f\n\
sparse_valleys=%d\n\
p_val_type=%d\n\
l_p=%d\n",
			signal_data_dir,
			max_signal_at_trough,
			min_signal_at_summit,
			min_summit2trough_ratio_per_trough,
			min_summit2trough_dist_in_bp,
			max_summit2trough_dist_in_bp,
			multi_mapp_signal_profile_dir,
			max_multimapp_signal_at_trough,
			genome_seq_dir,
			log_q_val_threshold,
			sparse_valleys,
			p_val_type,
			l_p);

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.txt", signal_data_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		if (chr_ids == NULL)
		{
			fprintf(stderr, "Could not load chromosome id's from %s.\n", chr_ids_list_fp);
			exit(0);
		}

		// set the extrema statistics to be used in computing significance.
		t_extrema_statistic_defition* extrema_statistic_defn = new t_extrema_statistic_defition();
		extrema_statistic_defn->max_signal_at_trough = max_signal_at_trough;
		extrema_statistic_defn->min_signal_at_summit = min_signal_at_summit;
		extrema_statistic_defn->min_summit2trough_ratio_per_trough = min_summit2trough_ratio_per_trough;
		extrema_statistic_defn->min_summit2trough_dist_in_bp = min_summit2trough_dist_in_bp;
		extrema_statistic_defn->max_summit2trough_dist_in_bp = max_summit2trough_dist_in_bp;
		extrema_statistic_defn->max_multimapp_signal_at_trough = max_multimapp_signal_at_trough;
		extrema_statistic_defn->log_q_val_threshold = log_q_val_threshold;
		extrema_statistic_defn->p_val_estimate_extrema_vic_window_length = l_p;
		extrema_statistic_defn->p_val_estimate_signal_scaling_factor = 1.0;
		extrema_statistic_defn->sparse_profile = sparse_valleys;
		extrema_statistic_defn->l_minima_vicinity_per_merging = 100;
		extrema_statistic_defn->p_val_type = p_val_type;

		//extrema_statistic_defn->hill_score_type = HEIGHT_BASED_HILL_SCORE;
		extrema_statistic_defn->hill_score_type = DIST_BASED_HILL_SCORE;

		if (sparse_valleys)
		{
			extrema_statistic_defn->p_val_estimate_signal_scaling_factor = 100;
		}

		// Process all the chromosomes.
		vector<t_annot_region*>* all_significant_valleys = new vector<t_annot_region*>();
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			char* chr_id = chr_ids->at(i_chr);

			if (sparse_valleys)
			{
				fprintf(stderr, "Detecting the features on %s (sparse)\n", chr_id);
			}
			else
			{
				fprintf(stderr, "Detecting the features on %s\n", chr_id);
			}

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

			// Load smoothed signal.
			int l_profile = 0;
			double* signal_profile = NULL;
			char bin_signal_fp[1000];
			sprintf(bin_signal_fp, "%s/spline_coded_%s.bin.gz", signal_data_dir, chr_id);
			if (check_file(bin_signal_fp))
			{
				signal_profile = load_per_nucleotide_binary_profile(bin_signal_fp, l_profile);
				fprintf(stderr, "Loaded %d long signal profile.\n", l_profile);
			}
			else
			{
				fprintf(stderr, "Could not load smoothed signal from %s.\n", bin_signal_fp);
				exit(0);
			}						
			
			// Load multimapp signal.
			int l_multimapp_profile = 0;
			unsigned char* multi_mapp_signal = NULL;
			char multi_mapp_signal_profile_fp[1000];
			sprintf(multi_mapp_signal_profile_fp, "%s/%s.bin", multi_mapp_signal_profile_dir, chr_id);
			if (check_file(multi_mapp_signal_profile_fp))
			{		
				multi_mapp_signal = load_normalized_multimappability_profile(multi_mapp_signal_profile_fp, l_multimapp_profile);
				fprintf(stderr, "Loaded %d long multimappability profile.\n", l_multimapp_profile);
			}
			else
			{
				fprintf(stderr, "Could not find multimapp file @ %s, skipping.\n", multi_mapp_signal_profile_fp);
				delete[] signal_profile;
				continue;
			}

			// Load genome sequence.			
			int l_chrom_seq = 0;
			char* chrom_seq = NULL;
			char bin_seq_fp[1000];
			sprintf(bin_seq_fp, "%s/%s.bin", genome_seq_dir, chr_ids->at(i_chr));			
			if (check_file(bin_seq_fp))
			{
				chrom_seq = load_binary_sequence_file(bin_seq_fp, l_chrom_seq);
				fprintf(stderr, "Loaded %d nucleotides for the sequence.\n", l_chrom_seq);
			}
			else
			{
				fprintf(stderr, "Could not load genome sequence @ %s.\n", bin_seq_fp);
				delete[] signal_profile;
				delete[] multi_mapp_signal;
				continue;
			}

			vector<t_annot_region*>* cur_chr_valleys = get_significant_extrema_per_signal_profile(signal_data_dir, chr_id, signal_profile, l_profile,
																								multi_mapp_signal, l_multimapp_profile, chrom_seq,
																								extrema_statistic_defn);

			all_significant_valleys->insert(all_significant_valleys->end(), cur_chr_valleys->begin(), cur_chr_valleys->end());

			// Free memory.
			delete[] signal_profile;
			if (multi_mapp_signal != NULL)
			{
				delete[] multi_mapp_signal;
			}

			if (chrom_seq != NULL)
			{
				delete[] chrom_seq;
			}
		} // i_chr loop.

		if (all_significant_valleys->size() == 0)
		{
			fprintf(stderr, "No valleys found.\n");
			exit(0);
		}

		// Estimate FDR.
		get_benjamini_hochberg_corrected_p_values_per_valleys(all_significant_valleys);

		// Filter valleys
		vector<t_annot_region*>* filtered_all_significant_valleys = new vector<t_annot_region*>();
		for (int i_val = 0; i_val < (int)(all_significant_valleys->size()); i_val++)
		{
			t_valley_significance_info* significance_info = all_significant_valleys->at(i_val)->significance_info;

			// Significance check.
			if (significance_info->log_q_val < log_q_val_threshold)
			{
				filtered_all_significant_valleys->push_back(all_significant_valleys->at(i_val));
			}
		} // i_val loop.

		fprintf(stderr, "Detected %d significant valleys, saving.\n", (int)filtered_all_significant_valleys->size());

		char op_fp[1000];
		sprintf(op_fp, "%s/significant_valleys.bed.gz", signal_data_dir);
		dump_valleys(filtered_all_significant_valleys, op_fp);
	} // -get_significant_extrema_per_signal_profile option.
	else if (t_string::compare_strings(argv[1], "-merge_valleys"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -merge_valleys [Valleys BED file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* valleys_BED_fp = argv[2];
		int l_minima_vicinity = atoi(argv[3]);
		char* op_fp = argv[4];

		if (!check_file(valleys_BED_fp))
		{
			fprintf(stderr, "%s does not exist.\n", valleys_BED_fp);
			exit(0);
		}

		vector<t_annot_region*>* merged_valleys = merge_overlapping_valleys_per_pval_minimization(valleys_BED_fp, l_minima_vicinity);

		// Read the header.
		FILE* f_valleys = open_f(valleys_BED_fp, "r");
		char* header_line = getline(f_valleys);
		fclose(f_valleys);

		// Save header then write merged regions.
		FILE* f_op = open_f(op_fp, "w");
		if (header_line[0] == '#')
		{
			fprintf(f_op, "%s\n", header_line);
		}
		for (int i_reg = 0; i_reg < (int)merged_valleys->size(); i_reg++)
		{
			fprintf(f_op, "%s\n", (char*)(merged_valleys->at(i_reg)->data));
		} // i_reg loop.
		fclose(f_op);
	} // merge_valleys option.
	else if (t_string::compare_strings(argv[1], "-bspline_encode"))
	{
		if (argc < 3)
		{
			//fprintf(stderr, "USAGE: %s -bspline_encode [bedGraph/processed reads directory path] [# Spline Coefficients] [Spline Order (>=2)] [Breakpoints type (Uniform=0;Derivative=1;Vicinity_Derivative=2;Random=3)] [Max max error] [Max avg error] [window length] [Sparse data flag (0/1)] [Post Median Filter Length]\n", argv[0]);
			fprintf(stderr, "USAGE: %s -bspline_encode [Options] [Arguments]\n\
	-signal_dir [bedGraph/processed reads directory path]\n\
	-n_spline_coeff [# Spline Coefficients]\n\
	-bspline_order [Spline Order (>=2)]\n\
	-max_max_err [Max max error]\n\
	-max_avg_err [Max avg error]\n\
	-l_win [Window length]\n\
	-l_step_win [Stepping length (<l_win)]\n\
	-min_pt2pt_distance_in_bps [Distance between consecutive POIs]\n\
	-l_frag [Fragment length]\n\
	-sparse_profile [Sparse data flag (0/1)]\n\
	-brkpts_type [Breakpoints type (Uniform=0;Derivative=1;Vicinity_Derivative=2;Random=3)]\n\
	-l_post_filt_win [Post Median Filter Length]\n\
	-target_regs [Target regions BED file path]\n", argv[0]);

			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "-");

		bool success = false;
		char* signal_dir = cli->get_value_by_option("-signal_dir", success);
		//char* target_regs_BED_fp = (cli->get_value_by_option("-target_regs_BED_fp", success));if (!success) { target_regs_BED_fp = NULL; }
		int n_spline_coeff = atoi_null(cli->get_value_by_option("-n_spline_coeff", success));if (!success) { n_spline_coeff = 10; }
		int bspline_order = atoi_null(cli->get_value_by_option("-bspline_order", success));if (!success) { bspline_order = 5; }
		double max_max_err = atof_null(cli->get_value_by_option("-max_max_err", success));if (!success) { max_max_err = 5; }
		double max_avg_err = atof_null(cli->get_value_by_option("-max_avg_err", success));if (!success) { max_avg_err = 3; }
		int l_win = atoi_null(cli->get_value_by_option("-l_win", success));if (!success) { l_win = 1000; }

		bool sparse_data_flag = false;
		if (!cli->get_value_by_option("-sparse_profile", success)) { sparse_data_flag = false; }
		else { sparse_data_flag = (cli->get_value_by_option("-sparse_profile", success)[0] == '1'); }

		int brkpts_type_index = atoi_null(cli->get_value_by_option("-brkpts_type", success));if (!success) { brkpts_type_index = 0; }
		int l_med_filt_win = atoi_null(cli->get_value_by_option("-l_post_filt_win", success));if (!success) { l_med_filt_win = 50; }
		int l_step_win = atoi_null(cli->get_value_by_option("-l_step_win", success));if (!success) { l_step_win = l_win; }
		int min_pt2pt_distance_in_bps = atoi_null(cli->get_value_by_option("-min_POI_distance", success));if (!success) { min_pt2pt_distance_in_bps = 50; }
		int l_frag = atoi_null(cli->get_value_by_option("-l_frag", success));if (!success) { l_frag = 200; }

		if (l_step_win > l_win)
		{
			fprintf(stderr, "Stepping window length cannot be longer than l_win.\n");
			l_step_win = l_win;
		}

		fprintf(stderr, "Parameters:\n\
signal_dir=%s\n\
n_spline_coeff=%d\n\
bspline_order=%d\n\
brkpts_type_index=%d\n\
max_max_err=%.3f\n\
max_avg_err=%.3f\n\
l_win=%d\n\
l_frag=%d\n\
min_pt2pt_distance_in_bps=%d\n\
sparse_data_flag=%d\n\
l_step_win=%d\n\
l_med_filt_win=%d\n",
			signal_dir,
			n_spline_coeff,
			bspline_order,
			brkpts_type_index,
			max_max_err,
			max_avg_err,
			l_win,
			l_frag,
			min_pt2pt_distance_in_bps,
			sparse_data_flag,
			l_step_win,
			l_med_filt_win);

		char brkpts_type = -1;
		if (brkpts_type_index == 0)
		{
			fprintf(stderr, "Uniform breakpoints.\n");
			brkpts_type = UNIFORM_BREAKPOINTS;
		}
		else if (brkpts_type_index == 1)
		{
			fprintf(stderr, "Derivative-based breakpoints.\n");
			brkpts_type = DERIVATIVE_NU_BREAKPOINTS;
		}
		else if (brkpts_type_index == 2)
		{
			fprintf(stderr, "Vicinity derivative-based breakpoints.\n");
			brkpts_type = VICINITY_DERIVATIVE_NU_BREAKPOINTS;
		}
		else if (brkpts_type_index == 3)
		{
			fprintf(stderr, "Random breakpoints.\n");
			brkpts_type = RANDOM_NU_BREAKPOINTS;
		}

		fprintf(stderr, "Breakpoints type: %d\n", brkpts_type);

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.txt", signal_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		if (chr_ids == NULL)
		{
			fprintf(stderr, "Could not load chromosome id's from %s.\n", chr_ids_list_fp);
			exit(0);
		}

		//int l_frag = 200;
		int min_n_pts_2_encode = 2;
		//int min_pt2pt_distance_in_bps = 50;

		double top_err_perc_frac = 0.90;

		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			if(sparse_data_flag)
			{ 
				fprintf(stderr, "Encoding %s (Sparse)\n", chr_ids->at(i_chr));
			}
			else
			{
				fprintf(stderr, "Encoding %s\n", chr_ids->at(i_chr));
			}
			
			//bspline_encode_mapped_read_profile_expanded_windows(signal_dir, chr_ids->at(i_chr), l_frag,
			bspline_encode_mapped_read_profile(signal_dir, chr_ids->at(i_chr), l_frag,
				n_spline_coeff, bspline_order,
				brkpts_type,
				min_n_pts_2_encode, min_pt2pt_distance_in_bps,
				max_max_err, max_avg_err,
				l_win, sparse_data_flag, top_err_perc_frac, l_step_win, l_med_filt_win);
		} // i_chr loop.
	} // -bspline_encode_mapped_read_profile option.
	else if (t_string::compare_strings(argv[1], "-append_signal_2_regions"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -append_signal_2_regions [Signals directory] [Fragment length] [BED file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* signal_directory = argv[2];
		int l_frag = atoi(argv[3]);
		char* bed_fp = argv[4];
		char* op_fp = argv[5];

		append_signal_2_regions(signal_directory, l_frag, bed_fp, op_fp);
	} // -append_signal_2_regions option.

	FILE* f_beacon = open_f("timing.log", "a");
	clock_t end_c = clock();
	fprintf(f_beacon, "%d\n", (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fprintf(stderr, "EpiSAFARI finished \"%s\" in %d seconds.\n", argv[1], (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fclose(f_beacon);
}

