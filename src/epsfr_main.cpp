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
	BedGraph Processing:\n\
		-separate_bedGraph_2_chromosomes [bedGraph file path] [Output directory]\n\
Signal Feature Detection:\n\
	-bspline_encode [bedGraph/processed reads directory path] [# Spline Coefficients] [Spline Order (>=2)] [Max max error] [Max avg error] [window length] [Sparse data flag (0/1)] [Post Median Filter Length]\n\
	-get_significant_extrema_per_signal_profile [Binary signal profile file path] [chr id] \
[Maximum signal at trough] [Minimum signal at summit] \
[Minimum summit2trough ratio per trough] [Maximum summit2trough distance in bps] \
[Multi-mappability profile directory] [Maximum multimapp signal at trough]\n\
Feature Annotation:\n\
	-annotate_features [Signal directory] [GFF file path] [Half promoter length] [Output file path]\n\n",	argv[0]);
}

int main(int argc, char* argv[])
{
	clock_t start_c = clock();

	if (argc < 3)
	{
		print_usage(argv);
		exit(0);
	}

	if (strcmp(argv[1], "-help") == 0 ||
		strcmp(argv[1], "-version") == 0 ||
		strcmp(argv[1], "-h") == 0 ||
		strcmp(argv[1], "-v") == 0)
	{
		print_usage(argv);
		exit(0);
	}
	else if (strcmp(argv[1], "-annotate_features") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s -annotate_features [Signal directory] [GFF file path] [Half promoter length] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* signal_directory = argv[2];
		char* gff_fp = argv[3];
		int l_prom = atoi(argv[4]);
		char* op_fp = argv[5];

		annotate_features(signal_directory, gff_fp, l_prom/2, op_fp);
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
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_chr));

			FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");

			// Load the entries, store them in a buffer.
			//FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");
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

			fclose(f_cur_chr_reads);

			for (int i_st = 0; i_st < (int)bucket_f_list->size(); i_st++)
			{
				fclose(bucket_f_list->at(i_st));
			} // i_st loop.
			delete(bucket_f_list);

			char cur_chr_sorted_reads_fp[1000];
			sprintf(cur_chr_sorted_reads_fp, "%s/%s_mapped_reads.txt", sorted_reads_op_dir, chr_ids->at(i_chr));
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
			fclose(f_cur_chr_sorted_reads);

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
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", sorted_preprocessed_reads_dir, chr_ids->at(i_chr));
			FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");

			char cur_chr_pruned_reads_fp[1000];
			sprintf(cur_chr_pruned_reads_fp, "%s/%s_mapped_reads.txt", pruned_reads_dir, chr_ids->at(i_chr));
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
			fclose(f_cur_chr_reads);
			fclose(f_cur_chr_pruned_reads);

			fprintf(stderr, "Processed %d reads, pruned to %d reads.\n", n_processed_reads, n_pruned_reads);
		} // i_chr loop.
	} // -remove_duplicates
	else if (strcmp(argv[1], "-separate_bedGraph_2_chromosomes") == 0)
	{
		if (argc != 4)
		{
			fprintf(stderr, "USAGE: %s -separate_bedGraph_2_chromosomes [bedGraph file path] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* bedgraph_fp = argv[2];
		char* op_dir = argv[3];

		vector<char*>* chr_ids = new vector<char*>();
		vector<FILE*>* f_bgr_op_list = new vector<FILE*>();

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
				sprintf(bgr_op_fp, "%s/%s.bgr", op_dir, copy_chr_id);

				if (check_file(bgr_op_fp))
				{
					fprintf(stderr, "%s exists, pooling on it.\n", bgr_op_fp);
				}

				FILE* cur_f_bgr_op = open_f(bgr_op_fp, "a");
				f_bgr_op_list->push_back(cur_f_bgr_op);
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
		fclose(f_bgr);

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
			fclose(f_bgr_op_list->at(i_chr));
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
	if (t_string::compare_strings(argv[1], "-get_significant_extrema"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s -get_significant_extrema [Encoded signals data directory path] \
[Maximum signal at trough] [Minimum signal at summit] \
[Minimum summit2trough ratio per trough] [Maximum summit2trough distance in bps] \
[Multi-mappability profile directory] [Maximum multimapp signal at trough]\n", argv[0]);
			exit(0);
		}

		char* signal_data_dir = argv[2];		
		double max_signal_at_trough = atof(argv[3]);
		double min_signal_at_summit = atof(argv[4]);
		double min_summit2trough_ratio_per_trough = atof(argv[5]);
		int max_summit2trough_dist_in_bp = atoi(argv[6]);
		char* multi_mapp_signal_profile_dir = argv[7];
		double max_multimapp_signal_at_trough = atof(argv[8]);

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
		extrema_statistic_defn->max_summit2trough_dist_in_bp = max_summit2trough_dist_in_bp;
		extrema_statistic_defn->max_multimapp_signal_at_trough = max_multimapp_signal_at_trough;

		// Process all the chromosomes.
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			char* chr_id = chr_ids->at(i_chr);

			fprintf(stderr, "Detecting the features on %s\n", chr_id);

			char bin_signal_fp[1000];
			sprintf(bin_signal_fp, "%s/spline_coded_%s.bin.gz", signal_data_dir, chr_id);

			int l_profile = 0;
			double* signal_profile = load_per_nucleotide_binary_profile(bin_signal_fp, l_profile);
			fprintf(stderr, "Loaded %d long signal profile.\n", l_profile);

			char multi_mapp_signal_profile_fp[1000];
			sprintf(multi_mapp_signal_profile_fp, "%s/%s.bin", multi_mapp_signal_profile_dir, chr_id);
			int l_multimapp_profile = 0;
			double* multi_mapp_signal = load_normalized_multimappability_profile(multi_mapp_signal_profile_fp, l_multimapp_profile);
			fprintf(stderr, "Loaded %d long multimappability profile.\n", l_multimapp_profile);

			get_significant_extrema_per_signal_profile(signal_data_dir, chr_id, signal_profile, l_profile,
				multi_mapp_signal, l_multimapp_profile,
				extrema_statistic_defn);

			delete[] multi_mapp_signal;
			delete[] signal_profile;
		} // i_chr loop.
	} // -get_significant_extrema_per_signal_profile option.
	else if (t_string::compare_strings(argv[1], "-bspline_encode"))
	{
		if (argc != 10)
		{
			fprintf(stderr, "USAGE: %s -bspline_encode [bedGraph/processed reads directory path] [# Spline Coefficients] [Spline Order (>=2)] [Max max error] [Max avg error] [window length] [Sparse data flag (0/1)] [Post Median Filter Length]\n", argv[0]);
			exit(0);
		}

		char* signal_dir = argv[2];
		int n_spline_coeff = atoi(argv[3]);
		int bspline_order = atoi(argv[4]);
		double max_max_err = atof(argv[5]);
		double max_avg_err = atof(argv[6]);
		int l_win = atoi(argv[7]);
		bool sparse_data_flag = (argv[8][0] == '1');
		int l_med_filt_win = atoi(argv[9]);

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.txt", signal_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		if (chr_ids == NULL)
		{
			fprintf(stderr, "Could not load chromosome id's from %s.\n", chr_ids_list_fp);
			exit(0);
		}

		int l_frag = 200;
		int min_n_pts_2_encode = 2;
		int min_pt2pt_distance_in_bps = 50;

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

			bspline_encode_mapped_read_profile(signal_dir, chr_ids->at(i_chr), l_frag,
				n_spline_coeff, bspline_order, min_n_pts_2_encode, min_pt2pt_distance_in_bps,
				max_max_err, max_avg_err,
				l_win, sparse_data_flag, top_err_perc_frac, l_med_filt_win);
		} // i_chr loop.
	} // -bspline_encode_mapped_read_profile option.

	//FILE* f_beacon = open_f("beacon.log", "a");
	clock_t end_c = clock();
	//fprintf(f_beacon, "EpiSAFARI finished in %d seconds.\n", (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fprintf(stderr, "EpiSAFARI finished in %d seconds.\n", (int)((end_c - start_c) / CLOCKS_PER_SEC));
	//fclose(f_beacon);
}

