#ifndef __MAPPED_READ_FILE_INTERFACE__
#define __MAPPED_READ_FILE_INTERFACE__

#include <vector>
using namespace std;

const int MEG_BASE = 1000 * 1000;
const int K_BASE = 1000;

class t_string;
struct t_annot_region;

struct t_read_line_sorting_info
{
	int start;
	char* read_line;
};

unsigned char* load_normalized_multimappability_profile(char* mapability_signal_profile_fp, int& l_mapability_profile);

#define MAX_N_PAIRS (10)

bool check_genome_index_update_per_CIGAR_entry(char entry_char);
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char);

void preprocess_mapped_reads_file(char* mrf_fp, char* parsed_reads_op_dir, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	bool dump_read_ids);
void preprocess_tagAlign_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_ELAND_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_bowtie_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_BED4_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_BED5_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_BED6_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str);
void preprocess_BED456_read_line(char* cur_line,
	char* read_id,
	char* chrom,
	int& chr_index, int& sequenced_length,
	char& strand_char,
	char* mapping_quality_str);

void preprocess_PE_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	bool& first_segment_in_template,
	bool& last_segment_in_template,
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	int& mapping_quality,
	char* cigar_str);

void buffer_per_nucleotide_profile_no_buffer(const char* sorted_read_fp, const int l_extended_tag, 
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int l_buffer, int& l_data);

void get_per_strand_read_stats(vector<t_annot_region*>* regions,
	char* preprocessed_reads_dir);

bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced);
void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										int& l_cur_entry,
										char& entry_type_char);

vector<char*>* sort_bucket_read_lines(char* bucket_fp);

#endif // __MAPPED_READ_FILE_INTERFACE__

