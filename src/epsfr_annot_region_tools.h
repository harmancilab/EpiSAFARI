#ifndef __ANNOT_REGION_TOOLS__
#define  __ANNOT_REGION_TOOLS__

#include <vector>

using namespace std;

class t_rng;

struct t_valley_significance_info;

// Stores the sorting information for the regions.
struct t_sorting_info
{
	// Following are the cumulative ends for the region: Basically, for a current read, sets the smallest and largest start and end at this read position.
	// This becomes useful when doing binary searches over the regions that overlap.
	int cumulative_sorted_start;
	int cumulative_sorted_end;	
};

// This is just a region on a chromosome.
struct t_annot_region
{
	char* chrom;
	char strand;
	int start; // If strand is '-', this coordinate is with respect to reverse complement on the forward strand.
	unsigned int score;
	double dbl_score;
	char* name;
	int end;

	// Load the remaining information if it is available.
	int n_exons;
	int thick_start;
	int thick_end;
	vector<t_annot_region*>* intervals; // These are the exons in the region. Some bed files include this information.

	t_valley_significance_info* significance_info;

	t_sorting_info* sort_info;
	void* annotation_info;
	void* data; // Extra data about the region, this is a pointer to a data structure that holds specific data about the region depending on the loaded data type: narrowPeak, GFF, ... Can also be used for other purposes for storing data about a region.
};

struct t_intersect_info
{
	t_annot_region* src_reg;
	t_annot_region* dest_reg;
	int l_overlap;
};

void sort_set_sorting_info(vector<t_annot_region*>* regions, bool (sort_regions_callback)(t_annot_region*, t_annot_region*));

struct t_restr_annot_region_list
{
	vector<t_annot_region*>** regions_per_chrom;
	vector<t_annot_region*>** pos_strand_regions_per_chrom;
	vector<t_annot_region*>** neg_strand_regions_per_chrom;
	vector<char*>* chr_ids; // These are the ids that have smae ordering in the region lists.
};

// Following is the extra information that is stored in a narrowPeak file. This information is stored for each narrowPeak entry.
#define NARROWPEAK_INFO (0x1234)
struct t_narrowPeak_info
{
	int info_type;
	int score;
	char strand;
	double signal_value;
	double p_value;
	double q_value;
	int peak_pos; // This is usually important. This is where the peak is.
};

t_restr_annot_region_list* restructure_annot_regions(vector<t_annot_region*>* regions); // Re-structure the list using the chromosome ids retreived from itself.
t_restr_annot_region_list* restructure_annot_regions(vector<t_annot_region*>* regions, vector<char*>* chr_ids); // Re-structure the list for a given list of chromosome ids.
void delete_restructured_annot_regions(t_restr_annot_region_list* restructured_region_lists);

vector<t_annot_region*>* load_BED(char* bed_fp);
vector<t_annot_region*>* load_BED12(char* bed_fp);
vector<t_annot_region*>* load_BED_with_line_information(char* bed_fp);
void dump_BED(const char* bed_fp, vector<t_annot_region*>* annot_regions);
vector<t_annot_region*>* load_narrowPeak(char* narrowPeak_fp); // This function also loads the broadPeak files.
void extend_BED(vector<t_annot_region*>* regions, int l_extend, bool extend_3p, bool extend_5p);

bool validate_region_coords(vector<t_annot_region*>* regions);

//t_annot_region* copy_region(t_annot_region* region);

void load_chromosome_lengths_per_tabbed_file(char* chr_lengths_fp, vector<char*>* chr_ids, vector<int>* chr_lengths);

t_annot_region* get_empty_region();

// Divide into strands first, then divide into chromosomes, then sort.
vector<t_annot_region*>* sort_regions_per_chromosome_per_strand(vector<t_annot_region*>* annot_region_list);

void delete_annot_regions(vector<t_annot_region*>* region_list);
void delete_annot_regions(t_annot_region* region);
void delete_annot_regions_with_line_information(vector<t_annot_region*>* region_list);

bool check_line_skip(char* cur_line);

vector<char*>* get_chr_ids(vector<t_annot_region*>* annot_regions);

vector<t_annot_region*>* get_regions_per_chromosome(vector<t_annot_region*>* annot_regions,
													char* chr_id);

vector<t_annot_region*>* get_regions_per_strand(vector<t_annot_region*>* annot_regions,
													char strand);

/*
Region intersection
*/
vector<t_annot_region*>* intersect_annot_regions(vector<t_annot_region*>* annot_regions1,
														vector<t_annot_region*>* annot_regions2,
														bool match_strands,
														bool find_all_overlaps);

vector<t_annot_region*>* intersect_annot_regions(vector<t_annot_region*>* annot_regions1,
													vector<t_annot_region*>* annot_regions2,
													bool find_all_overlaps);

void delete_intersect_info(vector<t_annot_region*>* regions);

/*
Region merging
*/
vector<t_annot_region*>* merge_annot_regions(vector<t_annot_region*>* annot_regions1,
												vector<t_annot_region*>* annot_regions2,
												int max_gap, 
												bool match_strands);

vector<t_annot_region*>* merge_annot_regions(vector<vector<t_annot_region*>*>* annot_regions_list,
												int max_gap, 
												bool match_strands);

vector<t_annot_region*>* merge_annot_regions(vector<t_annot_region*>* total_annot_regions,
											int max_gap, 
											bool match_strands);

vector<t_annot_region*>* merge_annot_regions(vector<t_annot_region*>* total_annot_regions,
											int max_gap);

/*
Region exclusion.
*/
vector<t_annot_region*>* exclude_annot_regions(vector<t_annot_region*>* annot_regions1,
													vector<t_annot_region*>* annot_regions2, 
													bool strand_specific);

// Backend functions.
vector<t_annot_region*>* exclude_annot_regions(vector<t_annot_region*>* annot_regions1,
													vector<t_annot_region*>* annot_regions2);

/*
Region subtraction.
*/
vector<t_annot_region*>* subtract_annot_regions(t_annot_region* region1, t_annot_region* region2, bool strand_specific);

// Necessary for fast comparison of the regions.
bool sort_regions(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_increasing_length(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_decreasing_length(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_ends(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_score(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_dbl_score_decreasing(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_dbl_score(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_increasing_p_value(t_annot_region* reg1, t_annot_region* reg2);
bool sort_regions_per_increasing_q_value(t_annot_region* reg1, t_annot_region* reg2);

// Duplicate a region.
t_annot_region* duplicate_region(t_annot_region* region_2_dup);

vector<t_annot_region*>* subtract_annot_regions(vector<t_annot_region*>* regions1, vector<t_annot_region*>* regions2);

// Locate the index of the region (out of a list of regions) whose start is just to the left of start_posn.
int locate_posn_region_per_region_starts(int start_posn, vector<t_annot_region*>* region_list, int i, int j);
int locate_posn_per_sorted_posn_list(int posn, vector<int>* sorted_posn_list, int i, int j);

int region_3p_accessor(void* void_ptr);
int region_5p_accessor(void* void_ptr);

#endif // __ANNOT_REGION_TOOLS__
