#ifndef __MIN_MAX_UTILS__
#define __MIN_MAX_UTILS__

#include <vector>
using namespace std;

struct t_annot_region;

struct t_plateau
{
	int start;
	int end;
};

enum{EXTREMA_MAX, EXTREMA_MIN};
struct t_extrema_node
{
	int extrema_posn;
	int extrema_type;

	// This is the scale at which this extrema is. Each extrema is found at a scale. 
	int scale;
	
	// Height at the extrema value.
	double height_at_extrema;
};

void delete_extrema_nodes(vector<t_extrema_node*>* extrema_nodes);

void get_extrema_per_plateaus(double* data, int l_signal,
	vector<t_extrema_node*>* maxes, 
	vector<t_extrema_node*>* mins,
	int* derivative_sign_map,
	int i_scale,
	double zero_deriv);

void get_extrema_per_plateaus(double* data, int l,
	vector<t_extrema_node*>* maxes, 
	vector<t_extrema_node*>* mins,
	int* derivative_sign_map,
	int i_scale);

bool sort_extremas_per_posn(t_extrema_node* node1, t_extrema_node* node2);
bool sort_extremas_per_height(t_extrema_node* node1, t_extrema_node* node2);
bool sort_extremas_per_decreasing_height(t_extrema_node* node1, t_extrema_node* node2);

#endif // __MIN_MAX_UTILS__
