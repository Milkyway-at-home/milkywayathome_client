#include "search_parameters.h"

typedef struct heap_node {
	double time;
	SEARCH_PARAMETERS *sp;
} HEAP_NODE;

void simulation_heap__init(int max_size, int number_parameters, char* time_template_filename);
void simulation_heap__insert(SEARCH_PARAMETERS *sp);
void simulation_heap__remove_min(SEARCH_PARAMETERS *sp);
