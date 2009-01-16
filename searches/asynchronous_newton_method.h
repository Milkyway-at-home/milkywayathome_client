#ifndef FGDO_ASYNCHRONOUS_NEWTON_H
#define FGDO_ASYNCHRONOUS_NEWTON_H

#include <stdio.h>

#include "../evaluation/search_manager.h"
#include "population.h"
#include "asynchronous_search.h"

#define NEWTON_FIXED_RANGE 0
#define NEWTON_ERROR_RANGE 1
#define NEWTON_UPDATE_RANGE 2
#define NEWTON_LS_FIXED_RANGE 3
#define NEWTON_LS_ERROR_RANGE 4
#define NEWTON_LS_UPDATE_RANGE 5

#define NMS_HESSIAN 1
#define NMS_LINE_SEARCH 2

typedef struct line_search {
	int iteration, max_iteration;
	int evaluation;

	double center, min_range, max_range;
} LINE_SEARCH;

typedef struct newton_method_search {
	int type, mode;

	int number_parameters;
	double *current_point, *previous_point;
	double *direction;
	double *parameter_range;
	double *min_bound, *max_bound;

	int current_iteration, maximum_iteration;
	int current_evaluation, evaluations_per_iteration;

	LINE_SEARCH* line_search;

	POPULATION *population;
} NEWTON_METHOD_SEARCH;

ASYNCHRONOUS_SEARCH* get_asynchronous_newton_method();

int create_newton_method(char* search_name, int type, int line_search, int maximum_iteration, int evaluations_per_iteration, int number_parameters, double *point, double *range, double* min_bound, double* max_bound);
int read_newton_method(char* search_name, void** search_data);
int checkpoint_newton_method(char* search_name, void* search_data);
int newton_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);
int newton_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);

#endif
