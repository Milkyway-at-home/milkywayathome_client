#ifndef FGDO_ASYNCHRONOUS_NEWTON_H
#define FGDO_ASYNCHRONOUS_NEWTON_H

#include <stdio.h>

#include "../evaluation/search_manager.h"
#include "population.h"

typedef struct newton_method_search {
	int number_parameters;
	double *parameters;
	double *parameter_range;
	double *min_parameters, *max_parameters;

	int current_iteration, maximum_iteration;
	int current_evaluation, evaluations_per_iteration;

	POPULATION *population;
} NEWTON_METHOD_SEARCH;

int create_newton_method(char* search_name, ...);
int read_newton_method(char* search_name, void** search_data);
int checkpoint_newton_method(char* search_name, void* search_data);
int newton_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);
int newton_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);

const extern ASYNCHRONOUS_SEARCH asynchronous_newton_method;

#endif
