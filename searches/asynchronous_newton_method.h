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

	POPULATION *current_population;
} NEWTON_METHOD_SEARCH;

int fwrite_newton_method(FILE* file, NEWTON_METHOD_SEARCH *nms);
int write_newton_method(char* file, NEWTON_METHOD_SEARCH *nms);

int fread_newton_method(FILE* file, NEWTON_METHOD_SEARCH **nms);
int read_newton_method(char* file, NEWTON_METHOD_SEARCH **nms);

int init_newton_method(char* search_name, SEARCH* search);

#endif
