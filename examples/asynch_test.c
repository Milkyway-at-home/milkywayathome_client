#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../searches/asynchronous_search.h"
#include "../searches/asynchronous_newton_method.h"

#include "../evaluation/search_manager.h"
#include "../evaluation/mpi_search_manager.h"
#include "../evaluation/simple_evaluator.h"


double *min_bound, *max_bound;
int number_parameters;


double sum_of_squares(double* parameters) {
	int i;
	double sum;

	sum = 0.0;
	for (i = 0; i < number_parameters; i++) {
		sum += parameters[i] * parameters[i];
	}
	return sum;
}

int main(int argc, char** argv) {
	int i;
	double* step;
	double* point;

	srand48(23);

	number_parameters = atoi(argv[1]);
	min_bound = (double*)malloc(sizeof(double) * number_parameters);
	max_bound = (double*)malloc(sizeof(double) * number_parameters);
	step = (double*)malloc(sizeof(double) * number_parameters);
	point = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		min_bound[i] = -100.0;
		max_bound[i] = 100.0;
		step[i] = 2.0;
		point[i] = (drand48() * (double)(max_bound[i] - min_bound[i])) + (double)min_bound[i];
	}

	init_simple_evaluator(sum_of_squares);
	register_search(get_asynchronous_newton_method());

	for (i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-s")) {
			char search_name[SEARCH_NAME_SIZE], *search_qualifier;
			strcpy(search_name, argv[++i]);
			get_qualifier_from_name(search_name, &search_qualifier);

			if (!strcmp(search_qualifier, "nm")) {
				printf("creating newton method...\n");
				create_newton_method(search_name, NEWTON_PARABOLIC_LINE_SEARCH, 3, 300, number_parameters, point, step, min_bound, max_bound);
				printf("created.\n");
			} else if (!strcmp(search_qualifier, "gs")) {
			} else if (!strcmp(search_qualifier, "de")) {
			} else if (!strcmp(search_qualifier, "pso")) {
			}
			free(search_qualifier);
		}
	}
	start_mpi_search_manager(argc, argv);

	return 0;
}
