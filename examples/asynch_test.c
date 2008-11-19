#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../searches/asynchronous_search.h"
#include "../searches/asynchronous_newton_method.h"

#include "../evaluation/search_manager.h"
#include "../evaluation/mpi_search_manager.h"
#include "../evaluation/simple_evaluator.h"


double *min_parameters, *max_parameters;
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
	min_parameters = (double*)malloc(sizeof(double) * number_parameters);
	max_parameters = (double*)malloc(sizeof(double) * number_parameters);
	step = (double*)malloc(sizeof(double) * number_parameters);
	point = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		min_parameters[i] = -100.0;
		max_parameters[i] = 100.0;
		step[i] = 2.0;
		point[i] = (drand48() * (double)(max_parameters[i] - min_parameters[i])) + (double)min_parameters[i];
	}

	init_simple_evaluator(sum_of_squares);
	register_search(asynchronous_newton_method);

	for (i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-s")) {
			char search_name[SEARCH_NAME_SIZE], *search_qualifier;
			strcpy(search_name, argv[++i]);
			get_qualifier_from_name(search_name, &search_qualifier);

			if (!strcmp(search_qualifier, "nm")) {
				printf("creating newton method...\n");
				create_newton_method(search_name, 10, 100, number_parameters, point, step);
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
