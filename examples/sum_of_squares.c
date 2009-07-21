#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../searches/bounds.h"
#include "../searches/asynchronous_genetic_search.h"
#include "../searches/asynchronous_particle_swarm.h"
#include "../searches/synchronous_gradient_descent.h"
#include "../searches/synchronous_newton_method.h"
#include "../searches/search_arguments.h"
#include "../evaluation/simple_evaluator.h"
#include "../evaluation/search_manager.h"

int number_parameters;
double *min_bound;
double *max_bound;

#define ERROR_RATE 0.0

double sum_of_squares(double* parameters) {
	int i;
	double sum;

	if (drand48() < ERROR_RATE) {
		sum = 0.0;
		for (i = 0; i < number_parameters; i++) {
			sum += (max_bound[i] - min_bound[i]) * drand48() * (max_bound[i] - min_bound[i]);
		}
	} else {
		sum = 0.0;
		for (i = 0; i < number_parameters; i++) {
			sum += parameters[i] * parameters[i];
		}
	}
	return -sum;
}

int main(int number_arguments, char** arguments) {
	int i;
	double* range;
	double* point;
	int* in_radians;
	BOUNDS* bounds;

	srand48(23);

	number_parameters = atoi(arguments[1]);
	in_radians = (int*)malloc(sizeof(int) * number_parameters);
	min_bound = (double*)malloc(sizeof(double) * number_parameters);
	max_bound = (double*)malloc(sizeof(double) * number_parameters);
	range = (double*)malloc(sizeof(double) * number_parameters);
	point = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		in_radians[i] = 0;
		min_bound[i] = -100.0;
		max_bound[i] = 100.0;
		range[i] = 2.0;
		point[i] = (drand48() * (double)(max_bound[i] - min_bound[i])) + (double)min_bound[i];
	}
	init_simple_evaluator(sum_of_squares);

	new_bounds(&bounds, number_parameters, min_bound, max_bound, in_radians);

	/********
		*       Start the search
	 ********/
	if (argument_exists("-asynch", number_arguments, arguments)) {
//		if (argument_exists("-nm", number_arguments, arguments))	register_search(get_asynchronous_newton_method());
		if (argument_exists("-gs", number_arguments, arguments))	register_search(get_asynchronous_genetic_search());
//		else if (argument_exists("-de", number_arguments, arguments))	register_search(get_asynchronous_differential_evolution());
		else if (argument_exists("-ps", number_arguments, arguments))	register_search(get_asynchronous_particle_swarm());
//		else if (argument_exists("-sx", number_arguments, arguments))	register_search(get_asynchronous_simplex());
		else {
			printf("Search not specified.\n");
			return 0;
		}
		asynchronous_search(number_arguments, arguments, number_parameters, point, range, bounds);
	} else {
		if (argument_exists("-nm", number_arguments, arguments))	synchronous_newton_method(number_arguments, arguments, number_parameters, point, range);
		else if (argument_exists("-gd", number_arguments, arguments))	synchronous_gradient_descent(number_arguments, arguments, number_parameters, point, range);
		else if (argument_exists("-cgd", number_arguments, arguments))	synchronous_conjugate_gradient_descent(number_arguments, arguments, number_parameters, point, range);
		else {
			return 0;
		}
	}

	return 0;
}
