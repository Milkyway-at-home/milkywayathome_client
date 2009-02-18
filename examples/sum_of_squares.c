#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../searches/asynchronous_newton_method.h"
#include "../searches/asynchronous_particle_swarm.h"
#include "../searches/synchronous_gradient_descent.h"
#include "../searches/synchronous_newton_method.h"
#include "../searches/search_arguments.h"
#include "../evaluation/simple_evaluator.h"
#include "../evaluation/search_manager.h"

double *min_bound, *max_bound;
int number_parameters;

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

	srand48(23);

	number_parameters = atoi(arguments[1]);
	min_bound = (double*)malloc(sizeof(double) * number_parameters);
	max_bound = (double*)malloc(sizeof(double) * number_parameters);
	range = (double*)malloc(sizeof(double) * number_parameters);
	point = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		min_bound[i] = -100.0;
		max_bound[i] = 100.0;
		range[i] = 2.0;
		point[i] = (drand48() * (double)(max_bound[i] - min_bound[i])) + (double)min_bound[i];
	}
	init_simple_evaluator(sum_of_squares);

	/********
		*       Start the search
	 ********/
	if (argument_exists("-asynch", number_arguments, arguments)) {
		if (argument_exists("-nm", number_arguments, arguments))	register_search(get_asynchronous_newton_method());
//		if (argument_exists("-gs", number_arguments, arguments))	register_search(get_asynchronous_genetic_search());
//		if (argument_exists("-de", number_arguments, arguments))	register_search(get_asynchronous_differential_evolution());
		if (argument_exists("-ps", number_arguments, arguments))	register_search(get_asynchronous_particle_swarm());
//		if (argument_exists("-sx", number_arguments, arguments))	register_search(get_asynchronous_simplex());
		else {
			printf("Search not specified.\n");
			return 0;
		}
		asynchronous_search(number_arguments, arguments, number_parameters, point, range, min_bound, max_bound);
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
