#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../searches/bounds.h"
#include "../searches/asynchronous_genetic_search.h"
#include "../searches/asynchronous_particle_swarm.h"
#include "../searches/asynchronous_differential_evolution.h"
#include "../searches/synchronous_gradient_descent.h"
#include "../searches/synchronous_newton_method.h"
#include "../searches/search_arguments.h"
#include "../evaluation/simple_evaluator.h"
#include "../evaluation/search_manager.h"

int number_parameters;
double *min_bound;
double *max_bound;

#define ERROR_RATE 0.0

double sphere(double* parameters) {
	int i;
	double sum;

	if (drand48() < ERROR_RATE) {
		sum = 0.0;
		for (i = 0; i < number_parameters; i++) {
			sum += (max_bound[i] - min_bound[i]) * (max_bound[i] - min_bound[i]) * drand48();
		}
	} else {
		sum = 0.0;
		for (i = 0; i < number_parameters; i++) {
			sum += parameters[i] * parameters[i];
		}
	}
	return -sum;
}

double ackley(double* parameters) {
	int i;
	double sum1, sum2;
	double *x;

	x = (double*)malloc(sizeof(double) * number_parameters);
	if (drand48() < ERROR_RATE) {
		for (i = 0; i < number_parameters; i++) {
			x[i] = (max_bound[i] - min_bound[i]) * drand48();
		}
	} else {
		for (i = 0; i < number_parameters; i++) {
			x[i] = parameters[i];
		}
	}

	sum1 = 0.0;
	for (i = 0; i < number_parameters; i++) {
		sum1 = x[i] * x[i];
	}
	sum1 /= (number_parameters);
	sum1 = -0.2 * sqrt(sum1);

	sum2 = 0.0;
	for (i = 0; i < number_parameters; i++) {
		sum2 += cos(2 * M_PI * x[i]);
	}
	sum2 /= (number_parameters);
	free(x);

	return -(-20 * exp(sum1) - exp(sum2) + 20 + M_E);
}

double griewank(double* parameters) {
	int i;
	double sum1, sum2;
	double *x;

	x = (double*)malloc(sizeof(double) * number_parameters);
	if (drand48() < ERROR_RATE) {
		for (i = 0; i < number_parameters; i++) {
			x[i] = (max_bound[i] - min_bound[i]) * drand48();
		}
	} else {
		for (i = 0; i < number_parameters; i++) {
			x[i] = parameters[i];
		}
	}

	sum1 = 0.0;
	for (i = 0; i < number_parameters; i++) {
		sum1 += x[i] * x[i];
	}
	sum1 /= 4000;

	sum2 = cos( x[0] / sqrt(1) );
	for (i = 1; i < number_parameters; i++) {
		sum2 *= cos( x[i] / sqrt(i+1) );
	}
	free(x);

	return -(sum1 - sum2 + 1);
}

double rosenbrock(double* parameters) {
	int i;
	double sum, tmp;
	double *x;

	x = (double*)malloc(sizeof(double) * number_parameters);
	if (drand48() < ERROR_RATE) {
		for (i = 0; i < number_parameters; i++) {
			x[i] = (max_bound[i] - min_bound[i]) * drand48();
		}
	} else {
		for (i = 0; i < number_parameters; i++) {
			x[i] = parameters[i];
		}
	}

	sum = 0.0;
	for (i = 0; i < number_parameters-1; i++) {
		tmp = (x[i+1] - (x[i] * x[i]));
		sum += (100 * tmp * tmp) + ((x[i] - 1) * (x[i] - 1));
	}

	return -sum;
}


double rastrigin(double* parameters) {
	int i;
	double sum;
	double xi;

	sum = 0.0;
	if (drand48() < ERROR_RATE) {
		for (i = 0; i < number_parameters; i++) {
			xi = drand48() * (max_bound[i] - min_bound[i]);
			sum += (xi * xi) - (10 * cos(2 * M_PI * xi)) + 10;
		}
	} else {
		for (i = 0; i < number_parameters; i++) {
			xi = parameters[i];
			sum += (xi * xi) - (10 * cos(2 * M_PI * xi)) + 10;
		}
	}
	return -sum;
}


int main(int number_arguments, char** arguments) {
	int i;
	double* range;
	double* point;
	double min_value, max_value;
	int* in_radians;
	BOUNDS* bounds;

	if (!strcmp(arguments[0], "-?") || !strcmp(arguments[0], "-help") || !strcmp(arguments[0], "-h")) {
		printf("usage: ./sum_of_squares <equation> <number_parameters> ...\n");
		exit(0);
	}

	srand48(23);


	if (!strcmp(arguments[1], "sphere")) {
		init_simple_evaluator(sphere);
		min_value = -100;
		max_value = 100;
	}
	else if (!strcmp(arguments[1], "ackley")) {
		init_simple_evaluator(ackley);
		min_value = -32;
		max_value = 32;
	} else if (!strcmp(arguments[1], "griewank")) {
		init_simple_evaluator(griewank);
		min_value = -600;
		max_value = 600;
	} else if (!strcmp(arguments[1], "rastrigin")) {
		init_simple_evaluator(rastrigin);
		min_value = -M_PI * 2;
		max_value = M_PI * 2;
	} else if (!strcmp(arguments[1], "rosenbrock")) {
		init_simple_evaluator(rosenbrock);
		min_value = -30;
		max_value = 30;
	} else {
		printf("search type specified: %s\n", arguments[1]);
		printf("unknown function type, valid functions: sphere, ackley, griewank, rastrigin, rosenbrock\n");
		exit(0);
	}

	number_parameters = atoi(arguments[2]);

	in_radians = (int*)malloc(sizeof(int) * number_parameters);
	min_bound = (double*)malloc(sizeof(double) * number_parameters);
	max_bound = (double*)malloc(sizeof(double) * number_parameters);
	range = (double*)malloc(sizeof(double) * number_parameters);
	point = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		in_radians[i] = 0;
		min_bound[i] = min_value;
		max_bound[i] = max_value;
		range[i] = 2.0;
		point[i] = (drand48() * (double)(max_bound[i] - min_bound[i])) + (double)min_bound[i];
	}

	new_bounds(&bounds, number_parameters, min_bound, max_bound, in_radians);

	/********
		*       Start the search
	 ********/
	if (argument_exists("-asynch", number_arguments, arguments)) {
//		if (argument_exists("-nm", number_arguments, arguments))	register_search(get_asynchronous_newton_method());
		if (argument_exists("-gs", number_arguments, arguments))	register_search(get_asynchronous_genetic_search());
		else if (argument_exists("-de", number_arguments, arguments))	register_search(get_asynchronous_differential_evolution());
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
