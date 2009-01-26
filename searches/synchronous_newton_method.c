#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "gradient.h"
#include "hessian.h"
#include "newton_method.h"
#include "line_search.h"
#include "regression.h"
#include "recombination.h"
#include "search_arguments.h"


#include "../evaluation/evaluator.h"

#include "../util/matrix.h"
#include "../util/io_util.h"

void synchronous_newton_method(int number_arguments, char** arguments, int number_parameters, double *point, double *range) {
	int i, j, k;
	double **hessian, *gradient, *step;
	double current_fitness;
	int iterations;

	iterations = get_int_arg("-nm_iterations", number_arguments, arguments);
	if (iterations <= 0) {
		printf("argument: '-nm_iterations #' not specified, quitting.\n");
		return;
	}

	current_fitness = evaluate(point);

	printf("initial [fitness : point] -- %.15lf :", current_fitness);
	for (i = 0; i < number_parameters; i++) printf(" %.15lf", point[i]);
	printf("\n");

	step = (double*)malloc(sizeof(double) * number_parameters);
	gradient = (double*)malloc(sizeof(double) * number_parameters);
	hessian = (double**)malloc(sizeof(double*) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		hessian[i] = (double*)malloc(sizeof(double) * number_parameters);
	}

	for (i = 0; i < iterations; i++) {
		get_gradient(number_parameters, point, range, gradient);
		printf("gradient:");
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", gradient[j]);
		printf("\n");

		get_hessian(number_parameters, point, range, hessian);
		printf("hessian:\n");
		for (j = 0; j < number_parameters; j++) {
			printf(" ");
			for (k = 0; k < number_parameters; k++) printf(" %.15lf", hessian[j][k]);
			printf("\n");
		}

		newton_step(number_parameters, hessian, gradient, step);
		for (j = 0; j < number_parameters; j++) point[j] -= step[j];
		current_fitness = evaluate(point);

		printf("iteration %d [fitness : point] -- %.15lf :", i, current_fitness);
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", point[j]);
		printf("\n");
	}
}
