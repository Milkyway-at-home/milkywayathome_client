#include "stdio.h"
#include "stdlib.h"

#include "gradient.h"
#include "hessian.h"
#include "line_search.h"
#include "regression.h"
#include "recombination.h"

#include "../evaluation/evaluator.h"

#include "../util/matrix.h"
#include "../util/io_util.h"


void newton_step(double* point, double *point_fitness, int number_parameters, double *gradient, double **hessian) {
	double *step, *new_point;
	double **inverse_hessian;
	int i, j, evaluations;

	matrix_invert(hessian, number_parameters, number_parameters, &inverse_hessian);
	printf("inverse hessian:\n ");
	for (i = 0; i < number_parameters; i++) {
		for (j = 0; j < number_parameters; j++) {
			printf(" %.15lf", inverse_hessian[i][j]);
		}
		printf("\n");
	}

	vector_matrix_multiply(gradient, number_parameters, inverse_hessian, number_parameters, number_parameters, &step);

	printf("step:");
	for (i = 0; i < number_parameters; i++) printf(" %.15lf", step[i]);
	printf("\n");

	evaluations = synchronous_line_search(point, (*point_fitness), step, number_parameters, &new_point, point_fitness);
	print_double_array(stdout, "\tnew point:", number_parameters, new_point);
	printf("\tline search took: %d evaluations for new fitness: %.15lf\n", evaluations, *(point_fitness));
	if (evaluations == -1) {
		printf("search completed\n");
		exit(0);
	}

	for (i = 0; i < number_parameters; i++) point[i] = new_point[i];
	free(new_point);
}

void newton_method(int number_parameters, double *point, double *step, int iterations) {
	int i, j, k;
	double **hessian, *gradient;
	double current_fitness;

	current_fitness = evaluate(point);

	printf("initial [fitness : point] -- %.15lf :", current_fitness);
	for (i = 0; i < number_parameters; i++) printf(" %.15lf", point[i]);
	printf("\n");

	gradient = (double*)malloc(sizeof(double) * number_parameters);
	hessian = (double**)malloc(sizeof(double*) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		hessian[i] = (double*)malloc(sizeof(double) * number_parameters);
	}

	for (i = 0; i < iterations; i++) {
		get_gradient(number_parameters, point, step, gradient);
		printf("gradient:");
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", gradient[j]);
		printf("\n");

		get_hessian(number_parameters, point, step, hessian);
		printf("hessian:\n");
		for (j = 0; j < number_parameters; j++) {
			printf(" ");
			for (k = 0; k < number_parameters; k++) printf(" %.15lf", hessian[j][k]);
			printf("\n");
		}

		newton_step(point, &current_fitness, number_parameters, gradient, hessian);

		printf("iteration %d [fitness : point] -- %.15lf :", i, current_fitness);
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", point[j]);
		printf("\n");
	}
}

void randomized_newton_method(int number_parameters, double *point, double *step, int evaluations_per_iteration, int iterations) {
	double *y, **x;
	double **hessian, **hessian_error;
	double *gradient, *gradient_error;
	double c, c_error;
	double *x_min, *x_max;
	double current_fitness;
	int i, j, k;

	current_fitness = evaluate(point);

	printf("initial [fitness : point] -- %.15lf :", current_fitness);
	for (i = 0; i < number_parameters; i++) printf(" %.15lf", point[i]);
	printf("\n");

	x_min = (double*)malloc(sizeof(double) * number_parameters);
	x_max = (double*)malloc(sizeof(double) * number_parameters);
	gradient = (double*)malloc(sizeof(double) * number_parameters);
	gradient_error = (double*)malloc(sizeof(double) * number_parameters);
	hessian = (double**)malloc(sizeof(double*) * number_parameters);
	hessian_error = (double**)malloc(sizeof(double*) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		hessian[i] = (double*)malloc(sizeof(double) * number_parameters);
		hessian_error[i] = (double*)malloc(sizeof(double) * number_parameters);
		x_min[i] = point[i] - step[i];
		x_max[i] = point[i] + step[i];
	}

	y = (double*)malloc(sizeof(double) * evaluations_per_iteration);
	x = (double**)malloc(sizeof(double*) * evaluations_per_iteration);
	for (i = 0; i < evaluations_per_iteration; i++) {
		x[i] = (double*)malloc(sizeof(double) * number_parameters);
	}

	for (i = 0; i < iterations; i++) {
		for (j = 0; j < evaluations_per_iteration; j++) {
			random_recombination(number_parameters, x_min, x_max, x[j]);
			y[j] = evaluate(x[j]);
		}

		parabolic2_regression(evaluations_per_iteration, number_parameters, x, y, hessian, hessian_error, gradient, gradient_error, &c, &c_error);
		printf("gradient:");
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", gradient[j]);
		printf("\n");

		printf("hessian:\n");
		for (j = 0; j < number_parameters; j++) {
			printf(" ");
			for (k = 0; k < number_parameters; k++) printf(" %.15lf", hessian[j][k]);
			printf("\n");
		}

		newton_step(point, &current_fitness, number_parameters, gradient, hessian);

		printf("iteration %d [fitness : point] -- %.15lf :", i, current_fitness);
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", point[j]);
		printf("\n");
	}
}
