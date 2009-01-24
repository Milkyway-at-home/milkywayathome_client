#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "gradient.h"
#include "hessian.h"
#include "newton_method.h"
#include "line_search.h"
#include "regression.h"
#include "recombination.h"


#include "../evaluation/evaluator.h"

#include "../util/matrix.h"
#include "../util/io_util.h"

void newton_method(int number_parameters, double *point, double *range, int iterations) {
	int i, j, k;
	double **hessian, *gradient, *step;
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

		newton_step(number_parameters, hessian, gradient, &step);
		for (j = 0; j < number_parameters; j++) point[j] -= step[j];
		current_fitness = evaluate(point);

		printf("iteration %d [fitness : point] -- %.15lf :", i, current_fitness);
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", point[j]);
		printf("\n");
		free(step);
	}
}

void randomized_newton_method(int number_parameters, double *point, double *range, int type, int evaluations_per_iteration, int iterations, int ls_evaluations, int ls_iterations) {
	double *y, **x;
	double **hessian, **hessian_error;
	double *gradient, *gradient_error;
	double c, c_error;
	double *x_min, *x_max, *step, *step_error;
	double *new_point;
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
		x_min[i] = point[i] - range[i];
		x_max[i] = point[i] + range[i];
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
			for (k = 0; k < number_parameters; k++) x[j][k] = x[j][k] - point[k];
		}

		parabolic_2d_regression_error(evaluations_per_iteration, number_parameters, x, y, hessian, hessian_error, gradient, gradient_error, &c, &c_error);
		if (type == 0) {
			newton_step(number_parameters, hessian, gradient, &step);
			newton_step(number_parameters, hessian_error, gradient_error, &step_error);
			for (k = 0; k < number_parameters; k++) {
				point[k] -= step[k];
				x_min[k] = point[k] - step_error[k];
				x_max[k] = point[k] + step_error[k];
			}
			free(step);
			free(step_error);
		} else if (type == 1) {
			newton_step(number_parameters, hessian, gradient, &step);
			randomized_line_search(number_parameters, point, step, ls_evaluations, ls_iterations, &new_point, &current_fitness);
			for (k = 0; k < number_parameters; k++) {
				point[k] = new_point[k];
				x_min[k] = point[k] - range[k];
				x_max[k] = point[k] + range[k];
			}
			free(new_point);
			free(step);
		}

		printf("iteration %d [fitness : point] -- %.15lf :", i, current_fitness);
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", point[j]);
		printf("\n");

	}
}
