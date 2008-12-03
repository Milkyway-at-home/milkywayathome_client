#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "gradient.h"
#include "hessian.h"
#include "line_search.h"
#include "regression.h"
#include "recombination.h"

#include "../evaluation/evaluator.h"

#include "../util/matrix.h"
#include "../util/io_util.h"

void error_newton_range(double *point, double *min_step, double *max_step, int number_parameters, double *gradient, double *gradient_error, double **hessian, double **hessian_error) {
	int i;
	double **inverse_hessian, **inverse_hessian_error;
	double *step, *step_error;

	matrix_invert(hessian, number_parameters, number_parameters, &inverse_hessian);
//	matrix_print(stdout, "inverse hessian", inverse_hessian, number_parameters, number_parameters);

	vector_matrix_multiply(gradient, number_parameters, inverse_hessian, number_parameters, number_parameters, &step);
//	print_double_array(stdout, "step", number_parameters, step);

	matrix_invert(hessian_error, number_parameters, number_parameters, &inverse_hessian_error);
//	matrix_print(stdout, "hessian error", hessian_error, number_parameters, number_parameters);
//	matrix_print(stdout, "inverse hessian error", inverse_hessian_error, number_parameters, number_parameters);

	vector_matrix_multiply(gradient_error, number_parameters, inverse_hessian_error, number_parameters, number_parameters, &step_error);
	print_double_array(stdout, "step error", number_parameters, step_error);

	for (i = 0; i < number_parameters; i++) {
		point[i] = point[i] - step[i];
		min_step[i] = point[i] - step_error[i];
		max_step[i] = point[i] + step_error[i];
	}
	print_double_array(stdout, "min step", number_parameters, min_step);
	print_double_array(stdout, "max step", number_parameters, max_step);

	free(step);
	free(step_error);
	for (i = 0; i < number_parameters; i++) {
		free(inverse_hessian[i]);
		free(inverse_hessian_error[i]);
	}
	free(inverse_hessian);
	free(inverse_hessian_error);
}

void randomized_newton_step(double *point, double *point_fitness, int number_parameters, double *gradient, double *gradient_error, double **hessian, double **hessian_error, int evaluations) {
	double *step, *step_error, *min_step, *max_step, *current;
	double current_fitness;
	double **inverse_hessian, **inverse_hessian_error;
	int i, j;

	matrix_invert(hessian, number_parameters, number_parameters, &inverse_hessian);
//	matrix_print(stdout, "inverse hessian", inverse_hessian, number_parameters, number_parameters);

	vector_matrix_multiply(gradient, number_parameters, inverse_hessian, number_parameters, number_parameters, &step);
//	print_double_array(stdout, "step", number_parameters, step);

	matrix_invert(hessian_error, number_parameters, number_parameters, &inverse_hessian_error);
//	matrix_print(stdout, "hessian error", hessian_error, number_parameters, number_parameters);
//	matrix_print(stdout, "inverse hessian error", inverse_hessian_error, number_parameters, number_parameters);

	vector_matrix_multiply(gradient_error, number_parameters, inverse_hessian_error, number_parameters, number_parameters, &step_error);
	print_double_array(stdout, "step error", number_parameters, step_error);

	min_step = (double*)malloc(sizeof(double) * number_parameters);
	max_step = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		point[i] = point[i] - step[i];
		min_step[i] = point[i] - step_error[i];
		max_step[i] = point[i] + step_error[i];
	}
	print_double_array(stdout, "min step", number_parameters, min_step);
	print_double_array(stdout, "max step", number_parameters, max_step);

	(*point_fitness) = evaluate(point);
	current = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < evaluations; i++) {
		random_recombination(number_parameters, min_step, max_step, current);
		current_fitness = evaluate(current);

		if (current_fitness < (*point_fitness)) {
			(*point_fitness) = current_fitness;
			memcpy(point, current, sizeof(double) * number_parameters);
			printf("iteration %d, improved fitness to: %.15lf, new point:", i, (*point_fitness));
			for (j = 0; j < number_parameters; j++) printf(" %.15lf", point[j]);
			printf("\n");
		}
	}

	free(current);
	free(min_step);
	free(max_step);
	free(step);
	free(step_error);
	for (i = 0; i < number_parameters; i++) {
		free(inverse_hessian[i]);
		free(inverse_hessian_error[i]);
	}
	free(inverse_hessian);
	free(inverse_hessian_error);
}

void newton_step(double* point, double *point_fitness, int number_parameters, double *gradient, double **hessian) {
	double *step, *new_point;
	double **inverse_hessian;
	int i, evaluations, retval;

	matrix_invert(hessian, number_parameters, number_parameters, &inverse_hessian);
	matrix_print(stdout, "inverse hessian", inverse_hessian, number_parameters, number_parameters);

	vector_matrix_multiply(gradient, number_parameters, inverse_hessian, number_parameters, number_parameters, &step);
	print_double_array(stdout, "step", number_parameters, step);

	retval = line_search(point, (*point_fitness), step, number_parameters, &new_point, point_fitness, &evaluations);
	print_double_array(stdout, "\tnew point:", number_parameters, new_point);

	printf("\tline search took: %d evaluations for new fitness: %.15lf, with result: [%s]\n", evaluations, *(point_fitness), LS_STR[retval]);
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
			for (k = 0; k < number_parameters; k++) x[j][k] = x[j][k] - point[k];
		}

		parabolic2_regression(evaluations_per_iteration, number_parameters, x, y, hessian, hessian_error, gradient, gradient_error, &c, &c_error);

//		print_double_array(stdout, "gradient", number_parameters, gradient);
//		matrix_print(stdout, "hessian", hessian, number_parameters, number_parameters);

//		print_double_array(stdout, "gradient_error", number_parameters, gradient_error);
//		matrix_print(stdout, "hessian_error", hessian_error, number_parameters, number_parameters);

//		newton_step(point, &current_fitness, number_parameters, gradient, hessian);
		randomized_newton_step(point, &current_fitness, number_parameters, gradient, gradient_error, hessian, hessian_error, 100);

		printf("iteration %d [fitness : point] -- %.15lf :", i, current_fitness);
		for (j = 0; j < number_parameters; j++) printf(" %.15lf", point[j]);
		printf("\n");

		for (j = 0; j < number_parameters; j++) {
			x_min[j] = point[j] - step[j];
			x_max[j] = point[j] + step[j];
		}
	}
}
