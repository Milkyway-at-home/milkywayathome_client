#include "stdio.h"
#include "stdlib.h"

#include "hessian.h"
#include "gradient.h"

#include "../evaluation/evaluator.h"
#include "../util/matrix.h"

void parse_newton_parameters(char* parameters, int *number_iterations) {
	sscanf(parameters, "newton/%d", number_iterations);
}


void newton_step(double* point, GRADIENT* gradient, HESSIAN* hessian) {
	double* step;
	double** inverse_hessian;
	int i, j;

	printf("inverting\n");
	matrix_invert(hessian->values, hessian->number_parameters, hessian->number_parameters, &inverse_hessian);
	printf("inverse hessian:\n");
	for (i = 0; i < hessian->number_parameters; i++) {
		for (j = 0; j < hessian->number_parameters; j++) {
			printf(" %lf", inverse_hessian[i][j]);
		}
		printf("\n");
	}

	printf("multiplying\n");
	matrix_vector_multiply(gradient->values, hessian->number_parameters, inverse_hessian, hessian->number_parameters, hessian->number_parameters, &step);

	printf("step:");
	for (i = 0; i < hessian->number_parameters; i++) printf(" %lf", step[i]);
	printf("\n");

	printf("updating point\n");
	for (i = 0; i < hessian->number_parameters; i++) {
		point[i] = point[i] + step[i];
	}
}

void synchronous_newton_method(char* search_path, char* search_parameters, double* point, double* step, int number_parameters) {
	HESSIAN* hessian;
	GRADIENT* gradient;
	char **metadata;
	double **individuals;
	double *fitness;
	int number_individuals;
	int number_iterations;
	int i, j;

	number_individuals = 1;

	printf("iteration: %d, current point:", 0);
	for (j = 0; j < number_parameters; j++) {
		printf(" %lf", point[j]);
	}
	printf("\n");

	parse_newton_parameters(search_parameters, &number_iterations);

	for (i = 0; i < number_iterations; i++) {
		create_hessian(point, step, i, number_parameters, &hessian);

		while (!hessian__complete(hessian)) {
			hessian__get_individuals(hessian, number_individuals, &individuals, &metadata);
			fitness = (double*)malloc(sizeof(double) * number_individuals);
			for (j = 0; j < number_individuals; j++) {
				fitness[j] = evaluate(individuals[j]);
			}
			hessian__insert_individuals(hessian, number_individuals, fitness, metadata);

			for (j = 0; j < number_individuals; j++) {
				free(individuals[j]);
				free(metadata[j]);
			}
			free(individuals);
			free(metadata);
			free(fitness);
		}
		fprintf_hessian(stdout, hessian);

		create_gradient(point, step, i, number_parameters, &gradient);
		while (!gradient__complete(gradient)) {
			gradient__get_individuals(gradient, number_individuals, &individuals, &metadata);
			fitness = (double*)malloc(sizeof(double) * number_individuals);
			for (j = 0; j < number_individuals; j++) {
				fitness[j] = evaluate(individuals[j]);
			}
			gradient__insert_individuals(gradient, number_individuals, fitness, metadata);

			for (j = 0; j < number_individuals; j++) {
				free(individuals[j]);
				free(metadata[j]);
			}
			free(individuals);
			free(metadata);
			free(fitness);
		}
		fprintf_gradient(stdout, gradient);

		newton_step(point, gradient, hessian);

		printf("iteration: %d, current point:", i);
		for (j = 0; j < number_parameters; j++) {
			printf(" %lf", point[j]);
		}
		printf("\n");

		free_gradient(gradient);
		free_hessian(hessian);
		free(gradient);
		free(hessian);
	}
}
