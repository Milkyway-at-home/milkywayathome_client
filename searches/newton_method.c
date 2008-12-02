#include "stdio.h"
#include "stdlib.h"

#include "hessian.h"
#include "gradient.h"

#include "../searches/line_search.h"
#include "../searches/population.h"
#include "../searches/recombination.h"
#include "../evaluation/evaluator.h"
#include "../util/matrix.h"
#include "../util/io_util.h"

void parse_newton_parameters(char* parameters, int *number_iterations) {
	sscanf(parameters, "newton/%d", number_iterations);
}

void parse_randomized_newton(char* parameters, int *population_size, int *max_evaluations, int *number_iterations) {
	sscanf(parameters, "randomized_newton/%d/%d/%d", population_size, max_evaluations, number_iterations);
}

void newton_step(double* point, double *point_fitness, GRADIENT* gradient, HESSIAN* hessian) {
	double *step, *new_point;
	double **inverse_hessian;
	int i, j, evaluations;

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

	printf("updating point with line search\n");
	evaluations = synchronous_line_search(point, (*point_fitness), step, hessian->number_parameters, &new_point, point_fitness);
	print_double_array(stdout, "\tnew point:", hessian->number_parameters, new_point);
	printf("\tline search took: %d evaluations for new fitness: %lf\n", evaluations, *(point_fitness));
	for (i = 0; i < hessian->number_parameters; i++) point[i] = new_point[i];
	free(new_point);
}

void synchronous_newton_method(char* search_path, char* search_parameters, double* point, double* step, int number_parameters) {
	HESSIAN* hessian;
	GRADIENT* gradient;
	char **metadata;
	double **individuals;
	double *fitness;
	double current_fitness;
	int number_individuals;
	int number_iterations;
	int i, j;

	number_individuals = 1;

	current_fitness = evaluate(point);
	printf("iteration: %d, fitness: %lf, current point:", 0, current_fitness);
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

		newton_step(point, &current_fitness, gradient, hessian);
//		current_fitness = evaluate(point);

		printf("iteration: %d, fitness: %lf, current point:", i, current_fitness);
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

void randomized_newton_method(char* search_path, char* search_parameters, double* point, double* step, int number_parameters) {
	POPULATION* p;
	double *current_point, *min_parameters, *max_parameters;
	double* parameters;
	double fitness;
	double** hessian;
	double* gradient;
	double** inverse_hessian;
	int population_size, max_evaluations, number_iterations;
	int i, j, k;

	parse_randomized_newton(search_parameters, &population_size, &max_evaluations, &number_iterations);

	current_point = (double*)malloc(sizeof(double) * number_parameters);
	min_parameters = (double*)malloc(sizeof(double) * number_parameters);
	max_parameters = (double*)malloc(sizeof(double) * number_parameters);
	parameters = (double*)malloc(sizeof(double) * number_parameters);

	for (i = 0; i < number_parameters; i++) current_point[i] = point[i];

	p = NULL;
	for (i = 0; i < number_iterations; i++) {
		for (j = 0; j < number_parameters; j++) {
			min_parameters[j] = current_point[j] - step[j];
			max_parameters[j] = current_point[j] + step[j];
		}
		print_double_array(stdout, "current point: ", number_parameters, current_point);
		print_double_array(stdout, "min range: ", number_parameters, min_parameters);
		print_double_array(stdout, "max range: ", number_parameters, max_parameters);
		
		if (p == NULL) {
			new_population(max_evaluations, number_parameters, &p);
		}

//		printf("\nevaluating individuals (fitness : parameters):\n");
		for (j = 0; j < population_size; j++) {
			random_recombination(min_parameters, max_parameters, number_parameters, parameters);
			fitness = evaluate(parameters);
			replace(p, j, parameters, fitness);

//			printf("\t%lf :", fitness);
//			for (k = 0; k < number_parameters; k++) printf(" %lf", parameters[k]);
//			printf("\n");
		}

		randomized_hessian(p->individuals, current_point, p->fitness, population_size, number_parameters, &hessian, &gradient);

		printf("\n");
		matrix_print(stdout, "hessian", hessian, number_parameters, number_parameters);
		printf("\n");
		print_double_array(stdout, "gradient: ", number_parameters, gradient);
		printf("\n");

		/********
			*	Take the newton step:  y = -hessian^-1 * gradient
		 ********/
		matrix_invert(hessian, number_parameters, number_parameters, &inverse_hessian);
		for (j = 0; j < number_parameters; j++) {
			current_point[j] = 0;
			for (k = 0; k < number_parameters; k++) current_point[j] += inverse_hessian[j][k] * gradient[j];
			current_point[j] = -current_point[j];
		}

		/********
			*	TODO
			*	Optional: use previous current_point with new point as a line search direction
		 ********/

		for (j = 0; j < number_parameters; j++) {
			free(inverse_hessian[j]);
			free(hessian[j]);
		}
		free(hessian);
		free(gradient);
	}
	print_double_array(stdout, "final point: ", number_parameters, current_point);

	free(current_point);
	free(min_parameters);
	free(max_parameters);
}
