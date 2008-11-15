#include "stdio.h"
#include "stdlib.h"

#include "hessian.h"
#include "gradient.h"

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
		point[i] = point[i] - step[i];
	}
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

		newton_step(point, gradient, hessian);

		current_fitness = evaluate(point);

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
			p = new_population(search_path, search_parameters, min_parameters, max_parameters, number_parameters, population_size, max_evaluations);
		} else {
			reset_population(p, min_parameters, max_parameters);
		}

		printf("\nevaluating individuals (fitness : parameters):\n");
		for (j = 0; j < population_size; j++) {
			parameters = random_recombination(min_parameters, max_parameters, number_parameters);
			fitness = evaluate(parameters);
			replace(p, j, parameters, fitness);

			printf("\t%lf :", fitness);
			for (k = 0; k < number_parameters; k++) printf(" %lf", parameters[k]);
			printf("\n");

			free(parameters);
		}

		randomized_hessian(p->individuals, p->fitness, population_size, number_parameters, &hessian, &gradient);

		printf("\n");
		print_double_array(stdout, "gradient: ", number_parameters, gradient);
		printf("\n");
		matrix_print(stdout, "hessian", hessian, number_parameters, number_parameters);
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


typedef struct newton_method_search {
	int current_iteration, maximum_iteration;
	int current_evaluation, evaluations_per_iteration;
	double *current_parameters;
	double *parameter_range;

	double best_fitness;
	double *best_parameters;

	POPULATION *current_population;
} NEWTON_METHOD_SEARCH;

int newton_generate_parameters(SEARCH* search, SEARCH_PARAMETERS* sp) {
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)(search->search_data);
	POPULATION *p = nms->current_population;

	if (current_iteration <= maximum_iteration) {
		char metadata[METADATA_SIZE];
		sprintf(metadata, "", );
		new_search_parameters(&sp, search->search_name, p->number_parameters, random_recombination(p->min_parameters, p->max_parameters, p->number_parameters), metadata);
	}
}

int newton_insert_parameters(SEARCH* search, SEARCH_PARAMETERS* sp) {
	NEWTON_METHOD_SEARCH *nms = (NEWTON_METHOD_SEARCH*)(search->search_data);

	/********
		*	Insert parameters into population.  If cutoff reached, calculate hessian
		*	and generate new population.
	 ********/
}

int init_newton_method(char* search_name, SEARCH* search) {



}
