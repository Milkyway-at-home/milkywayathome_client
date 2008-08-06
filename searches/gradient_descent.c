#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "gradient_descent.h"
#include "gradient.h"
#include "line_search.h"
#include "../evaluation/evaluator.h"
#include "../util/io_util.h"

void parse_gradient_parameters(char* parameters, int *number_iterations) {
	if (parameters[0] == 'g') {
		sscanf(parameters, "gd/%d", number_iterations);
	} else {
		sscanf(parameters, "cgd/%d", number_iterations);
	}
}

void synchronous_gradient_descent(char* search_path, char* search_parameters, double* point, double* step, int number_parameters) {
	GRADIENT *gradient;
	int i, evaluations, number_iterations;
	double point_fitness;
	double *new_point;

	parse_gradient_parameters(search_parameters, &number_iterations);

	point_fitness = evaluate(point);
	for (i = 0; i < number_iterations; i++) {
		printf("iteration %d:\n", i);
		synchronous_get_gradient(point, step, number_parameters, &gradient);

		evaluations = synchronous_line_search(point, point_fitness, gradient->values, number_parameters, &new_point, &point_fitness);
		print_double_array(stdout, "\tnew point:", number_parameters, new_point);
		printf("\tline search took: %d evaluations for new fitness: %lf\n", evaluations, point_fitness);

		free_gradient(gradient);
		free(gradient);

		if (evaluations < 0) break;
		memcpy(point, new_point, sizeof(double) * number_parameters);
	}
	free(new_point);
}

void synchronous_conjugate_gradient_descent(char* search_path, char* search_parameters, double* point, double* step, int number_parameters) {

}
