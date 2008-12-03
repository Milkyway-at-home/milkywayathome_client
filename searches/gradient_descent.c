#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "gradient_descent.h"
#include "gradient.h"
#include "line_search.h"
#include "../evaluation/evaluator.h"
#include "../util/io_util.h"

#define min_gradient_threshold 0.00001

void parse_gradient_parameters(char* parameters, int *number_iterations) {
	if (parameters[0] == 'g') {
		sscanf(parameters, "gd/%d", number_iterations);
	} else {
		sscanf(parameters, "cgd/%d", number_iterations);
	}
}

void synchronous_gradient_descent(char* search_path, char* search_parameters, double* point, double* step, int number_parameters) {
	GRADIENT *gradient;
	int i, evaluations, number_iterations, retval;
	double point_fitness;
	double *new_point;

	parse_gradient_parameters(search_parameters, &number_iterations);

	point_fitness = evaluate(point);
	for (i = 0; i < number_iterations; i++) {
		printf("iteration %d:\n", i);
		synchronous_get_gradient(point, step, number_parameters, &gradient);

		if (gradient_below_threshold(gradient, min_gradient_threshold)) {
			printf("Gradient dropped below threshold %.15lf\n", min_gradient_threshold);
			print_double_array(stdout, "\tgradient:", number_parameters, gradient->values);

			free_gradient(gradient);
			free(gradient);
			break;
		}

		retval = line_search(point, point_fitness, gradient->values, number_parameters, &new_point, &point_fitness, &evaluations);
		print_double_array(stdout, "\tnew point:", number_parameters, new_point);
		printf("\tline search took: %d evaluations for new fitness: %.15lf, with result: [%s]\n", evaluations, point_fitness, LS_STR[retval]);

		free_gradient(gradient);
		free(gradient);

		if (evaluations < 0) return;
		memcpy(point, new_point, sizeof(double) * number_parameters);
	}
	free(new_point);
}

void synchronous_conjugate_gradient_descent(char* search_path, char* search_parameters, double* point, double* step, int number_parameters) {
        GRADIENT *gradient;
	double *direction;
	double *previous_gradient;
	double *previous_direction;
        int i, j, evaluations, number_iterations, retval;
	int reset;
        double point_fitness;
        double *new_point;
	double bet, betdiv;

	reset = 8;
        parse_gradient_parameters(search_parameters, &number_iterations);

	previous_gradient = (double*)malloc(sizeof(double) * number_parameters);
	previous_direction = (double*)malloc(sizeof(double) * number_parameters);
	direction = (double*)malloc(sizeof(double) * number_parameters);
        point_fitness = evaluate(point);
        for (i = 0; i < number_iterations; i++) {
                printf("iteration %d:\n", i);
                synchronous_get_gradient(point, step, number_parameters, &gradient);
		if (gradient_below_threshold(gradient, min_gradient_threshold)) {
			printf("Gradient dropped below threshold %.15lf\n", min_gradient_threshold);
			print_double_array(stdout, "\tgradient:", number_parameters, gradient->values);

			free_gradient(gradient);
			free(gradient);
			break;
		}

		if (i > 0 && (i % reset) != 0) {
			// bet = gpres' * (gpres - gprev) / (gprev' * g_prev);
			bet = 0;
			betdiv = 0;
			for (j = 0; j < number_parameters; j++) {
				bet += (gradient->values[j] - previous_gradient[j]) * gradient->values[j];
				betdiv += previous_gradient[j] * previous_gradient[j];
			}
			bet /= betdiv;

			// dpres = -g_pres + bet * d_prev;
			for (j = 0; j < number_parameters; j++) {
				direction[j] = gradient->values[j] + bet * previous_direction[j];
			}
		} else {
			memcpy(direction, gradient->values, sizeof(double) * number_parameters);
		}
		memcpy(previous_direction, direction, sizeof(double) * number_parameters);
		memcpy(previous_gradient, gradient->values, sizeof(double) * number_parameters);

		print_double_array(stdout, "\tconjugate direction: ", number_parameters, direction);

		retval = line_search(point, point_fitness, gradient->values, number_parameters, &new_point, &point_fitness, &evaluations);
		print_double_array(stdout, "\tnew point:", number_parameters, new_point);
		printf("\tline search took: %d evaluations for new fitness: %.15lf, with result: [%s]\n", evaluations, point_fitness, LS_STR[retval]);

                free_gradient(gradient);
                free(gradient);

                if (evaluations < 0) break;
                memcpy(point, new_point, sizeof(double) * number_parameters);
        }
        free(new_point);

}
