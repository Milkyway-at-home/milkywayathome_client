/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "synchronous_gradient_descent.h"
#include "gradient.h"
#include "line_search.h"
#include "search_arguments.h"
#include "../evaluation/evaluator.h"
#include "../util/io_util.h"

void synchronous_gradient_descent(int number_arguments, char **arguments, int number_parameters, double *point, double *step) {
	int i, evaluations, number_iterations, retval;
	double point_fitness;
	double *new_point, *gradient;
	double min_gradient_threshold;

	number_iterations = get_int_arg("-gd_iterations", number_arguments, arguments);
	if (number_iterations <= 0) {
		printf("argument: '-gd_iterations #' not specified, quitting.\n");
		return;
	}

	min_gradient_threshold = get_double_arg("-gd_min_threshold", number_arguments, arguments);
	if (min_gradient_threshold < 0) {
		printf("argument: '-gd_min_threshold #' not specified, quitting.\n");
		return;
	}

	new_point = (double*)malloc(sizeof(double) * number_parameters);
	gradient = (double*)malloc(sizeof(double) * number_parameters);

	point_fitness = evaluate(point);
	for (i = 0; i < number_iterations; i++) {
                printf("iteration %d: current_fitness: %.20lf\n", i, point_fitness);
		get_gradient(number_parameters, point, step, gradient);

		if (gradient_below_threshold(number_parameters, gradient, min_gradient_threshold)) {
			printf("Gradient dropped below threshold %.15lf\n", min_gradient_threshold);
			fwrite_double_array(stdout, "\tgradient:", number_parameters, gradient);
			break;
		}

		retval = line_search(point, point_fitness, gradient, number_parameters, new_point, &point_fitness, &evaluations);
		fwrite_double_array(stdout, "\tnew point:", number_parameters, new_point);
		printf("\tline search took: %d evaluations for new fitness: %.15lf, with result: [%s]\n", evaluations, point_fitness, LS_STR[retval]);

		if (retval != LS_SUCCESS) break;
		if (evaluations < 0) break;
		memcpy(point, new_point, sizeof(double) * number_parameters);
	}
	free(new_point);
	free(gradient);
}

void synchronous_conjugate_gradient_descent(int number_arguments, char **arguments, int number_parameters, double *point, double *step) {
        double *gradient;
	double *direction;
	double *previous_gradient;
	double *previous_direction;
        int i, j, evaluations, number_iterations, retval;
	int reset, reset_count;
        double point_fitness;
        double *new_point;
	double bet, betdiv;
	double min_gradient_threshold;

	number_iterations = get_int_arg("-gd_iterations", number_arguments, arguments);
	if (number_iterations <= 0) {
		printf("argument: '-gd_iterations #' not specified, quitting.\n");
		return;
	}

	reset = 0;
	reset_count = get_int_arg("-gd_reset", number_arguments, arguments);
	if (reset < 0) {
		printf("argument: '-gd_reset #' not specified, quitting.\n");
		return;
	}

	min_gradient_threshold = get_double_arg("-gd_min_threshold", number_arguments, arguments);
	if (min_gradient_threshold < 0) {
		printf("argument: '-gd_min_threshold #' not specified, quitting.\n");
		return;
	}

	previous_gradient = (double*)malloc(sizeof(double) * number_parameters);
	previous_direction = (double*)malloc(sizeof(double) * number_parameters);
	direction = (double*)malloc(sizeof(double) * number_parameters);
	gradient = (double*)malloc(sizeof(double) * number_parameters);
	new_point = (double*)malloc(sizeof(double) * number_parameters);
        point_fitness = evaluate(point);

        for (i = 0; i < number_iterations; i++) {
                printf("iteration %d: current_fitness: %.20lf\n", i, point_fitness);
                get_gradient(number_parameters, point, step, gradient);
		if (gradient_below_threshold(number_parameters, gradient, min_gradient_threshold)) {
			printf("Gradient dropped below threshold %.15lf\n", min_gradient_threshold);
			fwrite_double_array(stdout, "\tgradient:", number_parameters, gradient);

			break;
		}

		if (i > 0 && reset != 0) {
			// bet = g_pres' * (g_pres - g_prev) / (g_prev' * g_prev);
			bet = 0;
			betdiv = 0;
			for (j = 0; j < number_parameters; j++) {
				bet += (gradient[j] - previous_gradient[j]) * gradient[j];
				betdiv += previous_gradient[j] * previous_gradient[j];
			}
			bet /= betdiv;

			// dpres = -g_pres + bet * d_prev;
			for (j = 0; j < number_parameters; j++) {
				direction[j] = gradient[j] + bet * previous_direction[j];
			}
		} else {
			memcpy(direction, gradient, sizeof(double) * number_parameters);
		}
		memcpy(previous_direction, direction, sizeof(double) * number_parameters);
		memcpy(previous_gradient, gradient, sizeof(double) * number_parameters);

		fwrite_double_array(stdout, "\tconjugate direction: ", number_parameters, direction);

		retval = line_search(point, point_fitness, direction, number_parameters, new_point, &point_fitness, &evaluations);
		fwrite_double_array(stdout, "\tnew point:", number_parameters, new_point);
		printf("\tline search took: %d evaluations for new fitness: %.15lf, with result: [%s]\n", evaluations, point_fitness, LS_STR[retval]);

		if ((retval != LS_SUCCESS || evaluations < 0) && reset == 0) break;
                memcpy(point, new_point, sizeof(double) * number_parameters);
		reset++;
		if (reset == reset_count) reset = 0;
        }
        free(new_point);
	free(gradient);
}

