/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

/********
	*	Includes for astronomy
 ********/
#include "mpi.h"
#include "math.h"

#include "astronomy_worker.h"
#include "parameters.h"
#include "star_points.h"
#include "evaluation.h"

#include "../evaluation/mpi_evaluator.h"
#include "../evaluation/evaluator.h"
#include "../searches/hessian.h"
#include "../util/matrix.h"

void errors(char* filename, double* point, double* step, int number_parameters, int number_stars) {
	double **hessian;
        double **sigma;
        int i, j;
        FILE *file;

        file = fopen(filename, "w");

        printf("Calculating hessian...\n");
	new_matrix(&hessian, number_parameters, number_parameters);
	get_hessian(number_parameters, point, step, hessian);

	fwrite_matrix(stdout, "hessian", hessian, number_parameters, number_parameters);
	fwrite_matrix(file, "hessian", hessian, number_parameters, number_parameters);
        fprintf(file, "\n");

        printf("\nInverting hessian to get SIGMA...\n");
	printf("Number of stars: %d\n", number_stars);
        matrix_invert__alloc(hessian, number_parameters, number_parameters, &sigma);
        printf("SIGMA:\n");
	fprintf(file, "SIGMA:\n");
        for (i = 0; i < number_parameters; i++) {
                for (j = 0; j < number_parameters; j++) {
                        sigma[i][j] /= (double)number_stars;
                        printf(" %lf", sigma[i][j]);
                        fprintf(file, " %lf", sigma[i][j]);
                }
                printf("\n");
		fprintf(file, "\n");
        }
	printf("\nCalculating standard deviations (errors)...\n");
	double stdDevs[number_parameters];
	printf("\nErrors: q, r0, epsilon, mu, r, theta, phi, sigma\n");
	fprintf(file,"\nErrors: q, r0, epsilon, mu, r, theta, phi, sigma\n");
	for (i = 0; i < number_parameters; i++) {
		stdDevs[i] = sqrt(fabs(sigma[i][i]));
		printf(" %lf", stdDevs[i]);
                fprintf(file, " %lf", stdDevs[i]);
	}
	printf("\n");
	fprintf(file, "\n");

}

int main(int number_arguments, char **arguments){
	double *point, *step;
	int number_parameters, integral_parameter_length, integral_results_length, likelihood_parameter_length, likelihood_results_length;
	int i;

	/********
		*	Initialize the mpi_evaluator
	 ********/
	mpi_evaluator__init(&number_arguments, &arguments);

	sprintf(star_points_file, "stars.txt");
	sprintf(astronomy_parameters_file, "parameters.txt");
	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-stars")) {
			sprintf(star_points_file, "%s", arguments[++i]);
		} else if (!strcmp(arguments[i], "-parameters")) {
			sprintf(astronomy_parameters_file, "%s", arguments[++i]);
		}
	}

	mpi_evaluator__read_data(read_data);

	number_parameters = get_optimized_parameter_count(ap);
	integral_parameter_length = number_parameters;
	integral_results_length = 1 + ap->number_streams;
	likelihood_parameter_length = 1 + ap->number_streams;
	likelihood_results_length = 2;

	printf("number_parameters: %d\n", number_parameters);
	printf("integral_parameter_length: %d, integral_results_length: %d\n", integral_parameter_length, integral_results_length);
	printf("likelihood_parameter_length: %d, likelihood_results_length: %d\n", likelihood_parameter_length, likelihood_results_length);

	evaluator__init_integral(integral_f, integral_parameter_length, integral_compose, integral_results_length);
	evaluator__init_likelihood(likelihood_f, likelihood_parameter_length, likelihood_compose, likelihood_results_length);

	/********
		*	Start the mpi_evaluator
	 ********/

	printf("starting...\n");
	mpi_evaluator__start();

	printf("getting parameters...\n");
	get_search_parameters(ap, &point);
	get_step(ap, &step);

      	errors(arguments[1], point, step, get_optimized_parameter_count(ap), sp->number_stars);
	return 0;
}
