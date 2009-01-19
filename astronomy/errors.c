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

#include "parameters.h"
#include "star_points.h"
#include "evaluation.h"

#include "../evaluation/mpi_evaluator.h"
#include "../evaluation/evaluator.h"
#include "../searches/hessian.h"
#include "../util/matrix.h"

#define max_iterations			35000
#define astronomy_parameters_file	"parameters.txt"
#define star_points_file		"stars_unconvolved_82.txt"
#define population_file_name		"population.txt"

ASTRONOMY_PARAMETERS *ap;
STAR_POINTS *sp;
EVALUATION_STATE *es;
int total_number_stars;

void read_data(int rank, int max_rank) {
	ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	/********
		*	READ THE ASTRONOMY PARAMETERS
	 ********/
	printf("reading parameters...\n");
	int retval = read_astronomy_parameters(astronomy_parameters_file, ap);
	if (retval) {
		fprintf(stderr, "APP: error reading astronomy parameters: %d\n", retval);
		exit(1);
	}
	fwrite_astronomy_parameters(stdout, ap);
	printf("splitting parameters...\n");
	split_astronomy_parameters(ap, rank, max_rank);
	/********
		*	READ THE STAR POINTS
	 ********/
	printf("reading star points...\n");
	sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
	retval = read_star_points(star_points_file, sp);
	if (retval) {
		fprintf(stderr, "APP: error reading star points: %d\n", retval);
		exit(1);
	}
	printf("read %d stars.\n", sp->number_stars);
	total_number_stars = sp->number_stars;
	split_star_points(sp, rank, max_rank);

	/********
		*	INITIALIZE THE EVALUATION STATE
	 ********/
	printf("initializing state...\n");
	es = (EVALUATION_STATE*)malloc(sizeof(EVALUATION_STATE));
	initialize_state(es, ap->number_streams);
}

void integral_f(double* parameters, double** results) {
	int i;
	/********
		*	CALCULATE THE INTEGRALS
	 ********/
	set_astronomy_parameters(ap, parameters);

	es->r_step_current = 0;
	es->mu_step_current = 0;
	es->nu_step_current = 0;

	es->background_integral = 0.0;
	for (i = 0; i < es->number_streams; i++) {
		es->stream_integrals[i] = 0.0;
	}
	es->current_star_point = 0;
	es->num_zero = 0;
	es->bad_jacobians = 0;
	es->prob_sum = 0.0;

	int retval = calculate_integrals(ap, es, sp);
	if (retval) {
		fprintf(stderr, "APP: error calculating integrals: %d\n", retval);
		exit(retval);
	}
	(*results) = (double*)malloc(sizeof(double) * (1 + ap->number_streams));
	(*results)[0] = es->background_integral;
	printf("background integral: %lf, stream integrals: ", (*results[0]));
	for (i = 1; i <= ap->number_streams; i++) {
		(*results)[i] = es->stream_integrals[i-1];
		printf(" %lf", (*results)[i]);
	}
	printf("\n");
}

void integral_compose(double* integral_results, int num_results, double** results) {
	int i, j;
	(*results) = (double*)malloc(sizeof(double) * (1 + ap->number_streams));
	(*results)[0] = 0.0;
	for (i = 0; i < ap->number_streams; i++) {
		(*results)[i] = 0.0;
	}

	for (i = 0; i < num_results; i++) {
		(*results)[0] += integral_results[((ap->number_streams+1)*i)];
		for (j = 0; j < ap->number_streams; j++) {
			(*results)[j+1] += integral_results[((ap->number_streams+1)*i)+j];
		}
	}
	printf("background integral: %lf, stream integrals:", (*results)[0]);
	for (i = 0; i < ap->number_streams; i++) printf(" %lf", (*results)[i+1]);
	printf("\n");
}

void likelihood_f(double* integrals, double** results) {
	int i;

	es->background_integral = integrals[0];
	for (i = 0; i < ap->number_streams; i++) es->stream_integrals[i] = integrals[i+1];

	/********
		*	CALCULATE THE LIKELIHOOD
	 ********/
	int retval = calculate_likelihood(ap, es, sp);
	if (retval) {
		fprintf(stderr, "APP: error calculating likelihood: %d\n", retval);
		exit(retval);
	}
	(*results) = (double*)malloc(sizeof(double) * 2);
	(*results)[0] = es->prob_sum;
	(*results)[1] = es->bad_jacobians;
//	printf("calculated likelihood: %lf, bad_jacobs: %lf\n", (*results)[0], (*results)[1]);
}

double likelihood_compose(double* results, int num_results) {
	double prob_sum = 0.0;
	double bad_jacobians = 0.0;
	int i;
	for (i = 0; i < num_results; i++) {
		prob_sum += results[(2*i)];
		bad_jacobians += results[(2*i)+1];
	}
	prob_sum /= (total_number_stars - bad_jacobians);
	prob_sum *= -1;
	printf("composed likelihood: %lf\n", prob_sum);
	return prob_sum;
}

void errors(char* filename, double* point, double* step, int number_parameters, int number_stars) {
        HESSIAN* hessian;
        char **metadata;
        double **individuals;
        double *fitness;
        double **sigma;
        int i, j;
        FILE *file;

        file = fopen(filename, "w");

        int number_individuals = 1;

int count = 0;
        printf("Calculating hessian...\n");
        create_hessian(point, step, 0, number_parameters, &hessian);
        while (!hessian__complete(hessian)) {
                hessian__get_individuals(hessian, number_individuals, &individuals, &metadata);
                fitness = (double*)malloc(sizeof(double) * number_individuals);
                for (j = 0; j < number_individuals; j++) {
                        fitness[j] = evaluate(individuals[j]);
			printf("fitness[%d/%d] = %lf\n", count, hessian->number_parameters*hessian->number_parameters*4, fitness[j]);
			count++;
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
        fprintf_hessian(file, hessian);

        fprintf(file, "\n");

        printf("\nInverting hessian to get SIGMA...\n");
	printf("Number of stars: %d\n", number_stars);
        matrix_invert(hessian->values, hessian->number_parameters, hessian->number_parameters, &sigma);
        printf("SIGMA:\n");
	fprintf(file, "SIGMA:\n");
        for (i = 0; i < hessian->number_parameters; i++) {
                for (j = 0; j < hessian->number_parameters; j++) {
                        sigma[i][j] /= (double)number_stars;
                        printf(" %lf", sigma[i][j]);
                        fprintf(file, " %lf", sigma[i][j]);
                }
                printf("\n");
		fprintf(file, "\n");
        }
	printf("\nCalculating standard deviations (errors)...\n");
	double stdDevs[hessian->number_parameters];
	printf("\nErrors: q, r0, epsilon, mu, r, theta, phi, sigma\n");
	fprintf(file,"\nErrors: q, r0, epsilon, mu, r, theta, phi, sigma\n");
	for (i = 0; i < hessian->number_parameters; i++) {
		stdDevs[i] = sqrt(sigma[i][i]);
		printf(" %lf", stdDevs[i]);
                fprintf(file, " %lf", stdDevs[i]);
	}
	printf("\n");
	fprintf(file, "\n");

}

int main(int number_arguments, char **arguments){
	int integral_parameter_length, integral_results_length;
	int likelihood_parameter_length, likelihood_results_length;
	double *point;
	double *step;

	printf("init data...\n");
	evaluator__init(&number_arguments, &arguments, read_data);

	integral_parameter_length = ap->number_parameters;
	integral_results_length = 1 + ap->number_streams;
	likelihood_parameter_length = 1 + ap->number_streams;
	likelihood_results_length = 2;

	printf("init integral...\n");
	evaluator__init_integral(integral_f, integral_parameter_length, integral_compose, integral_results_length);
	printf("init likelihood...\n");
	evaluator__init_likelihood(likelihood_f, likelihood_parameter_length, likelihood_compose, likelihood_results_length);
	printf("starting...\n");
	mpi_evaluator__start();

	printf("getting parameters...\n");
	get_search_parameters(ap, &point);
	get_step(ap, &step);

      	errors(arguments[1], point, step, ap->number_parameters, sp->number_stars);

	return 0;
}
