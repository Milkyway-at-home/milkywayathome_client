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
	*	Includes for FGDO
 ********/
//#include "../evaluation/mpi_evaluator.h"

/********
	*	Includes for astronomy
 ********/
#include "astronomy_worker.h"
#include "parameters.h"
#include "star_points.h"
#include "evaluation_optimized.h"

char astronomy_parameters_file[1024];
char star_points_file[1024];

ASTRONOMY_PARAMETERS *ap;
STAR_POINTS *sp;
EVALUATION_STATE *es;
int total_number_stars;

void read_data(int rank, int max_rank) {
	/********
		*	READ THE ASTRONOMY PARAMETERS
	 ********/
	ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	int retval = read_astronomy_parameters(astronomy_parameters_file, ap);
	if (retval) {
		fprintf(stderr, "APP: error reading astronomy parameters: %d\n", retval);
		exit(1);
	}
	split_astronomy_parameters(ap, rank, max_rank);

	/********
		*	READ THE STAR POINTS
	 ********/
	sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
	retval = read_star_points(star_points_file, sp);
	if (retval) {
		fprintf(stderr, "APP: error reading star points: %d\n", retval);
		exit(1);
	}
	total_number_stars = sp->number_stars;
	split_star_points(sp, rank, max_rank);

	/********
		*	INITIALIZE THE EVALUATION STATE
	 ********/
	es = (EVALUATION_STATE*)malloc(sizeof(EVALUATION_STATE));
	initialize_state(ap, sp, es);
}

void integral_f(double* parameters, double* results) {
	int i;
	/********
		*	CALCULATE THE INTEGRALS
	 ********/
	set_astronomy_parameters(ap, parameters);
	reset_evaluation_state(es);

	int retval = calculate_integrals(ap, es, sp);
	if (retval) {
		fprintf(stderr, "APP: error calculating integrals: %d\n", retval);
		exit(retval);
	}

	results[0] = es->background_integral;
	for (i = 0; i < ap->number_streams; i++) {
		results[i+1] = es->stream_integrals[i];
	}
//	printf("[worker: %d] background integral: %lf, stream integrals:", get_mpi_rank(), results[0]);
//	for (i = 0; i < ap->number_streams; i++) printf(" %lf", results[i+1]);
//	printf("\n");
}

void integral_compose(double* integral_results, int num_results, double* results) {
	int i, j, current;
	results[0] = 0.0;
	for (i = 0; i < ap->number_streams; i++) {
		results[i+1] = 0.0;
	}

	for (i = 0; i < num_results; i++) {
		current = (ap->number_streams + 1) * i;
		results[0] += integral_results[current];
		for (j = 0; j < ap->number_streams; j++) {
			results[j+1] += integral_results[current + j + 1];
		}
	}
//	printf("[worker: %d] [compose] background integral: %lf, stream integrals:", get_mpi_rank(), results[0]);
//	for (i = 0; i < ap->number_streams; i++) printf(" %lf", results[i+1]);
//	printf("\n");
}

void likelihood_f(double* integrals, double* results) {
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
	results[0] = es->prob_sum;
	results[1] = es->bad_jacobians;
//	printf("[worker: %d] calculated likelihood: %lf, bad_jacobs: %lf\n", get_mpi_rank(), results[0], results[1]);
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
//	printf("[worker: %d] [compose] likelihood: %.10lf\n", get_mpi_rank(), prob_sum);
	return prob_sum - 3.0;
}
