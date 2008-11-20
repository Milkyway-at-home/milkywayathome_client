#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/********
	*	FGDO includes
 ********/

/********
	*	Includes for astronomy
 ********/
#include "mpi.h"

#include "parameters.h"
#include "star_points.h"
#include "evaluation_optimized.h"

#include "../evaluation/mpi_evaluator.h"
#include "../evaluation/evaluator.h"
#include "../evaluation/mpi_search_manager.h"
#include "../evaluation/search_manager.h"
#include "../searches/asynchronous_search.h"
#include "../searches/asynchronous_newton_method.h"
//#include "../searches/asynchronous_genetic_search.h"
#include "../searches/search_parameters.h"

#include "../util/io_util.h"

#define max_iterations			35000
#define astronomy_parameters_file	"parameters.txt"
#define star_points_file		"86-cut.txt"
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
	printf("[worker: %d] reading parameters...\n", rank);
	int retval = read_astronomy_parameters(astronomy_parameters_file, ap);
	if (retval) {
		fprintf(stderr, "APP: error reading astronomy parameters: %d\n", retval);
		exit(1);
	}
//	fwrite_astronomy_parameters(stdout, ap);
	printf("[worker: %d] splitting parameters...\n", rank);
	split_astronomy_parameters(ap, rank, max_rank);
	/********
		*	READ THE STAR POINTS
	 ********/
	printf("[worker: %d] reading star points...\n", rank);
	sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
	retval = read_star_points(star_points_file, sp);
	if (retval) {
		fprintf(stderr, "APP: error reading star points: %d\n", retval);
		exit(1);
	}
	printf("[worker: %d] read %d stars.\n", rank, sp->number_stars);
	total_number_stars = sp->number_stars;
	split_star_points(sp, rank, max_rank);

	/********
		*	INITIALIZE THE EVALUATION STATE
	 ********/
	printf("[worker: %d] initializing state...\n", rank);
	es = (EVALUATION_STATE*)malloc(sizeof(EVALUATION_STATE));
	initialize_state(ap, es);
}

void integral_f(double* parameters, double** results) {
	int i;
	/********
		*	CALCULATE THE INTEGRALS
	 ********/
//	printf("[worker: %d] ", get_mpi_rank());
//	print_double_array(stdout, "parameters", ap->number_parameters, parameters);
	set_astronomy_parameters(ap, parameters);
	reset_evaluation_state(ap, es);

	int retval = calculate_integrals(ap, es, sp);
	if (retval) {
		fprintf(stderr, "APP: error calculating integrals: %d\n", retval);
		exit(retval);
	}

	(*results) = (double*)malloc(sizeof(double) * (1 + ap->number_streams));
	(*results)[0] = es->background_integral;
	for (i = 0; i < ap->number_streams; i++) {
		(*results)[i+1] = es->stream_integrals[i];
	}
//	printf("[worker] background integral: %lf, stream integrals:", (*results)[0]);
//	for (i = 0; i < ap->number_streams; i++) printf(" %lf", (*results)[i+1]);
//	printf("\n");
}

void integral_compose(double* integral_results, int num_results, double** results) {
	int i, j, current;
	(*results) = (double*)malloc(sizeof(double) * (1 + ap->number_streams));
	(*results)[0] = 0.0;
	for (i = 0; i < ap->number_streams; i++) {
		(*results)[i+1] = 0.0;
	}

	for (i = 0; i < num_results; i++) {
		current = (ap->number_streams + 1) * i;
		(*results)[0] += integral_results[current];
		for (j = 0; j < ap->number_streams; j++) {
			(*results)[j+1] += integral_results[current + j + 1];
		}
	}
//	printf("[compose] background integral: %lf, stream integrals:", (*results)[0]);
//	for (i = 0; i < ap->number_streams; i++) printf(" %lf", (*results)[i+1]);
//	printf("\n");
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
	printf("[worker: %d] composed likelihood: %.10lf\n", get_mpi_rank(), prob_sum);
	return prob_sum;
}

int main(int argc, char **argv) {
	int i;
	int integral_parameter_length, integral_results_length;
	int likelihood_parameter_length, likelihood_results_length;
	double *min_parameters;
	double *max_parameters;
	double *point;
	double *step;

	printf("initializing mpi evaluator\n");
	evaluator__init(&argc, &argv, read_data);

	integral_parameter_length = ap->number_parameters;
	integral_results_length = 1 + ap->number_streams;
	likelihood_parameter_length = 1 + ap->number_streams;
	likelihood_results_length = 2;

	evaluator__init_integral(integral_f, integral_parameter_length, integral_compose, integral_results_length);
	evaluator__init_likelihood(likelihood_f, likelihood_parameter_length, likelihood_compose, likelihood_results_length);

	printf("starting evaluator...\n");
	mpi_evaluator__start();

	get_min_parameters(ap, &min_parameters);
	get_max_parameters(ap, &max_parameters);
	get_search_parameters(ap, &point);
	get_step(ap, &step);

	printf("searching...\n");
	register_search(asynchronous_newton_method);
//	register_search("gs", init_genetic_search);
//	register_search("de", init_differential_evolution);
//	register_search("pso", init_particle_swarm);

	for (i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-s")) {
			char search_name[SEARCH_NAME_SIZE], *search_qualifier;
			strcpy(search_name, argv[++i]);
			get_qualifier_from_name(search_name, &search_qualifier);

			if (!strcmp(search_qualifier, "nm")) {
				printf("creating newton method...\n");
				create_newton_method(search_name, 10, 300, ap->number_parameters, point, step);
				printf("created.\n");
			} else if (!strcmp(search_qualifier, "gs")) {
//				printf("creating genetic search...\n");
//				create_genetic_search(search_name, ap->number_parameters, min_parameters, max_parameters, 20000, 500);
//				pritnf("created.\n");
			} else if (!strcmp(search_qualifier, "de")) {
			} else if (!strcmp(search_qualifier, "pso")) {
			}
			free(search_qualifier);
		}
	}

	start_mpi_search_manager(argc, argv);

	return 0;
}
