/********
	*	Includes for astronomy
 ********/
#include "mpi.h"

#include "parameters.h"
#include "star_points.h"
#include "evaluation.h"

#include "../searches/gradient_descent.h"
#include "../searches/synchronous_search.h"
#include "../searches/genetic_search.h"
#include "../searches/differential_evolution.h"
#include "../searches/particle_swarm.h"
#include "../searches/newton_method.h"
#include "../evaluation/mpi_evaluator.h"
#include "../evaluation/evaluator.h"

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
	initialize_state(es, ap->number_streams);

}

void integral_f(double* parameters, double** results) {
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
	printf("[worker: %d] composed likelihood: %lf\n", get_mpi_rank(), prob_sum);
	return prob_sum;
}

int main(int number_arguments, char **arguments){
	int integral_parameter_length, integral_results_length;
	int likelihood_parameter_length, likelihood_results_length;
	double *min_parameters;
	double *max_parameters;
	double *point;
	double *step;

	evaluator__init(&number_arguments, &arguments, read_data);

	integral_parameter_length = ap->number_parameters;
	integral_results_length = 1 + ap->number_streams;
	likelihood_parameter_length = 1 + ap->number_streams;
	likelihood_results_length = 2;

	evaluator__init_integral(integral_f, integral_parameter_length, integral_compose, integral_results_length);
	evaluator__init_likelihood(likelihood_f, likelihood_parameter_length, likelihood_compose, likelihood_results_length);
	printf("starting...\n");
	mpi_evaluator__start();

	get_min_parameters(ap, &min_parameters);
	get_max_parameters(ap, &max_parameters);
	get_search_parameters(ap, &point);
	get_step(ap, &step);

	printf("searching...\n");
	if (arguments[2][0] == 'g' && arguments[2][1] == 'd') {
		synchronous_gradient_descent(arguments[1], arguments[2], point, step, ap->number_parameters);
	} else if (arguments[2][0] == 'c') {
		synchronous_conjugate_gradient_descent(arguments[1], arguments[2], point, step, ap->number_parameters);
	} else if (arguments[2][0] == 'n') {
		synchronous_newton_method(arguments[1], arguments[2], point, step, ap->number_parameters);
	} else if (arguments[2][0] == 'g' && arguments[2][1] == 's') {
		synchronous_search(arguments[1], arguments[2], min_parameters, max_parameters, ap->number_parameters, start_genetic_search);
	} else if (arguments[2][0] == 'd') {
		synchronous_search(arguments[1], arguments[2], min_parameters, max_parameters, ap->number_parameters, start_differential_evolution);
	} else if (arguments[2][0] == 'p') {
		synchronous_search(arguments[1], arguments[2], min_parameters, max_parameters, ap->number_parameters, start_particle_swarm);
	}
	return 0;
}
