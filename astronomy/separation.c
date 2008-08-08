/********
	*	Includes for astronomy
 ********/
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "parameters.h"
#include "star_points.h"
#include "evaluation.h"
#include "stVector.h"
#include "probability.h"
#include "stCoords.h"

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
	(*results) = (double*)malloc(sizeof(double) * 2);
	(*results)[0] = es->background_integral;
	(*results)[1] = es->stream_integrals[0];
//	printf("calculated integrals: %lf, %lf\n", (*results)[0], (*results)[1]);
}

double integral_compose(double* integral_results, int num_results) {
	int i;
	es->background_integral = 0.0;
	es->stream_integrals[0] = 0.0;
	for (i = 0; i < num_results; i++) {
		es->background_integral += integral_results[(2*i)];
		es->stream_integrals[0] += integral_results[(2*i)+1];
	}
	printf("composed integrals: %lf, %lf\n", es->background_integral, es->stream_integrals[0]);
	return -1;
}


void separation(char* filename, double background_integral, double* stream_integrals) {
	int q;
	double nstars;
	int total;
	double sprob;
	double prob_s;
	double prob_b;
	double pbx;
	double psg;
	double d;
	int twoPanel;
	double **cmatrix;
	double dnormal[3];
	double dortho[3];
	double xsun[3];
	double epsilon_s, epsilon_b;
	double *star_coords;
	double starxyz[3];
	double starxyzTransform[3];
	int s_ok;
	int i;
	FILE *file;

	twoPanel = 1;
	nstars = 0;
	total = 0;
	q = 0;

	printf("Integral complete.\n Beginning probability calculations...\n");
	file = fopen(filename, "w");

	stripe_normal(ap->wedge, dnormal);

	cmatrix = (double**)malloc(sizeof(double*) * 3);
	for (i = 0; i < 3; i++) cmatrix[i] = (double*)malloc(sizeof(double) * 3);
	dortho[0] = 0.0;
	dortho[1] = 0.0;
	dortho[2] = 1.0;
	get_transform(dnormal, dortho, cmatrix);
	
	printf("\nTransformation matrix:\n");
	printf("\t%lf %lf %lf\n", cmatrix[0][0], cmatrix[0][1], cmatrix[0][2]);
	printf("\t%lf %lf %lf\n", cmatrix[1][0], cmatrix[1][1], cmatrix[1][2]);
	printf("\t%lf %lf %lf\n", cmatrix[2][0], cmatrix[2][1], cmatrix[2][2]);

	xsun[0] = -8.5;
	xsun[1] = 0.0;
	xsun[2] = 0.0;
	d = dotp(dnormal, xsun);

	printf("bint: %lf, sint: %lf\n", background_integral, stream_integrals[0]);

	/*get stream & background weight constants*/
	epsilon_s = exp(ap->stream_weights[0]) / (1.0 + exp(ap->stream_weights[0]));
	epsilon_b = 1.0 / (1 + exp(ap->stream_weights[0]));
	printf("epsilon_s: %lf\n", epsilon_s);
	printf("epsilon_b: %lf\n", epsilon_b);
	
	for (i = 0; i < sp->number_stars; i++) {
		star_coords = sp->stars[i];

		if (twoPanel == 1) {
			if (ap->convolve != 0) {
				prob_s = stPsgConvolved(star_coords, ap->stream_parameters[0], ap->wedge, ap->convolve, ap->sgr_coordinates);
				prob_b = stPbxConvolved(star_coords, ap->background_parameters, ap->convolve, ap->wedge);
			} else {
                        	prob_s = stPsg(star_coords, ap->stream_parameters[0], ap->wedge, ap->sgr_coordinates);
				prob_b = stPbx(star_coords, ap->background_parameters);
                   	}

			pbx = epsilon_b * prob_b / background_integral;
			psg = epsilon_s * prob_s / stream_integrals[0]; 

			sprob = psg / (psg + pbx);
			nstars += sprob;
		} else {	
			sprob = 1.0;
			nstars += sprob;
		}


		/*determine if star with sprob should be put into stream*/	
		s_ok = prob_ok(sprob);
		if (s_ok == 1) q++;

		lbr2xyz(star_coords, starxyz);
		transform_point(starxyz, cmatrix, xsun, starxyzTransform);	

		fprintf(file, "%d %lf %lf %lf\n", s_ok, starxyzTransform[0], starxyzTransform[1], starxyzTransform[2]);
		//free(starxyz);
		//free(starxyzTransform);

		total += 1;


		if( (total % 10000) == 0 ) printf("%d\n", total);
	}
		
	printf("%d total stars, %lf in stream (%lf %% )\n", total, nstars, (nstars/total*100));
	printf("%d stars separated into stream\n", q);
	fclose(file);
	printf("Output written to: %s\n", filename);
}

int main(int number_arguments, char **arguments){
	int integral_parameter_length, integral_results_length;
	double *point;

	printf("init data...\n");
	evaluator__init(&number_arguments, &arguments, read_data);

	integral_parameter_length = ap->number_parameters;
	integral_results_length = 1 + ap->number_streams;

	printf("init integral...\n");
	evaluator__init_likelihood(integral_f, integral_parameter_length, integral_compose, integral_results_length);
	printf("starting...\n");
	mpi_evaluator__start();

	printf("getting parameters...\n");
	get_search_parameters(ap, &point);

	evaluate(point);
	separation(arguments[1], es->background_integral, es->stream_integrals);

	return 0;
}
