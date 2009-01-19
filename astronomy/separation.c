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
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "parameters.h"
#include "star_points.h"
#include "evaluation.h"
#include "stVector.h"
#include "probability.h"
#include "stCoords.h"

#include "../evaluation/mpi_evaluator.h"
#include "../evaluation/evaluator.h"

#define max_iterations			35000
#define astronomy_parameters_file	"parameters-20.txt"
#define star_points_file		"20-f-cut-noMon.txt"
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
	(*results) = (double*)malloc(sizeof(double) * 1+ap->number_streams);
	(*results)[0] = es->background_integral;
	for(i = 0; i < ap->number_streams; i++) {
		(*results)[i+1] = es->stream_integrals[i];
	}
//	printf("calculated integrals: %lf, %lf\n", (*results)[0], (*results)[1]);
}

double integral_compose(double* integral_results, int num_results) {
	int i, j;

	es->background_integral = 0.0;
	for (i = 0; i < ap->number_streams; i++) {
		es->stream_integrals[i] = 0.0;
	}
	for (i = 0; i < num_results; i++) {
		es->background_integral += integral_results[((ap->number_streams+1)*i)];
		for (j = 0; j < ap->number_streams; j++) {
			es->stream_integrals[j] += integral_results[((ap->number_streams+1)*i)+j+1];
		}
	}
        printf("background integral: %lf, stream integrals:", es->background_integral);
        for (i = 0; i < ap->number_streams; i++) printf(" %lf", es->stream_integrals[i]);
        printf("\n");
	return -1;
}


void separation(char* filename, double background_integral, double* stream_integrals) {
	int q[ap->number_streams];
	double nstars[ap->number_streams];
	int total;
	double sprob[ap->number_streams];
	double prob_s[ap->number_streams];
	double prob_b;
	double pbx;
	double psg[ap->number_streams];
	double d;
	int twoPanel;
	double **cmatrix;
	double dnormal[3];
	double dortho[3];
	double xsun[3];
	double epsilon_s[ap->number_streams];
	double epsilon_b;
	double *star_coords;
	double starxyz[3];
	double starxyzTransform[3];
	int s_ok;
	int i, j, retval;
	FILE *file;

	twoPanel = 1;
	for(j = 0; j < ap->number_streams; j++) {
		nstars[j] = 0;
		q[j] = 0;
	}
	total = 0;
	prob_ok_init();

	printf("Integral complete.\n Beginning probability calculations...\n");
	file = fopen(filename, "w");

        if(ap->sgr_coordinates == 0){
               stripe_normal(ap->wedge, dnormal);
        }else if (ap->sgr_coordinates == 1) {
               sgr_stripe_normal(ap->wedge, dnormal);
        } else {
               printf("Error: ap->sgr_coordinates not valid");
        }

	free_star_points(sp);
	free(sp);
        sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
        retval = read_star_points(star_points_file, sp);
        if (retval) {
                fprintf(stderr, "APP: error reading star points: %d\n", retval);
                exit(1);
        }
        printf("read %d stars.\n", sp->number_stars);


	cmatrix = (double**)malloc(sizeof(double*) * 3);
	for (i = 0; i < 3; i++) cmatrix[i] = (double*)malloc(sizeof(double) * 3);
	dortho[0] = 0.0;
	dortho[1] = 0.0;
	dortho[2] = 1.0;
	get_transform(dnormal, dortho, cmatrix);
	
	printf("\nTransformation matrix:\n");
	printf("\t%lf %lf %lf\n", cmatrix[0][0], cmatrix[0][1], cmatrix[0][2]);
	printf("\t%lf %lf %lf\n", cmatrix[1][0], cmatrix[1][1], cmatrix[1][2]);
	printf("\t%lf %lf %lf\n\n", cmatrix[2][0], cmatrix[2][1], cmatrix[2][2]);

	xsun[0] = -8.5;
	xsun[1] = 0.0;
	xsun[2] = 0.0;
	d = dotp(dnormal, xsun);

	printf("==============================================\n");
	printf("bint: %lf", background_integral);
	for(j = 0; j < ap->number_streams; j++) { 
		printf(", ");
		printf("sint[%d]: %lf", j, stream_integrals[j]);
	} 
	printf("\n");

	/*get stream & background weight constants*/
	double denom = 1.0;
	for(j = 0; j < ap->number_streams; j++) {
		denom += exp(ap->stream_weights[j]);
	}
	for(j = 0; j < ap->number_streams; j++) {
		epsilon_s[j] = exp(ap->stream_weights[j]) / denom;
		printf("epsilon_s[%d]: %lf\n", j, epsilon_s[j]);
	}
	epsilon_b = 1.0 / denom;
	printf("epsilon_b:    %lf\n", epsilon_b);
	
	for (i = 0; i < sp->number_stars; i++) {
		star_coords = sp->stars[i];
		//printf("star_coords: %g %g %g\n", star_coords[0], star_coords[1], star_coords[2]);
	
		//printf("twoPanel: %d\n", twoPanel);

		if (twoPanel == 1) {
			if (ap->convolve != 0) {
				for(j = 0; j < ap->number_streams; j++) {
					prob_s[j] = stPsgConvolved(star_coords, ap->stream_parameters[j], ap->wedge, ap->convolve, ap->sgr_coordinates);
				}
				prob_b = stPbxConvolved(star_coords, ap->background_parameters, ap->convolve, ap->wedge);
			} else {
				for(j = 0; j < ap->number_streams; j++) {
                        		prob_s[j] = stPsg(star_coords, ap->stream_parameters[j], ap->wedge, ap->sgr_coordinates);
				}
				prob_b = stPbx(star_coords, ap->background_parameters);
                   	}
		
			//printf("prob_s: %g\n", prob_s);
			//printf("prob_b: %g\n", prob_b);
	
			pbx = epsilon_b * prob_b / background_integral;
			
			for(j = 0; j < ap->number_streams; j++) {
				psg[j] = epsilon_s[j] * prob_s[j] / stream_integrals[j]; 
			}

			//printf("pbx: %g\n", pbx);
			//printf("psg: %g\n", psg);
				
			double psgSum = 0;
			for(j = 0; j < ap->number_streams; j++) {
				psgSum += psg[j];
			}
			for(j = 0; j < ap->number_streams; j++) {
				sprob[j] = psg[j] / (psgSum + pbx);
			}

			//printf("sprob: %g\n", sprob);

			for(j = 0; j < ap->number_streams; j++) {
				nstars[j] += sprob[j];
			}
		} else {	
			for(j = 0; j < ap->number_streams; j++) {
				sprob[j] = 1.0;
				nstars[j] += 1.0;
			}
		}


		/*determine if star with sprob should be put into stream*/
		for(j = 0; j < ap->number_streams; j++) {
			s_ok = prob_ok(sprob[j]);
			if (s_ok == 1) {	
				s_ok += j;
				break;
			}
		}

		//printf("s_ok: %d\n", s_ok);

		if (s_ok >= 1) {
			q[s_ok-1]++;
		}

		lbr2xyz(star_coords, starxyz);
		transform_point(starxyz, cmatrix, xsun, starxyzTransform);	

		fprintf(file, "%d %lf %lf %lf\n", s_ok, starxyzTransform[0], starxyzTransform[1], starxyzTransform[2]);
		//free(starxyz);
		//free(starxyzTransform);

		total += 1;


		if( (total % 10000) == 0 ) printf("%d\n", total);
	}
		
	printf("%d total stars\n", total); 
	for(j=0; j<ap->number_streams;j++) {
		printf("%lf in stream[%d] (%lf%%)\n", nstars[j], j, (nstars[j]/total*100));
	}

	for(j=0; j<ap->number_streams;j++) {
		printf("%d stars separated into stream\n", q[j]);
	}
	fclose(file);
	printf("Output written to: %s\n", filename);
}

int main(int number_arguments, char **arguments){
	int integral_parameter_length, integral_results_length;
	double *point;

	printf("init data...\n");
	mpi_evaluator__init(&number_arguments, &arguments);
	mpi_evaluator__read_data(read_data);

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
