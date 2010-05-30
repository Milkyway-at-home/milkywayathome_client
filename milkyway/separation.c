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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "astronomy_worker.h"
#include "evaluation_optimized.h"
#include "parameters.h"
#include "probability.h"
#include "star_points.h"
#include "stCoords.h"
#include "stVector.h"

#include "../evaluation/evaluator.h"
#include "../evaluation/mpi_evaluator.h"


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
	double star_coords[3];
	double starxyz[3];
	double starxyzTransform[3];
	int s_ok = 0;
	int i, j, retval;
	FILE *file;
	double reff_xr_rp3, *qw_r3_N, *r_point;

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

	r_point = (double*)malloc(sizeof(double) * ap->convolve);
	qw_r3_N = (double*)malloc(sizeof(double) * ap->convolve);

	init_constants(ap);

	printf("initialized constants\n");

	for (i = 0; i < sp->number_stars; i++) {
//		printf("[%d/%d] setting star coords\n", i, sp->number_stars);
		star_coords[0] = sp->stars[i][0];
		star_coords[1] = sp->stars[i][1];
		star_coords[2] = sp->stars[i][2];
//		printf("star_coords: %g %g %g\n", star_coords[0], star_coords[1], star_coords[2]);

//		printf("twoPanel: %d\n", twoPanel);

		if (twoPanel == 1) {
//			printf("setting probability constants\n");
			set_probability_constants(ap->convolve, star_coords[2], r_point, qw_r3_N, &reff_xr_rp3);
//			printf("calculating probabilities\n");
			calculate_probabilities(r_point, qw_r3_N, reff_xr_rp3, star_coords, ap, &prob_b, prob_s);
//			printf("calculated probabilities\n");

//			printf("prob_s: %lf\n", prob_s[0]);
//			printf("prob_b: %lf\n", prob_b);

			pbx = epsilon_b * prob_b / background_integral;

			for(j = 0; j < ap->number_streams; j++) {
				psg[j] = epsilon_s[j] * prob_s[j] / stream_integrals[j];
			}

//			printf("pbx: %g\n", pbx);
//			printf("psg: %g\n", psg[0]);

			double psgSum = 0;
			for(j = 0; j < ap->number_streams; j++) {
				psgSum += psg[j];
			}
			for(j = 0; j < ap->number_streams; j++) {
				sprob[j] = psg[j] / (psgSum + pbx);
			}

//			printf("sprob: %g\n", sprob[0]);

			for(j = 0; j < ap->number_streams; j++) {
				nstars[j] += sprob[j];
			}
//			printf("nstars: %g\n", nstars[0]);
		} else {
			for(j = 0; j < ap->number_streams; j++) {
				sprob[j] = 1.0;
				nstars[j] += 1.0;
			}
		}


		/*determine if star with sprob should be put into stream*/
		//for(j = 0; j < ap->number_streams; j++) {
			s_ok = prob_ok(ap->number_streams, sprob);
		//	if (s_ok == 1) {
		//		s_ok += j;
		//		break;
		//	}
		//}

//		printf("s_ok: %d\n", s_ok);

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
	free(r_point);
	free(qw_r3_N);
	free_constants(ap);
}

int main(int number_arguments, char **arguments){
	int integral_parameter_length, integral_results_length, likelihood_parameter_length, likelihood_results_length, i, number_parameters;
	double *point;
	char outfile[512];

	printf("init data...\n");
	mpi_evaluator__init(&number_arguments, &arguments);

	sprintf(star_points_file, "stars.txt");
	sprintf(astronomy_parameters_file, "parameters.txt");
	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-stars")) {
			sprintf(star_points_file, "%s", arguments[++i]);
		} else if (!strcmp(arguments[i], "-parameters")) {
			sprintf(astronomy_parameters_file, "%s", arguments[++i]);
		} else if (!strcmp(arguments[i], "-outfile")) {
			sprintf(outfile, "%s", arguments[++i]);
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

	printf("init integral...\n");
	evaluator__init_integral(integral_f, integral_parameter_length, integral_compose, integral_results_length);
	printf("init likelihood...\n");
	evaluator__init_likelihood(likelihood_f, likelihood_parameter_length, likelihood_compose, likelihood_results_length);

	printf("starting...\n");
	mpi_evaluator__start();

	printf("getting parameters...\n");
	get_search_parameters(ap, &point);

	printf("evaluating point...\n");
	evaluate(point);

	printf("performing separation...\n");
	separation(outfile, es->background_integral, es->stream_integrals);

	return 0;
}

