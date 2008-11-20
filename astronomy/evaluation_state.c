/*
 *  integral_function.c
 *  Astronomy
 *
 *  Created by Travis Desell on 2/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

/****
        *       BOINC includes
*****/

#ifdef GMLE_BOINC
	#ifdef _WIN32
		#include "boinc_win.h"
	#else
		#include "config.h"
	#endif

	#ifndef _WIN32
		#include <cstdio>
		#include <cctype>
		#include <cstring>
		#include <cstdlib>
		#include <csignal>
	#endif

	#ifdef BOINC_APP_GRAPHICS
		#include "graphics_api.h"
		#include "graphics_lib.h"
	#endif

	#include "diagnostics.h"
	#include "util.h"
	#include "filesys.h"
	#include "boinc_api.h"
	#include "mfile.h"
#endif

#define CHECKPOINT_FILE "astronomy_checkpoint"

/****
	*	Astronomy includes
*****/
#include <math.h>
#include <time.h>
#include <stdio.h>

#include "evaluation_optimized.h"
#include "parameters.h"
#include "probability.h"
#include "stCoords.h"
#include "atSurveyGeometry.h"
#include "star_points.h"
#include "numericalIntegration.h"
#include "../evaluation/mpi_evaluator.h"

#ifndef _WIN32
	#define pi M_PI
#else
	#define pi 3.14159265358979323846
#endif

#define deg (180.0/pi)

void print_integral_area(FILE *file, INTEGRAL_AREA *ia) {
	int i;
	fprintf(file, "mu[min,steps,step_size]: %.10lf, %d, %.10lf\n", ia->mu_min, ia->mu_steps, ia->mu_step_size);
	fprintf(file, "nu[min,steps,step_size]: %.10lf, %d, %.10lf\n", ia->nu_min, ia->nu_steps, ia->nu_step_size);
	fprintf(file, " r[min,steps,step_size]: %.10lf, %d, %.10lf\n", ia->r_min, ia->r_steps, ia->r_step_size);
	fprintf(file, "background integral: %.10lf\n", ia->background_integral);
	for(i = 0; i < sizeof(ia->stream_integrals)/sizeof(double); i++) {
		fprintf(file, "stream integral[%d]: %.10lf\n", i, ia->stream_integrals[i]);
	}
}

void initialize_integal_area(INTEGRAL_AREA *ia, double mu_min, int mu_steps, double mu_step_size, double nu_min, int nu_steps, double nu_step_size, double r_min, int r_steps, double r_step_size, int number_streams) {
	int i;

	ia->mu_min		= mu_min;
	ia->mu_steps		= mu_steps;
	ia->mu_step_size	= mu_step_size;
	ia->mu_step_current	= 0;
	ia->nu_min		= nu_min;
	ia->nu_steps		= nu_steps;
	ia->nu_step_size	= nu_step_size;
	ia->nu_step_current	= 0;
	ia->r_min		= r_min;
	ia->r_steps		= r_steps;
	ia->r_step_size		= r_step_size;
	ia->r_step_current	= 0;

	ia->background_integral	= 0;
	ia->stream_integrals	= (double*)malloc(sizeof(double) * number_streams);
	for (i = 0; i < number_streams; i++) {
		ia->stream_integrals[i] = 0;
	}
}

void initialize_state(ASTRONOMY_PARAMETERS *ap, EVALUATION_STATE *es) {
	int i;

	es->current_cut = -1;
	es->background_integral = 0;
	es->stream_integrals = (double*)malloc(sizeof(double) * ap->number_streams);
	for (i = 0; i < ap->number_streams; i++) es->stream_integrals[i] = 0;

	es->current_star_point = 0;
	es->num_zero = 0;
	es->bad_jacobians = 0;
	es->prob_sum = 0;

	es->main_integral = (INTEGRAL_AREA*)malloc(sizeof(INTEGRAL_AREA));
	initialize_integal_area(es->main_integral, ap->mu_min, ap->mu_steps, ap->mu_step_size, ap->nu_min, ap->nu_steps, ap->nu_step_size, ap->r_min, ap->r_steps, ap->r_step_size, ap->number_streams);

	es->cuts = (INTEGRAL_AREA**)malloc(sizeof(INTEGRAL_AREA*) * ap->number_cuts);
	for (i = 0; i < ap->number_cuts; i++) {
		es->cuts[i] = (INTEGRAL_AREA*)malloc(sizeof(INTEGRAL_AREA));
		initialize_integal_area(es->cuts[i], ap->mu_cut[i][0], ap->mu_cut[i][2], ap->mu_cut_step_size[i], ap->nu_cut[i][0], ap->nu_cut[i][2], ap->nu_cut_step_size[i], ap->r_cut[i][0], ap->r_cut[i][2], ap->r_cut_step_size[i], ap->number_streams);
	}
}

void reset_evaluation_state(ASTRONOMY_PARAMETERS *ap, EVALUATION_STATE *es) {
	int i, j;

	es->current_cut = -1;
	es->background_integral = 0;
	es->current_star_point = 0;
	es->num_zero = 0;
	es->bad_jacobians = 0;
	es->prob_sum = 0;

	es->main_integral->mu_step_current = 0;
	es->main_integral->nu_step_current = 0;
	es->main_integral->r_step_current = 0;
	es->main_integral->background_integral = 0;
	for (i = 0; i < ap->number_streams; i++) {
		es->stream_integrals[i] = 0;
		es->main_integral->stream_integrals[i] = 0;
	}

	for (i = 0; i < ap->number_cuts; i++) {
		es->cuts[i]->mu_step_current = 0;
		es->cuts[i]->nu_step_current = 0;
		es->cuts[i]->r_step_current = 0;
		es->cuts[i]->background_integral = 0;
		for (j = 0; j < ap->number_streams; j++) {
			es->cuts[i]->stream_integrals[i] = 0;
		}
	}
}

void free_integral_area(INTEGRAL_AREA *ia) {
	free(ia->stream_integrals);
}

void free_state(ASTRONOMY_PARAMETERS *ap, EVALUATION_STATE* es) {
	int i;
	free(es->stream_integrals);
	free_integral_area(es->main_integral);
	for (i = 0; i < ap->number_cuts; i++) {
		free_integral_area(es->cuts[i]);
	}
	free(es->main_integral);
	free(es->cuts);
}

#ifdef GMLE_BOINC
	int write_checkpoint(EVALUATION_STATE* es) {
		int i;
		char output_path[512];
		boinc_resolve_filename(CHECKPOINT_FILE, output_path, sizeof(output_path));

		MFILE cp;
		int retval = cp.open(output_path, "w");
		if (retval) {
	                fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file) %d\n", retval);
	                return retval;
		}

		cp.printf("r_step_current: %d\n", es->r_step_current);
		cp.printf("mu_step_current: %d\n", es->mu_step_current);
		cp.printf("nu_step_current: %d\n", es->nu_step_current);

		cp.printf("main_integral_calculated: %d\n", es->main_integral_calculated);
		cp.printf("current_cut: %d\n", es->current_cut);
		cp.printf("background_integral: %lf\n", es->background_integral);
		cp.printf("stream_integrals [%d]:", es->number_streams);
		for (i = 0; i < es->number_streams; i++) cp.printf(" %lf", es->stream_integrals[i]);
		cp.printf("\ncurrent_star_point: %d\n", es->current_star_point);
		cp.printf("num_zero: %d\n", es->num_zero);
		cp.printf("bad_jacobians: %d\n", es->bad_jacobians);
		cp.printf("prob_sum: %lf\n", es->prob_sum);

		retval = cp.flush();
		if (retval) {
	                fprintf(stderr, "APP: error writing checkpoint (flushing checkpoint file) %d\n", retval);
	                return retval;
		}
		retval = cp.close();
		if (retval) {
	                fprintf(stderr, "APP: error writing checkpoint (closing checkpoint file) %d\n", retval);
	                return retval;
		}
		return 0;
	}

	int read_checkpoint(EVALUATION_STATE* es) {
		int i;
		char input_path[512];
		int retval = boinc_resolve_filename(CHECKPOINT_FILE, input_path, sizeof(input_path));
		if (retval) {
			return 0;
		}

	        FILE* data_file = boinc_fopen(input_path, "r");
	        if (!data_file) {
	                fprintf(stderr, "APP: error reading checkpoint (opening file)\n");
	                return 1;
	        }
		fscanf(data_file, "r_step_current: %d\n", &es->r_step_current);
		fscanf(data_file, "mu_step_current: %d\n", &es->mu_step_current);
		fscanf(data_file, "nu_step_current: %d\n", &es->nu_step_current);
	
		fscanf(data_file, "main_integral_calculated: %d\n", &es->main_integral_calculated);
		fscanf(data_file, "current_cut: %d\n", &es->current_cut);
		fscanf(data_file, "background_integral: %lf\n", &es->background_integral);
		fscanf(data_file, "stream_integrals [%d]:", &es->number_streams);
		for (i = 0; i < es->number_streams; i++) fscanf(data_file, " %lf", &es->stream_integrals[i]);

		fscanf(data_file, "\ncurrent_star_point: %d\n", &es->current_star_point);
		fscanf(data_file, "num_zero: %d\n", &es->num_zero);
		fscanf(data_file, "bad_jacobians: %d\n", &es->bad_jacobians);
		fscanf(data_file, "prob_sum: %lf\n", &es->prob_sum);

	        fclose(data_file);
	        return 0;
	}
#endif
