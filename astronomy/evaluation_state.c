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
#include "../util/io_util.h"

#ifndef _WIN32
	#define pi M_PI
#else
	#define pi 3.14159265358979323846
#endif

#define deg (180.0/pi)


void fwrite_integral_area(FILE *file, INTEGRAL_AREA *ia) {
	fprintf(file, "mu[min,max,steps]: %.3lf, %.3lf, %d\n", ia->mu_min, ia->mu_max, ia->mu_steps);
	fprintf(file, "nu[min,max,steps]: %.3lf, %.3lf, %d\n", ia->nu_min, ia->nu_max, ia->nu_steps);
	fprintf(file, " r[min,max,steps]: %.3lf, %.3lf, %d\n", ia->r_min, ia->r_max, ia->r_steps);
	fprintf(file, "min_calculation: %ld, max_calculation: %ld, current_calculation: %ld\n", ia->min_calculation, ia->max_calculation, ia->current_calculation);
	fprintf(file, "background_integral: %.20lf\n", ia->background_integral);
	print_double_array(file, "stream_integrals", ia->number_streams, ia->stream_integrals);
}

void fread_integral_area(FILE *file, INTEGRAL_AREA *ia) {
	int i;

	fscanf(file, "mu[min,max,steps]: %lf, %lf, %d\n", &(ia->mu_min), &(ia->mu_max), &(ia->mu_steps));
	fscanf(file, "nu[min,max,steps]: %lf, %lf, %d\n", &(ia->nu_min), &(ia->nu_max), &(ia->nu_steps));
	fscanf(file, " r[min,max,steps]: %lf, %lf, %d\n", &(ia->r_min), &(ia->r_max), &(ia->r_steps));
	ia->mu_step_size = (ia->mu_max - ia->mu_min) / ia->mu_steps;
	ia->nu_step_size = (ia->nu_max - ia->nu_min) / ia->nu_steps;
	ia->r_step_size = (ia->r_max - ia->r_min) / ia->r_steps;
	fscanf(file, "min_calculation: %ld, max_calculation: %ld, current_calculation: %ld\n", &(ia->min_calculation), &(ia->max_calculation), &(ia->current_calculation));
	fscanf(file, "background_integral: %lf\n", &(ia->background_integral));
	fscanf(file, "stream_integrals[%d]: ", &(ia->number_streams));
	for (i = 0; i < ia->number_streams; i++) {
		fscanf(file, "%lf", &(ia->stream_integrals[i]));
		if (i != ia->number_streams-1) fscanf(file, ", ");
	}
}

void initialize_integral_area(INTEGRAL_AREA *ia, INTEGRAL *integral, int number_streams) {
	int i;

	ia->mu_min		= integral->mu_min;
	ia->mu_max		= integral->mu_max;
	ia->mu_steps		= integral->mu_steps;
	ia->nu_min		= integral->nu_min;
	ia->nu_max		= integral->nu_max;
	ia->nu_steps		= integral->nu_steps;
	ia->r_min		= integral->r_min;
	ia->r_max		= integral->r_max;
	ia->r_steps		= integral->r_steps;

	ia->min_calculation	= integral->min_calculation;
	ia->max_calculation	= integral->max_calculation;
	ia->current_calculation	= integral->min_calculation;

	ia->mu_step_size = (ia->mu_max - ia->mu_min) / ia->mu_steps;
	ia->nu_step_size = (ia->nu_max - ia->nu_min) / ia->nu_steps;
	ia->r_step_size = (ia->r_max - ia->r_min) / ia->r_steps;

	ia->number_streams = number_streams;
	ia->background_integral	= 0;
	ia->stream_integrals	= (double*)malloc(sizeof(double) * number_streams);
	for (i = 0; i < number_streams; i++) {
		ia->stream_integrals[i] = 0;
	}
}

void get_steps(INTEGRAL_AREA *ia, int *mu_step, int *nu_step, int *r_step) {
	(*r_step) = ia->current_calculation % ia->r_steps;
//	printf("r_step = [current: %ld] mod [r_steps: %d] = %d\n", ia->current_calculation, ia->r_steps, (*r_step));

	(*nu_step) = (ia->current_calculation / ia->r_steps) % ia->nu_steps;
//	printf("nu_step = [current: %ld] / [r_steps: %d] mod [nu_steps: %d] = %d\n", ia->current_calculation, ia->r_steps, ia->nu_steps, (*nu_step));

	(*mu_step) = ia->current_calculation / (ia->r_steps * ia->nu_steps);
//	printf("mu_step = [current: %ld] / [r_steps: %d * nu_steps: %d] = %d\n", ia->current_calculation, ia->r_steps, ia->nu_steps, (*mu_step));
}

void initialize_state(ASTRONOMY_PARAMETERS *ap, STAR_POINTS *sp, EVALUATION_STATE *es) {
	int i;

	es->current_integral = 0;
	es->background_integral = 0;
	es->stream_integrals = (double*)malloc(sizeof(double) * ap->number_streams);
	for (i = 0; i < ap->number_streams; i++) es->stream_integrals[i] = 0;

	es->number_streams = ap->number_streams;
	es->total_stars = sp->number_stars;
	es->current_star_point = 0;
	es->num_zero = 0;
	es->bad_jacobians = 0;
	es->prob_sum = 0;

	es->number_integrals = ap->number_integrals;
	es->integral = (INTEGRAL_AREA**)malloc(sizeof(INTEGRAL_AREA*) * ap->number_integrals);
	for (i = 0; i < ap->number_integrals; i++) {
		es->integral[i] = (INTEGRAL_AREA*)malloc(sizeof(INTEGRAL_AREA));
		initialize_integral_area(es->integral[i], ap->integral[i], ap->number_streams);
	}
}

void reset_evaluation_state(EVALUATION_STATE *es) {
	int i, j;

	es->current_integral = 0;
	es->background_integral = 0;
	for (i = 0; i < es->number_streams; i++) es->stream_integrals[i] = 0;
	es->current_star_point = 0;
	es->num_zero = 0;
	es->bad_jacobians = 0;
	es->prob_sum = 0;

	for (i = 0; i < es->number_integrals; i++) {
		es->integral[i]->background_integral = 0;
		for (j = 0; j < es->integral[i]->number_streams; j++) {
			es->integral[i]->stream_integrals[j] = 0;
		}
		es->integral[i]->current_calculation = es->integral[i]->min_calculation;
	}
}

void free_integral_area(INTEGRAL_AREA *ia) {
	free(ia->stream_integrals);
}

void free_state(EVALUATION_STATE* es) {
	int i;
	free(es->stream_integrals);
	for (i = 0; i < es->number_integrals; i++) {
		free_integral_area(es->integral[i]);
	}
	free(es->integral);
}

#ifdef GMLE_BOINC
	int write_checkpoint(EVALUATION_STATE* es) {
		int i, retval;
		char output_path[512];
		FILE *file;

		boinc_resolve_filename(CHECKPOINT_FILE, output_path, sizeof(output_path));

		file = boinc_fopen(output_path, "w+");
		if (!file) {
			fprintf(stderr, "APP: error writing checkpoint (opening checkpoint file)\n");
			return 1;
		}

		fprintf(file, "background_integral: %.20lf\n", es->background_integral);
		fprintf(file, "stream_integrals[%d]: ", es->number_streams);
		for (i = 0; i < es->number_streams; i++) {
			fprintf(file, "%.20lf", es->stream_integrals[i]);
			if (i != (es->number_streams-1)) fprintf(file, ", ");
		}
		fprintf(file, "\n");

		fprintf(file, "prob_sum: %.20lf, num_zero: %d, bad_jacobians: %d\n", es->prob_sum, es->num_zero, es->bad_jacobians);
		fprintf(file, "current_star_point: %d\n", es->current_star_point);
		fprintf(file, "current_integral: %d\n", es->current_integral);
		fprintf(file, "number_integrals: %d\n", es->number_integrals);
		for (i = 0; i < es->number_integrals; i++) {
			fwrite_integral_area(file, es->integral[i]);
		}

		if ((retval = fclose(file))) {
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

	        FILE* file = boinc_fopen(input_path, "r");
	        if (file == NULL) {
	                return 0;
	        }

		if (1 > fscanf(file, "background_integral: %lf\n", &(es->background_integral))) return 1;
		es->number_streams = read_double_array(file, "stream_integrals", &(es->stream_integrals));

		fscanf(file, "prob_sum: %lf, num_zero: %d, bad_jacobians: %d\n", &(es->prob_sum), &(es->num_zero), &(es->bad_jacobians));
		fscanf(file, "current_star_point: %d\n", &(es->current_star_point));
		fscanf(file, "current_integral: %d\n", &(es->current_integral));
		fscanf(file, "number_integrals: %d\n", &(es->number_integrals));
		for (i = 0; i < es->number_integrals; i++) {
			fread_integral_area(file, es->integral[i]);
		}

	        fclose(file);
	        return 0;
	}
#endif
