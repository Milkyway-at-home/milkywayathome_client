/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gradient.h"
#include "../evaluation/evaluator.h"
#include "../util/settings.h"
#include "../util/io_util.h"

#ifdef BOINC_APPLICATION
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

int checkpoint_gradient(char *checkpoint_file, int number_parameters, int j, double *gradient) {
	FILE *file;
	#ifdef BOINC_APPLICATION 
		char input_path[512];
		int retval = boinc_resolve_filename(checkpoint_file, input_path, sizeof(input_path));
		if (retval) {
			fprintf(stderr, "APP: error resolving hessian checkpoint file (for write): %d\n", retval);
			fprintf(stderr, "\tfilename: %s\n", checkpoint_file);
			fprintf(stderr, "\tresolved input path: %s\n", input_path);
			return retval;
		}

		file = boinc_fopen(input_path, "w");
	#else
		file = fopen(checkpoint_file, "w");
	#endif
	if (file == NULL) {
		fprintf(stderr, "APP: error reading hessian checkpoint file (for write): data_file == NULL\n");
		return 1;
	}

	fprintf(file, "n: %d, j: %d\n", number_parameters, j);
	fwrite_double_array(file, "gradient", number_parameters, gradient);
	fclose(file);
	return 0;
}

int read_gradient_checkpoint(char *checkpoint_file, int number_parameters, int *j, double *gradient) {
	FILE *file;
	#ifdef BOINC_APPLICATION 
		char input_path[512];
		int retval = boinc_resolve_filename(checkpoint_file, input_path, sizeof(input_path));
		if (retval) {
			fprintf(stderr, "APP: error resolving hessian checkpoint file (for write): %d\n", retval);
			fprintf(stderr, "\tfilename: %s\n", checkpoint_file);
			fprintf(stderr, "\tresolved input path: %s\n", input_path);
			return retval;
		}

		file = boinc_fopen(input_path, "r");
	#else
		file = fopen(checkpoint_file, "r");
	#endif
	if (file == NULL) {
		fprintf(stderr, "APP: error reading hessian checkpoint file (for write): data_file == NULL\n");
		return 1;
	}

	fscanf(file, "n: %d, j: %d\n", &number_parameters, j);
	fread_double_array__no_alloc(file, "gradient", number_parameters, gradient);
	fclose(file);
	return 0;
}

void get_gradient__checkpointed(int number_parameters, double *point, double *step, double *gradient, char *checkpoint_file) {
	int j;
	double e1, e2, pj;

	read_gradient_checkpoint(checkpoint_file, number_parameters, &j, gradient);

	for (j = 0; j < number_parameters; j++) {
		pj = point[j];
		point[j] = pj + step[j];
		e1 = evaluate(point);
		point[j] = pj - step[j];
		e2 = evaluate(point);
		point[j] = pj;

		gradient[j] = (e1 - e2)/(step[j] + step[j]);
		printf("\t\tgradient[%d]: %.20lf, (%.20lf - %.20lf)/(2 * %lf)\n", j, gradient[j], e1, e2, step[j]);
		
		#ifdef BOINC_APPLICATION
			if (boinc_time_to_checkpoint()) {
				checkpoint_gradient(checkpoint_file, number_parameters, j, gradient);
			}
		#else
			checkpoint_gradient(checkpoint_file, number_parameters, j, gradient);
		#endif
	}
}

void get_gradient(int number_parameters, double *point, double *step, double *gradient) {
	int j;
	double e1, e2, pj;
	for (j = 0; j < number_parameters; j++) {
		pj = point[j];
		point[j] = pj + step[j];
		e1 = evaluate(point);
		point[j] = pj - step[j];
		e2 = evaluate(point);
		point[j] = pj;

		gradient[j] = (e1 - e2)/(step[j] + step[j]);
		printf("\t\tgradient[%d]: %.20lf, (%.20lf - %.20lf)/(2 * %lf)\n", j, gradient[j], e1, e2, step[j]);
	}
}

int gradient_below_threshold(int number_parameters, double* gradient, double threshold) {
	int i;
	for (i = 0; i < number_parameters; i++) {
		if (gradient[i] > threshold || gradient[i] < -threshold) return 0;
	}
	return 1;
}
