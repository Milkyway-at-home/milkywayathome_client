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
         *      BOINC includes
****/
#ifdef BOINC_APPLICATION 
	#ifdef _WIN32
		#include "boinc_win.h"
	#else
		#include "config.h"
	#endif

	#ifndef _WIN32
		#include <cstdio>
		#include <cctype>
		#include <ctime>
		#include <cstring>
			#include <cstdlib>
		#include <csignal>
		#include <unistd.h>
	#endif

	#include "diagnostics.h"
	#include "util.h"
	#include "filesys.h"
	#include "boinc_api.h"
	#include "mfile.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/****
	*	Astronomy includes
*****/
#include "parameters.h"
#include "../util/io_util.h"

void free_parameters(ASTRONOMY_PARAMETERS* ap) {
	int i;
	free(ap->background_parameters);
	free(ap->background_step);
	free(ap->background_min);
	free(ap->background_max);
	free(ap->background_optimize);

	free(ap->stream_weights);
	free(ap->stream_weight_step);
	free(ap->stream_weight_min);
	free(ap->stream_weight_max);
	free(ap->stream_weight_optimize);

	for (i = 0; i < ap->number_streams; i++) {
		free(ap->stream_parameters[i]);
		free(ap->stream_step[i]);
		free(ap->stream_min[i]);
		free(ap->stream_max[i]);
		free(ap->stream_optimize[i]);
	}
	free(ap->stream_parameters);
	free(ap->stream_step);
	free(ap->stream_min);
	free(ap->stream_max);
	free(ap->stream_optimize);

	for (i = 0; i < ap->number_integrals; i++) {
		free(ap->integral[i]);
	}
	free(ap->integral);
}

int read_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS *ap) {
#ifdef BOINC_APPLICATION 
	char input_path[512];
	int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));

	if (retval) {
		fprintf(stderr, "APP: error resolving parameters file [%s], %d\n", filename, retval);
		return retval;
	}

	FILE* data_file = boinc_fopen(input_path, "r");
#else
	FILE* data_file = fopen(filename, "r");
#endif
	if (!data_file) {
		fprintf(stderr, "Couldn't find input file [%s] to read astronomy parameters.\n", filename);
		return 1;
	}

	fread_astronomy_parameters(data_file, ap);
	fclose(data_file);
	return 0;
}

int write_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS *ap) {
	FILE* data_file = fopen(filename, "w");
	if (!data_file) {
		fprintf(stderr, "Couldn't find output file [%s] to write astronomy parameters.\n", filename);
		return 1;
	}

	fwrite_astronomy_parameters(data_file, ap);
	fclose(data_file);
	return 0;
}

void fread_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS *ap) {
	int i;

	fscanf(file, "number_parameters: %d\n", &ap->number_background_parameters);
	fscanf(file, "background_weight: %lf\n", &ap->background_weight);

	fread_double_array(file, "background_parameters", &ap->background_parameters);
	fread_double_array(file, "background_step", &ap->background_step);
	fread_double_array(file, "background_min", &ap->background_min);
	fread_double_array(file, "background_max", &ap->background_max);
	fread_int_array(file, "optimize_parameter", &ap->background_optimize);

	fscanf(file, "number_streams: %d, %d\n", &ap->number_streams, &ap->number_stream_parameters);
	ap->stream_weights				= (double*)malloc(sizeof(double) * ap->number_streams);
	ap->stream_weight_step				= (double*)malloc(sizeof(double) * ap->number_streams);
	ap->stream_weight_min				= (double*)malloc(sizeof(double) * ap->number_streams);
	ap->stream_weight_max				= (double*)malloc(sizeof(double) * ap->number_streams);
	ap->stream_weight_optimize			=   (int*)malloc(sizeof(int)   * ap->number_streams);

	ap->stream_parameters				= (double**)malloc(sizeof(double*) * ap->number_streams);
	ap->stream_step					= (double**)malloc(sizeof(double*) * ap->number_streams);
	ap->stream_min					= (double**)malloc(sizeof(double*) * ap->number_streams);
	ap->stream_max					= (double**)malloc(sizeof(double*) * ap->number_streams);
	ap->stream_optimize				=   (int**)malloc(sizeof(int*)   * ap->number_streams);

	for (i = 0; i < ap->number_streams; i++) {
		fscanf(file, "stream_weight: %lf\n", &ap->stream_weights[i]);
		fscanf(file, "stream_weight_step: %lf\n", &ap->stream_weight_step[i]);
		fscanf(file, "stream_weight_min: %lf\n", &ap->stream_weight_min[i]);
		fscanf(file, "stream_weight_max: %lf\n", &ap->stream_weight_max[i]);
		fscanf(file, "optimize_weight: %d\n", &ap->stream_weight_optimize[i]);

		fread_double_array(file, "stream_parameters", &ap->stream_parameters[i]);
		fread_double_array(file, "stream_step", &ap->stream_step[i]);
		fread_double_array(file, "stream_min", &ap->stream_min[i]);
		fread_double_array(file, "stream_max", &ap->stream_max[i]);
		fread_int_array(file, "optimize_parameter", &ap->stream_optimize[i]);
	}

	fscanf(file, "convolve: %d\n", &ap->convolve);
	fscanf(file, "sgr_coordinates: %d\n", &ap->sgr_coordinates);
	fscanf(file, "aux_bg_profile: %d\n", &ap->aux_bg_profile); //vickej2_bg 
	fscanf(file, "wedge: %d\n", &ap->wedge);

	ap->integral = (INTEGRAL**)malloc(sizeof(INTEGRAL*));
	ap->integral[0] = (INTEGRAL*)malloc(sizeof(INTEGRAL));

	fscanf(file, "r[min,max,steps]: %lf, %lf, %d\n", &ap->integral[0]->r_min, &ap->integral[0]->r_max, &ap->integral[0]->r_steps);
	fscanf(file, "mu[min,max,steps]: %lf, %lf, %d\n", &ap->integral[0]->mu_min, &ap->integral[0]->mu_max, &ap->integral[0]->mu_steps);
	fscanf(file, "nu[min,max,steps]: %lf, %lf, %d\n", &ap->integral[0]->nu_min, &ap->integral[0]->nu_max, &ap->integral[0]->nu_steps);
	ap->integral[0]->r_step_size = (ap->integral[0]->r_max - ap->integral[0]->r_min)/(double)ap->integral[0]->r_steps;
	ap->integral[0]->mu_step_size = (ap->integral[0]->mu_max - ap->integral[0]->mu_min)/(double)ap->integral[0]->mu_steps;
	ap->integral[0]->nu_step_size = (ap->integral[0]->nu_max - ap->integral[0]->nu_min)/(double)ap->integral[0]->nu_steps;
	ap->integral[0]->min_calculation = 0;
	ap->integral[0]->max_calculation = ap->integral[0]->r_steps * ap->integral[0]->mu_steps * ap->integral[0]->nu_steps;

	fscanf(file, "number_cuts: %d\n", &ap->number_integrals);
	ap->number_integrals++;
	if (ap->number_integrals > 1) {
		ap->integral = (INTEGRAL**)realloc(ap->integral, sizeof(INTEGRAL*) * ap->number_integrals);
		for (i = 1; i < ap->number_integrals; i++) {
			int temp;
			ap->integral[i] = (INTEGRAL*)malloc(sizeof(INTEGRAL));
			fscanf(file, "r_cut[min,max,steps][%d]: %lf, %lf, %d\n", &temp, &ap->integral[i]->r_min, &ap->integral[i]->r_max, &ap->integral[i]->r_steps);
			fscanf(file, "mu_cut[min,max,steps][%d]: %lf, %lf, %d\n", &temp, &ap->integral[i]->mu_min, &ap->integral[i]->mu_max, &ap->integral[i]->mu_steps);
			fscanf(file, "nu_cut[min,max,steps][%d]: %lf, %lf, %d\n", &temp, &ap->integral[i]->nu_min, &ap->integral[i]->nu_max, &ap->integral[i]->nu_steps);
			ap->integral[i]->r_step_size = (ap->integral[i]->r_max - ap->integral[i]->r_min)/(double)ap->integral[i]->r_steps;
			ap->integral[i]->mu_step_size = (ap->integral[i]->mu_max - ap->integral[i]->mu_min)/(double)ap->integral[i]->mu_steps;
			ap->integral[i]->nu_step_size = (ap->integral[i]->nu_max - ap->integral[i]->nu_min)/(double)ap->integral[i]->nu_steps;
			ap->integral[i]->min_calculation = 0;
			ap->integral[i]->max_calculation = ap->integral[i]->r_steps * ap->integral[i]->mu_steps * ap->integral[i]->nu_steps;
		}
	}
}

void fwrite_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS *ap) {
	int i;

	fprintf(file, "number_parameters: %d\n", ap->number_background_parameters);
	fprintf(file, "background_weight: %lf\n", ap->background_weight);
	fwrite_double_array(file, "background_parameters", ap->number_background_parameters, ap->background_parameters);
	fwrite_double_array(file, "background_step", ap->number_background_parameters, ap->background_step);
	fwrite_double_array(file, "background_min", ap->number_background_parameters, ap->background_min);
	fwrite_double_array(file, "background_max", ap->number_background_parameters, ap->background_max);
	fwrite_int_array(file, "optimize_parameter", ap->number_background_parameters, ap->background_optimize);
	
	fprintf(file, "number_streams: %d, %d\n", ap->number_streams, ap->number_stream_parameters);
	for (i = 0; i < ap->number_streams; i++) {
		fprintf(file, "stream_weight: %lf\n", ap->stream_weights[i]);
		fprintf(file, "stream_weight_step: %lf\n", ap->stream_weight_step[i]);
		fprintf(file, "stream_weight_min: %lf\n", ap->stream_weight_min[i]);
		fprintf(file, "stream_weight_max: %lf\n", ap->stream_weight_max[i]);
		fprintf(file, "optimize_weight: %d\n", ap->stream_weight_optimize[i]);

		fwrite_double_array(file, "stream_parameters", ap->number_stream_parameters, ap->stream_parameters[i]);
		fwrite_double_array(file, "stream_step", ap->number_stream_parameters, ap->stream_step[i]);
		fwrite_double_array(file, "stream_min", ap->number_stream_parameters, ap->stream_min[i]);
		fwrite_double_array(file, "stream_max", ap->number_stream_parameters, ap->stream_max[i]);
		fwrite_int_array(file, "optimize_parameter", ap->number_stream_parameters, ap->stream_optimize[i]);
	}

        fprintf(file, "convolve: %d\n", ap->convolve);
	fprintf(file, "sgr_coordinates: %d\n", ap->sgr_coordinates);
        fprintf(file, "aux_bg_profile: %d\n", ap->aux_bg_profile); //vickej2_bg 
	fprintf(file, "wedge: %d\n", ap->wedge);

	fprintf(file, "r[min,max,steps]: %lf, %lf, %d\n", ap->integral[0]->r_min, ap->integral[0]->r_max, ap->integral[0]->r_steps);
	fprintf(file, "mu[min,max,steps]: %lf, %lf, %d\n", ap->integral[0]->mu_min, ap->integral[0]->mu_max, ap->integral[0]->mu_steps);
	fprintf(file, "nu[min,max,steps]: %lf, %lf, %d\n", ap->integral[0]->nu_min, ap->integral[0]->nu_max, ap->integral[0]->nu_steps);

	fprintf(file, "number_cuts: %d\n", ap->number_integrals-1);
	for (i = 1; i < ap->number_integrals; i++) {
		fprintf(file, "r_cut[min,max,steps][3]: %lf, %lf, %d\n", ap->integral[i]->r_min, ap->integral[i]->r_max, ap->integral[i]->r_steps);
		fprintf(file, "mu_cut[min,max,steps][3]: %lf, %lf, %d\n", ap->integral[i]->mu_min, ap->integral[i]->mu_max, ap->integral[i]->mu_steps);
		fprintf(file, "nu_cut[min,max,steps][3]: %lf, %lf, %d\n", ap->integral[i]->nu_min, ap->integral[i]->nu_max, ap->integral[i]->nu_steps);
	}
}

int get_optimized_parameter_count(ASTRONOMY_PARAMETERS *ap) {
	int i, j, count;
	count = 0;
	for (i = 0; i < ap->number_background_parameters; i++) {
		if (ap->background_optimize[i]) count++;
	}
	for (i = 0; i < ap->number_streams; i++) {
		if (ap->stream_weight_optimize[i]) count++;
		for (j = 0; j < ap->number_stream_parameters; j++) {
			if (ap->stream_optimize[i][j]) count++;
		}
	}
	return count;
}

void set_astronomy_parameters(ASTRONOMY_PARAMETERS *ap, double* parameters) {
	int i, j;
	int current;
	current = 0;
	for (i = 0; i < ap->number_background_parameters; i++) {
		if (ap->background_optimize[i]) {
			ap->background_parameters[i] = parameters[current];
			current++;
		}
	}

	for (i = 0; i < ap->number_streams; i++) {
		if (ap->stream_weight_optimize[i]) {
			ap->stream_weights[i] = parameters[current];
			current++;
		}

		for (j = 0; j < ap->number_stream_parameters; j++) {
			if (ap->stream_optimize[i][j]) {
				ap->stream_parameters[i][j] = parameters[current];
				current++;
			}
		}
	}
}

void get_search_parameters(ASTRONOMY_PARAMETERS *ap, double** result) {
	int i, j, current;
	current = 0;
	(*result) = (double*)malloc(sizeof(double) * get_optimized_parameter_count(ap));

	for (i = 0; i < ap->number_background_parameters; i++) {
		if (ap->background_optimize[i]) {
			(*result)[current] = ap->background_parameters[i];
			current++;
		}
	}

	for (i = 0; i < ap->number_streams; i++) {
		if (ap->stream_weight_optimize[i]) {
			(*result)[current] = ap->stream_weights[i];
			current++;
		}

		for (j = 0; j < ap->number_stream_parameters; j++) {
			if (ap->stream_optimize[i][j]) {
				(*result)[current] = ap->stream_parameters[i][j];
				current++;
			}
		}
	}
}

void get_step(ASTRONOMY_PARAMETERS *ap, double** result) {
	int i, j, current;
	current = 0;
	(*result) = (double*)malloc(sizeof(double) * get_optimized_parameter_count(ap));

	for (i = 0; i < ap->number_background_parameters; i++) {
		if (ap->background_optimize[i]) {
			(*result)[current] = ap->background_step[i];
			current++;
		}
	}

	for (i = 0; i < ap->number_streams; i++) {
		if (ap->stream_weight_optimize[i]) {
			(*result)[current] = ap->stream_weight_step[i];
			current++;
		}

		for (j = 0; j < ap->number_stream_parameters; j++) {
			if (ap->stream_optimize[i][j]) {
				(*result)[current] = ap->stream_step[i][j];
				current++;
			}
		}
	}
}


void get_min_parameters(ASTRONOMY_PARAMETERS *ap, double** result) {
	int i, j, current;
	current = 0;
	(*result) = (double*)malloc(sizeof(double) * get_optimized_parameter_count(ap));

	for (i = 0; i < ap->number_background_parameters; i++) {
		if (ap->background_optimize[i]) {
			(*result)[current] = ap->background_min[i];
			current++;
		}
	}

	for (i = 0; i < ap->number_streams; i++) {
		if (ap->stream_weight_optimize[i]) {
			(*result)[current] = ap->stream_weight_min[i];
			current++;
		}

		for (j = 0; j < ap->number_stream_parameters; j++) {
			if (ap->stream_optimize[i][j]) {
				(*result)[current] = ap->stream_min[i][j];
				current++;
			}
		}
	}
}

void get_max_parameters(ASTRONOMY_PARAMETERS *ap, double** result) {
	int i, j, current;
	current = 0;
	(*result) = (double*)malloc(sizeof(double) * get_optimized_parameter_count(ap));

	for (i = 0; i < ap->number_background_parameters; i++) {
		if (ap->background_optimize[i]) {
			(*result)[current] = ap->background_max[i];
			current++;
		}
	}

	for (i = 0; i < ap->number_streams; i++) {
		if (ap->stream_weight_optimize[i]) {
			(*result)[current] = ap->stream_weight_max[i];
			current++;
		}

		for (j = 0; j < ap->number_stream_parameters; j++) {
			if (ap->stream_optimize[i][j]) {
				(*result)[current] = ap->stream_max[i][j];
				current++;
			}
		}
	}
}

void split_astronomy_parameters(ASTRONOMY_PARAMETERS *ap, int rank, int max_rank) {
	long total_calculations;
	int i;

	for (i = 0; i < ap->number_integrals; i++) {
		total_calculations = ap->integral[i]->r_steps * ap->integral[i]->nu_steps * ap->integral[i]->mu_steps;
		ap->integral[i]->min_calculation = (long) (total_calculations * ((double)rank / (double)max_rank));
		ap->integral[i]->max_calculation = (long) (total_calculations * ((double)(rank + 1) / (double)max_rank));

		printf("[worker: %d] [integral: %d] min_calculation: %ld / max_calculation: %ld, total_calculations: %ld\n", rank, i, ap->integral[i]->min_calculation, ap->integral[i]->max_calculation, total_calculations);
	}	
}


#ifdef ASTRONOMY_PARAMETERS_MAIN

int main(int argc, char **argv) {
	ASTRONOMY_PARAMETERS *ap;

	printf("in: %s, out: %s\n", argv[1], argv[2]);

	ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	read_astronomy_parameters(argv[1], ap);
	fwrite_astronomy_parameters(stdout, ap);
	write_astronomy_parameters(argv[2], ap);

	return 0;
}

#endif


#ifdef ASTRONOMY_PARAMETERS_PERMUTE
#include <time.h>

int main(int argc, char **argv) {
	ASTRONOMY_PARAMETERS *ap;
	double *point, *min_bound, *max_bound, range, current_range, test;
	int i, number_parameters;
	time_t t1;

	time(&t1);

	range = atof(argv[3]);
	printf("in: %s, out: %s, range: %lf\n", argv[1], argv[2], range);

	srand48((long) t1);

	ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	read_astronomy_parameters(argv[1], ap);
	fwrite_astronomy_parameters(stdout, ap);

	printf("\n\n");

	number_parameters = get_optimized_parameter_count(ap);
	get_search_parameters(ap, &point);
	get_min_parameters(ap, &min_bound);
	get_max_parameters(ap, &max_bound);

	for (i = 0; i < number_parameters; i++) {
		current_range = (max_bound[i] - min_bound[i]) * range;
		test = 0.0;
		while (test == 0.0 || test > max_bound[i] || test < min_bound[i]) {
			test = point[i] + ((drand48() * current_range * 2.0) - current_range);
		}
		point[i] = test;
	}
	set_astronomy_parameters(ap, point);

	fwrite_astronomy_parameters(stdout, ap);
	write_astronomy_parameters(argv[2], ap);

	return 0;
}

#endif
