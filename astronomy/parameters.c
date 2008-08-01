/*
 *  parameters.c
 *  Astronomy
 *
 *  Created by Travis Desell on 2/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

/****
         *      BOINC includes
****/
#ifdef GMLE_BOINC
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
}

int read_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS *ap) {
	FILE* data_file = fopen(filename, "r");
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
	int i, j;

	fscanf(file, "number_parameters: %d\n", &ap->number_background_parameters);
	fscanf(file, "background_weight: %lf\n", &ap->background_weight);

	read_double_array(file, "background_parameters: ", ap->number_background_parameters, &ap->background_parameters);
	read_double_array(file, "background_step: ", ap->number_background_parameters, &ap->background_step);
	read_double_array(file, "background_min: ", ap->number_background_parameters, &ap->background_min);
	read_double_array(file, "background_max: ", ap->number_background_parameters, &ap->background_max);
	read_int_array(file, "optimize_parameter: ", ap->number_background_parameters, &ap->background_optimize);

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

		read_double_array(file, "stream_parameters: ", ap->number_stream_parameters, &ap->stream_parameters[i]);
		read_double_array(file, "stream_step: ", ap->number_stream_parameters, &ap->stream_step[i]);
		read_double_array(file, "stream_min: ", ap->number_stream_parameters, &ap->stream_min[i]);
		read_double_array(file, "stream_max: ", ap->number_stream_parameters, &ap->stream_max[i]);
		read_int_array(file, "optimize_parameter: ", ap->number_stream_parameters, &ap->stream_optimize[i]);
	}

	fscanf(file, "convolve: %d\n", &ap->convolve);
	fscanf(file, "wedge: %d\n", &ap->wedge);
	fscanf(file, "r[min,max,steps]: %lf, %lf, %d\n", &ap->r_min, &ap->r_max, &ap->r_steps);
	fscanf(file, "mu[min,max,steps]: %lf, %lf, %d\n", &ap->mu_min, &ap->mu_max, &ap->mu_steps);
	fscanf(file, "nu[min,max,steps]: %lf, %lf, %d\n", &ap->nu_min, &ap->nu_max, &ap->nu_steps);
	ap->r_step_size = (ap->r_max - ap->r_min)/(double)ap->r_steps;
	ap->mu_step_size = (ap->mu_max - ap->mu_min)/ap->mu_steps;
	ap->nu_step_size = (ap->nu_max - ap->nu_min)/ap->nu_steps;

	ap->number_parameters = 0;
	for (i = 0; i < ap->number_background_parameters; i++) {
		if (ap->background_optimize[i]) {
			ap->number_parameters++;
		}
	}

	for (i = 0; i < ap->number_streams; i++) {
		if (ap->stream_weight_optimize[i]) {
			ap->number_parameters++;
		}

		for (j = 0; j < ap->number_stream_parameters; j++) {
			if (ap->stream_optimize[i][j]) {
				ap->number_parameters++;
			}
		}
	}
}

void fwrite_astronomy_parameters(FILE* file, ASTRONOMY_PARAMETERS *ap) {
	int i;

	fprintf(file, "number_parameters: %d\n", ap->number_background_parameters);
	fprintf(file, "background_weight: %lf\n", ap->background_weight);
	print_double_array(file, "background_parameters: ", ap->number_background_parameters, ap->background_parameters);
	print_double_array(file, "background_step: ", ap->number_background_parameters, ap->background_step);
	print_double_array(file, "background_min: ", ap->number_background_parameters, ap->background_min);
	print_double_array(file, "background_max: ", ap->number_background_parameters, ap->background_max);
	print_int_array(file, "optimize_parameter: ", ap->number_background_parameters, ap->background_optimize);
	
	fprintf(file, "number_streams: %d, %d\n", ap->number_streams, ap->number_stream_parameters);
	for (i = 0; i < ap->number_streams; i++) {
		fprintf(file, "stream_weight: %lf\n", ap->stream_weights[i]);
		fprintf(file, "stream_weight_step: %lf\n", ap->stream_weight_step[i]);
		fprintf(file, "stream_weight_min: %lf\n", ap->stream_weight_min[i]);
		fprintf(file, "stream_weight_max: %lf\n", ap->stream_weight_max[i]);
		fprintf(file, "optimize_weight: %d\n", ap->stream_weight_optimize[i]);

		print_double_array(file, "stream_parameters: ", ap->number_stream_parameters, ap->stream_parameters[i]);
		print_double_array(file, "stream_step: ", ap->number_stream_parameters, ap->stream_step[i]);
		print_double_array(file, "stream_min: ", ap->number_stream_parameters, ap->stream_min[i]);
		print_double_array(file, "stream_max: ", ap->number_stream_parameters, ap->stream_max[i]);
		print_int_array(file, "optimize_parameter: ", ap->number_stream_parameters, ap->stream_optimize[i]);
	}

        fprintf(file, "convolve: %d\n", ap->convolve);
        fprintf(file, "wedge: %d\n", ap->wedge);
	fprintf(file, "r[min,max,steps]: %lf, %lf, %d\n", ap->r_min, ap->r_max, ap->r_steps);
	fprintf(file, "mu[min,max,steps]: %lf, %lf, %d\n", ap->mu_min, ap->mu_max, ap->mu_steps);
	fprintf(file, "nu[min,max,steps]: %lf, %lf, %d\n", ap->nu_min, ap->nu_max, ap->nu_steps);
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
	(*result) = (double*)malloc(sizeof(double) * ap->number_parameters);

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
	(*result) = (double*)malloc(sizeof(double) * ap->number_parameters);

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
	(*result) = (double*)malloc(sizeof(double) * ap->number_parameters);

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
	(*result) = (double*)malloc(sizeof(double) * ap->number_parameters);

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



void split(int rank, int max_rank, int divisor, double *min, double *max, int *steps, double step) {
	int sub_rank, next_sub_rank, extra;

	printf("rank: %d, max_rank: %d, divisor: %d, min: %lf, max: %lf, steps: %d, step: %lf\n", rank, max_rank, divisor, (*min), (*max), (*steps), step);
	sub_rank = rank % divisor;
	next_sub_rank = sub_rank + 1;
	extra = (*steps) % divisor;

	(*steps) = (*steps) / divisor;

	int min_slice, max_slice;
	int difference = -1;
	if (sub_rank < extra) {
		(*steps) = (*steps) + 1;
		min_slice = sub_rank * (*steps);
		max_slice = next_sub_rank * (*steps);
	} else {
		difference = sub_rank - extra;

		min_slice = extra * ((*steps) + 1) + difference * (*steps);
		max_slice = extra * ((*steps) + 1) + (difference+1) * (*steps);
	}
	(*max) = (*min) + max_slice * step;
	(*min) = (*min) + min_slice * step;

//	printf("splitting: sub_rank: %d, next_sub_rank: %d, extra: %d, difference: %d, min_slice: %d, max_slice: %d, [rank: %d, max_rank: %d]\n", sub_rank, next_sub_rank, extra, difference, min_slice, max_slice, rank, max_rank);
}

void split_astronomy_parameters(ASTRONOMY_PARAMETERS *ap, int rank, int max_rank) {
	int r_divisor = max_rank;
	int mu_divisor = 1;
	int nu_divisor = 1;

	if (rank == 0 && max_rank == 0) return;

	while (r_divisor > ap->r_steps) {
		if (r_divisor > ap->r_steps && (r_divisor % 2) != 0) {
			printf("ERROR cannot split data,  slices > ap->r_steps and not divisible by 2");
		}

		r_divisor /= 2;
		mu_divisor *= 2;
	}
	while (mu_divisor > ap->mu_steps) {
		if (mu_divisor > ap->mu_steps && (mu_divisor % 2) != 0) {
			printf("ERROR cannot split data,  slices > ap->r_steps and not divisible by 2");
		}

		mu_divisor /= 2;
		nu_divisor *= 2;
	}

//	printf("r_divisor: %d, mu_divisor: %d, nu_divisor: %d\n", r_divisor, mu_divisor, nu_divisor);

	printf("splitting ap->r\n");
	/********
		*	Split ap->r
	 ********/
	split(rank, max_rank, r_divisor, &(ap->r_min), &(ap->r_max), &(ap->r_steps), ap->r_step_size);


	printf("splitting ap->mu\n");
	/********
		*	Split ap->mu
	 ********/
	split((rank/r_divisor), max_rank, mu_divisor, &(ap->mu_min), &(ap->mu_max), &(ap->mu_steps), ap->mu_step_size);


	printf("splitting ap->nu\n");
	/********
		*	Split ap->nu
	 ********/
	split(rank/(r_divisor*mu_divisor), max_rank, nu_divisor, &(ap->nu_min), &(ap->nu_max), &(ap->nu_steps), ap->nu_step_size);

//	printf("div: [%d %d %d] rank: [%d / %d] [r_min: %lf, r_max: %lf, r_steps: %d], [mu_min: %lf, mu_max: %lf, mu_steps: %d], [nu_min: %lf, nu_max: %lf, nu_steps: %d]\n", r_divisor, mu_divisor, nu_divisor, rank, max_rank, ap->r_min, ap->r_max, ap->r_steps, ap->mu_min, ap->mu_max, ap->mu_steps, ap->nu_min, ap->nu_max, ap->nu_steps);
}

#ifdef GMLE_BOINC
	int boinc_read_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS *ap) {
	char input_path[512];
		int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));

		if (retval) {
			fprintf(stderr, "APP: error resolving parameters file [%s], %d\n", filename, retval);
			return retval;
		}

		FILE* data_file = boinc_fopen(input_path, "r");
		if (!data_file) {
			fprintf(stderr, "Couldn't find input file [%s] to read astronomy parameters.\n", filename);
			return 1;
		}

		fread_astronomy_parameters(data_file, ap);
		fclose(data_file);
		return 0;
	}

	int boinc_write_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS *ap) {
		char input_path[512];
		int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));

		if (retval) {
			fprintf(stderr, "APP: error writing astronomy parameters [%s], %d\n", filename, retval);
			return retval;
		}

		FILE* data_file = boinc_fopen(input_path, "w");
		if (!data_file) {
			fprintf(stderr, "Couldn't find output file [%s] to write astronomy parameters.\n", filename);
			return 1;
		}

		fwrite_astronomy_parameters(data_file, ap);
		fclose(data_file);
		return 0;
	}

#endif
