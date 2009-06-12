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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "search_parameters.h"
#include "../util/settings.h"
#include "../util/matrix.h"
#include "../util/io_util.h"

/****
         *     BOINC includes
*****/
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

	using std::string;
#endif


int fwrite_gpu_result(FILE *file, int number_parameters, double **hessian, double *gradient, double initial_fitness, double *initial_parameters, double result_fitness, double *result_parameters, int number_evaluations, char *metadata) {
	fwrite_matrix(file, "hessian", hessian, number_parameters, number_parameters);
	fwrite_double_array(file, "gradient", number_parameters, gradient);
	fprintf(file, "initial_fitness: %.20lf\n", initial_fitness);
	fwrite_double_array(file, "inital_parameters", number_parameters, initial_parameters);
	fprintf(file, "result_fitness: %.20lf\n", result_fitness);
	fwrite_double_array(file, "result_parameters", number_parameters, result_parameters);
	fprintf(file, "number_evaluations: %d\n", number_evaluations);
	fprintf(file, "metadata: %s\n", metadata);
	return 0;
}

int write_gpu_result(char *filename, int number_parameters, double **hessian, double *gradient, double initial_fitness, double *initial_parameters, double result_fitness, double *result_parameters, int number_evaluations, char *metadata) {
#ifdef BOINC_APPLICATION
	char output_path[512];
	int retval = boinc_resolve_filename(filename, output_path, sizeof(output_path));
	if (retval) {
		fprintf(stderr, "APP: error resolving search parameters file (for write): %d\n", retval);
		fprintf(stderr, "\tfilename: %s\n", filename);
		fprintf(stderr, "\tresolved output path: %s\n", output_path);
		return retval;
	}

	FILE* data_file = boinc_fopen(output_path, "w");
#else
	FILE* data_file = fopen(filename, "w");
#endif
	return fwrite_gpu_result(data_file, number_parameters, hessian, gradient, initial_fitness, initial_parameters, result_fitness, result_parameters, number_evaluations, metadata);
}

int fread_gpu_result(FILE *file, int number_parameters, double **hessian, double *gradient, double *initial_fitness, double *initial_parameters, double *result_fitness, double *result_parameters, int *number_evaluations, char *metadata) {
	int retval;
	retval = fread_matrix(file, "hessian", hessian, number_parameters, number_parameters);
	if (retval == MATRIX__READ_ERROR) {
		fprintf(stderr, "error reading hessian from GPU result: [%s]\n", MATRIX__ERROR_MSG);
		return 1;
	}

	retval = fread_double_array__no_alloc(file, "gradient", number_parameters, gradient);
	if (retval < number_parameters) {
		fprintf(stderr, "error reading gradient from GPU result: %d\n", retval);
		return 1;
	}

	retval = fscanf(file, "initial_fitness: %lf\n", initial_fitness);
	if (retval != 1) {
		fprintf(stderr, "error reading initial fitness from GPU result\n");
		return 1;
	}

	retval = fread_double_array__no_alloc(file, "initial_parameters", number_parameters, initial_parameters);
	if (retval < number_parameters) {
		fprintf(stderr, "error reading initial parameters from GPU result: %d\n", retval);
		return 1;
	}

	retval = fscanf(file, "result_fitness: %lf\n", result_fitness);
	if (retval != 1) {
		fprintf(stderr, "error reading result fitness from GPU result\n");
		return 1;
	}

	retval = fread_double_array__no_alloc(file, "result_parameters", number_parameters, result_parameters);
	if (retval < number_parameters) {
		fprintf(stderr, "error reading result parameters from GPU result: %d\n", retval);
		return 1;
	}

	retval = fscanf(file, "number_evaluations: %d\n", number_evaluations);
	if (retval != 1) {
		fprintf(stderr, "error reading number evaluations from GPU result\n");
		return 1;
	}

	retval = fread_metadata(file, metadata);
	if (retval != 0) {
		fprintf(stderr, "error reading metadata from GPU result\n");
		return 1;
	}
	return 0;
}

int read_gpu_result(char *filename, int number_parameters, double **hessian, double *gradient, double *initial_fitness, double *initial_parameters, double *result_fitness, double *result_parameters, int *number_evaluations, char *metadata) {
	FILE* data_file = fopen(filename, "r");
	if (data_file == NULL) {
		fprintf(stderr, "APP: error reading gpu result file: data_file == NULL\n");
		return 1;
	}
	return fread_gpu_result(data_file, number_parameters, hessian, gradient, initial_fitness, initial_parameters, result_fitness, result_parameters, number_evaluations, metadata);
}

int fread_cpu_result(FILE *file, int number_parameters, double *parameters, double *fitness, char *metadata) {
	fread_double_array__no_alloc(file, "parameters", number_parameters, parameters);
	fscanf(file, "fitness: %lf\n", fitness);
	fread_metadata(file, metadata);
	return 0;
}

int read_cpu_result(char *filename, int number_parameters, double *parameters, double *fitness, char *metadata) {
	int retval;
	FILE *data_file = fopen(filename, "r");
	if (data_file == NULL) {
		fprintf(stderr, "APP: error reading cpu result file: data_file == NULL\n");
		return 1;
	}
	retval = fread_cpu_result(data_file, number_parameters, parameters, fitness, metadata);
	fclose(data_file);
	return retval;
}

int fread_cpu_result__realloc(FILE *file, char *search_name, int *number_parameters, double **parameters, double *fitness, char *metadata, char *app_version) {
	int i, count;

	fscanf(file, "%s\n", search_name);
	fread_double_array__realloc(file, "parameters", number_parameters, parameters);
	fread_metadata(file, metadata);
	fscanf(file, "fitness: %lf\n", fitness);

	memset(app_version, '\0', sizeof(char) * 128);
	if (fgets(app_version, 128, file) == NULL || strlen(app_version) < 5) {
		sprintf(app_version, "?");
	} else {
		count = 0;
		for (i = strlen(app_version); i >= 0; i--) {
			if (app_version[i] == '\0') {
			} else if (app_version[i] == ' ' || app_version[i] == 10 || app_version[i] == 13) {
				app_version[i] = '\0';
//				count++;
			} else {
//				printf("breaked after count: %d\n", count);
				break;
			}
		}
	}

	return 0;
}

int read_cpu_result__realloc(const char *filename, char *search_name, int *number_parameters, double **parameters, double *fitness, char *metadata, char *app_version) {
	int retval;
	FILE *data_file = fopen(filename, "r");
	if (data_file == NULL) {
		fprintf(stderr, "APP: error reading cpu result file: data_file == NULL\n");
		return 1;
	}
	retval = fread_cpu_result__realloc(data_file, search_name, number_parameters, parameters, fitness, metadata, app_version);
	fclose(data_file);
	return retval;
}

int fwrite_cpu_result(FILE *file, char *search_name, int number_parameters, double *parameters, double fitness, char *metadata, char *app_version, char *precision) {
	fprintf(file, "%s\n", search_name);
	fwrite_double_array(file, "parameters", number_parameters, parameters);
	fprintf(file, "metadata: %s\n", metadata);
	fprintf(file, "fitness: %.20lf\n", fitness);
	fprintf(file, "%s %s\n", app_version, precision);
	return 0;
}

int write_cpu_result(char *filename, char *search_name, int number_parameters, double *parameters, double fitness, char *metadata, char *app_version, char *precision) {
#ifdef BOINC_APPLICATION
	char output_path[512];
	int retval = boinc_resolve_filename(filename, output_path, sizeof(output_path));
	if (retval) {
		fprintf(stderr, "APP: error resolving search parameters file (for write): %d\n", retval);
		fprintf(stderr, "\tfilename: %s\n", filename);
		fprintf(stderr, "\tresolved output path: %s\n", output_path);
		return retval;
	}

	FILE* data_file = boinc_fopen(output_path, "w");
#else
	FILE* data_file = fopen(filename, "w");
#endif
	return fwrite_cpu_result(data_file, search_name, number_parameters, parameters, fitness, metadata, app_version, precision); 
}
