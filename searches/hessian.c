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

#include "hessian.h"
#include "../evaluation/evaluator.h"
#include "../util/settings.h"
#include "../util/matrix.h"
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

int checkpoint_hessian(char *checkpoint_file, int number_parameters, int i, int j, double **hessian) {
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

	fprintf(file, "n: %d, i: %d, j: %d\n", number_parameters, i, j);
	fwrite_matrix(file, "hessian", hessian, number_parameters, number_parameters);
	fclose(file);

	return 0;
}

int read_hessian_checkpoint(char *checkpoint_file, int *number_parameters, int *i, int *j, double **hessian) {
	FILE *file;
	#ifdef BOINC_APPLICATION 
		char input_path[512];
		int retval = boinc_resolve_filename(checkpoint_file, input_path, sizeof(input_path));
		if (retval) {
			fprintf(stderr, "APP: error resolving hessian checkpoint file (for read): %d\n", retval);
			fprintf(stderr, "\tfilename: %s\n", checkpoint_file);
			fprintf(stderr, "\tresolved input path: %s\n", input_path);
			return retval;
		}

		file = boinc_fopen(input_path, "r");
	#else
		file = fopen(checkpoint_file, "r");
	#endif
	if (file == NULL) {
		fprintf(stderr, "APP: error reading hessian checkpoint file (for read): data_file == NULL\n");
		return 1;
	}

	fscanf(file, "n: %d, i: %d, j: %d\n", number_parameters, i, j);
	fread_matrix(file, "hessian", hessian, *number_parameters, *number_parameters);
	fclose(file);

	return 0;
}

void get_hessian__checkpointed(int number_parameters, double *point, double *step, double **hessian, char *checkpoint_file) {
	int i, j, np;
	double e1, e2, e3, e4;
	double pi, pj;

	read_hessian_checkpoint(checkpoint_file, &np, &i, &j, hessian);

	for (i = 0; i < number_parameters; i++) {
		pi = point[i];
		pj = point[j];
		point[i] = pi + step[i] + step[i];
		e1 = evaluate(point);
		point[i] = pi;
		e2 = e3 = evaluate(point);
		point[i] = pi - (step[i] + step[i]); 
		e4 = evaluate(point);
		point[i] = pi;
		point[j] = pj;

		hessian[i][i] = (e1 - e3 - e2 + e4)/(4 * step[i] * step[i]);
		printf("\t\thessian[%d][%d] = %.20lf, (%.20lf - %.20lf - %.20lf + %.20lf)/(4 * %.20lf * %.20lf)\n", i, i, hessian[i][i], e1, e3, e2, e4, step[i], step[j]);

		for (j = i+1; j < number_parameters; j++) {
			pi = point[i];
			pj = point[j];
			point[i] = pi + step[i];
			point[j] = pj + step[j];
			e1 = evaluate(point);

			point[i] = pi - step[i];
			e2 = evaluate(point);

			point[i] = pi + step[i];
			point[j] = pj - step[j];
			e3 = evaluate(point);

			point[i] = pi - step[i];
			e4 = evaluate(point);
			point[i] = pi;
			point[j] = pj;

			hessian[i][j] = (e1 - e3 - e2 + e4)/(4 * step[i] * step[j]);
			hessian[j][i] = hessian[i][j];
			printf("\t\thessian[%d][%d] = hessian[%d][%d] = %.20lf, (%.20lf - %.20lf - %.20lf + %.20lf)/(4 * %.20lf * %.20lf)\n", i, j, j, i, hessian[i][j], e1, e3, e2, e4, step[i], step[j]);

			#ifdef BOINC_APPLICATION
				if (boinc_time_to_checkpoint()) {
					checkpoint_hessian(checkpoint_file, number_parameters, i, j, hessian);
				}
			#else
				checkpoint_hessian(checkpoint_file, number_parameters, i, j, hessian);
			#endif
		}
	}
}

void get_hessian(int number_parameters, double *point, double *step, double **hessian) {
	int i, j;
	double e1, e2, e3, e4;
	double pi, pj;

	for (i = 0; i < number_parameters; i++) {
		for (j = 0; j < number_parameters; j++) {
			pi = point[i];
			pj = point[j];
			if (i == j) {
				point[i] = pi + step[i] + step[i];
				e1 = evaluate(point);
				point[i] = pi;
				e2 = e3 = evaluate(point);
				point[i] = pi - (step[i] + step[i]); 
				e4 = evaluate(point);
			} else {
				point[i] = pi + step[i];
				point[j] = pj + step[j];
				e1 = evaluate(point);

				point[i] = pi - step[i];
				e2 = evaluate(point);

				point[i] = pi + step[i];
				point[j] = pj - step[j];
				e3 = evaluate(point);

				point[i] = pi - step[i];
				e4 = evaluate(point);
			}
			point[i] = pi;
			point[j] = pj;

			hessian[i][j] = (e1 - e3 - e2 + e4)/(4 * step[i] * step[j]);
			printf("\t\thessian[%d][%d]: %.20lf\n", i, j, hessian[i][j]);
		}
	}
}

void randomized_hessian(double** actual_points, double* center, double* fitness, int number_points, int number_parameters, double*** hessian, double** gradient) {
	double** Y;
	double** X;
	double** X2;
	double** X3;
	double** X_transpose;
	double** X_inverse;
	double** W;
	double** points;
	int x_len, i, j, k, current;

	points = (double**)malloc(sizeof(double*) * number_points);
	for (i = 0; i < number_points; i++) {
		points[i] = (double*)malloc(sizeof(double) * number_parameters);
		for (j = 0; j < number_parameters; j++) {
			points[i][j] = center[j] - actual_points[i][j];
		}
	}


	/********
		*	X = [1, x1, ... xn, 0.5*x1^2, ... 0.5*xn^2, x1*x2, ..., x1*xn, x2*x3, ..., x2*xn, ...]
	 ********/
	x_len = 1 + number_parameters + number_parameters;
	for (i = number_parameters - 1; i > 0; i--) x_len += i;


	Y = (double**)malloc(sizeof(double*) * number_points);
	X = (double**)malloc(sizeof(double) * number_points);
	for (i = 0; i < number_points; i++) {
		Y[i] = (double*)malloc(sizeof(double) * 1);
		Y[i][0] = fitness[i];
		X[i] = (double*)malloc(sizeof(double) * x_len);
		X[i][0] = 1;
		for (j = 0; j < number_parameters; j++) {
			X[i][1+j] = points[i][j];
			X[i][1+number_parameters+j] = 0.5 * points[i][j] * points[i][j];
		}
		current = 0;
		for (j = 0; j < number_parameters; j++) {
			for (k = j+1; k < number_parameters; k++) {
				X[i][1+number_parameters+number_parameters+current] = points[i][j] * points[i][k];
				current++;
			}
		}
	}

        matrix_transpose__alloc(X, number_points, x_len, &X_transpose);
        matrix_multiply__alloc(X_transpose, x_len, number_points, X, number_points, x_len, &X2);
	matrix_invert__alloc(X2, x_len, x_len, &X_inverse);
        matrix_multiply__alloc(X_inverse, x_len, x_len, X_transpose, x_len, number_points, &X3);
	matrix_multiply__alloc(X3, x_len, number_points, Y, number_points, 1, &W);


	(*gradient) = (double*)malloc(sizeof(double) * number_parameters);
	(*hessian) = (double**)malloc(sizeof(double*) * number_parameters);
	for (i = 0; i < number_parameters; i++) (*hessian)[i] = (double*)malloc(sizeof(double) * number_parameters);

	for (i = 0; i < number_parameters; i++) {
		(*gradient)[i] = W[1+i][0];
		(*hessian)[i][i] = W[1 + number_parameters + i][0];

		current = 0;
		for (j = i; j < number_parameters; j++) {
			for (k = j+1; k < number_parameters; k++) {
				(*hessian)[j][k] = W[1 + number_parameters + number_parameters + current][0];
				(*hessian)[k][j] = W[1 + number_parameters + number_parameters + current][0];
				current++;
			}
		}
	}

	/********
		*	Y = number_points by 1
		*	X = number_points by x_len
		*	X2 = x_len by x_len
		*	X3 = x_len by number_points
		*	X_transpose = x_len by number_points
		*	X_inverse = x_len by x_len
		*	W = x_len by 1
	 ********/
	for (i = 0; i < number_points; i++) {
		free(Y[i]);
		free(X[i]);
		free(points[i]);
	}
	for (i = 0; i < x_len; i++) {
		free(X2[i]);
		free(X3[i]);
		free(X_transpose[i]);
		free(X_inverse[i]);
		free(W[i]);
	}
	free(points);
	free(Y);
	free(X);
	free(X2);
	free(X3);
	free(X_transpose);
	free(X_inverse);
	free(W);
}
