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

#include <stdio.h>
#include <stdlib.h>

void init_matrix(double ***m, int r, int c) {
	int i;
	(*m) = (double**)malloc(sizeof(double*) * r);
	for (i = 0; i < r; i++) (*m)[i] = (double*)malloc(sizeof(double) * c);
}

void free_matrix(double ***m, int r, int c) {
	int i;
	for (i = 0; i < r; i++) free( (*m)[i] );
	free( (*m) );
}


void matrix_print(FILE* file, const char* name, double** m, int r, int c) {
	int i, j;
	fprintf(file, "%s [%d x %d]:\n", name, r, c);
	for (i = 0; i < r; i++) {
		for (j = 0; j < c; j++) {
			fprintf(file, " %.10lf", m[i][j]);
		}
		fprintf(file, "\n");
	}
}


void matrix_transpose(double** m, int r, int c, double*** result) {
	int i, j;

	(*result) = (double**)malloc(sizeof(double*) * c);

	for (i = 0; i < c; i++) {
		(*result)[i] = (double*)malloc(sizeof(double) * r);
		for (j = 0; j < r; j++) {
			(*result)[i][j] = m[j][i];
		}
	}
}

void matrix_multiply(double** m1, int r1, int c1, double** m2, int r2, int c2, double*** result) {
	int i, j, k;
	(*result) = (double**)malloc(sizeof(double) * r1);
	for (i = 0; i < r1; i++) {
		(*result)[i] = (double*)malloc(sizeof(double) * c2);
		for (j = 0; j < c2; j++) {
			(*result)[i][j] = 0;
			for (k = 0; k < c1; k++) {
				(*result)[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}

void matrix_vector_multiply(double** m2, int r2, int c2, double* m1, int c1, double** result) {
	int j, k;

	(*result) = (double*)malloc(sizeof(double) * r2);
	for (j = 0; j < r2; j++) {
		(*result)[j] = 0;
		for (k = 0; k < c1; k++) {
			(*result)[j] += m1[k] * m2[j][k];
		}
	}
}

void vector_matrix_multiply(double* m1, int c1, double** m2, int r2, int c2, double** result) {
	int j, k;

	(*result) = (double*)malloc(sizeof(double) * r2);
	for (j = 0; j < r2; j++) {
		(*result)[j] = 0;
		for (k = 0; k < c1; k++) {
			(*result)[j] += m1[k] * m2[k][j];
		}
	}
}
 
void matrix_invert(double** initial, int rows, int cols, double*** inverse) {
	int i, j, next;
	double** matrix;
	double* temp;
	double working_value;

	matrix = (double**)malloc(sizeof(double*) * rows);
	(*inverse) = (double**)malloc(sizeof(double*) * cols);
	for (i = 0; i < rows; i++) {
		matrix[i] = (double*)malloc(sizeof(double) * cols);
		(*inverse)[i] = (double*)malloc(sizeof(double) * rows);
		for (j = 0; j < cols; j++) {
			matrix[i][j] = initial[i][j];
			(*inverse)[i][j] = 0;
			if (j == i) (*inverse)[i][j] = 1;
		}
	}

//	matrix_print(stdout, "initial for invert", initial, rows, cols);

	for (i = 0; i < rows; i++) {
		working_value = matrix[i][i];
		next = i + 1;
		while (working_value == 0) {
			if (next == rows - 1 || i == rows - 1) {
				fprintf(stderr, "ERROR: matrix non-invertible\n");
				return;
			}
			temp = matrix[i];
			matrix[i] = matrix[next];
			matrix[next] = temp;
//			printf("swapped matrix[%d] with matrix[%d]: %lf %lf\n", i, next, matrix[i][i], matrix[next][i]);
			temp = (*inverse)[i];
			(*inverse)[i] = (*inverse)[next];
			(*inverse)[next] = temp;
//			printf("swapped inverse[%d] with inverse[%d]: %lf %lf\n", i, next, (*inverse)[i][i], (*inverse)[next][i]);
			working_value = matrix[i][i];
			next++;
		}
//		printf("working on row %d with value: %lf\n", i, working_value);

		for (j = 0; j < cols; j++) {
			matrix[i][j] /= working_value;
			(*inverse)[i][j] /= working_value;
//			printf("divided matrix[%d][%d] and inverse[%d][%d]: %lf %lf\n", i, j, i, j, matrix[i][j], (*inverse)[i][j]);
		}

		for (next = 0; next < rows; next++) {
			working_value = matrix[next][i];
			if (next == i || working_value == 0) continue;
			for (j = 0; j < cols; j++) {
				matrix[next][j] -= matrix[i][j] * working_value;
				(*inverse)[next][j] -= (*inverse)[i][j] * working_value;
			}
		}
//		matrix_print(stdout, "matrix", matrix, rows, cols);
//		matrix_print(stdout, "inverse", (*inverse), rows, cols);
//		printf("\n\n");
	}

	for (i = 0; i < rows; i++) {
		free(matrix[i]);
	}
	free(matrix);
}
