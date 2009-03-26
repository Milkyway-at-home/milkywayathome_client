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
#include <math.h>

#define swap(type, i, j) {type t = i; i = j; j = t;}

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

void copy_matrix(double **m, int r, int c, double **result) {
	int i, j;
	for (i = 0; i < c; i++) {
		for (j = 0; j < r; j++) {
			result[i][j] = m[i][j];
		}
	}
}

void copy_matrix__alloc(double **m, int r, int c, double ***result) {
	init_matrix(result, r, c);
	copy_matrix(m, r, c, *result);
}

void identity_matrix(double **m, int r, int c) {
	int i, j;
	for (i = 0; i < c; i++) {
		for (j = 0; j < r; j++) {
			if (i == j) m[i][j] = 1;
			else m[i][j] = 0;
		}
	}
}

void identity_matrix__alloc(double ***m, int r, int c) {
	init_matrix(m, r, c);
	identity_matrix(*m, r, c);
}

void fwrite_matrix(FILE* file, const char* name, double** m, int r, int c) {
	int i, j;
	fprintf(file, "%s [%d x %d]:\n", name, r, c);
	for (i = 0; i < c; i++) {
		for (j = 0; j < r; j++) {
			fprintf(file, " %.10lf", m[i][j]);
		}
		fprintf(file, "\n");
	}
}

void matrix_transpose__alloc(double** m, int r, int c, double*** result) {
	int i, j;

	(*result) = (double**)malloc(sizeof(double*) * c);

	for (i = 0; i < c; i++) {
		(*result)[i] = (double*)malloc(sizeof(double) * r);
		for (j = 0; j < r; j++) {
			(*result)[i][j] = m[j][i];
		}
	}
}

void matrix_transpose__inline(double** m, int r, int c) {
	int i, j;
	for (i = 0; i < c; i++) {
		for (j = i+1; j < c; j++) {
			swap(double, m[i][j], m[j][i]);
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


int LUP_decomposition(int length, double **A, double **LU, int *P) {
	int i, j, k, k_prime;
	double p, divisor;

	k_prime = 0;
	copy_matrix(A, length, length, LU);

	for (i = 0; i < length; i++) P[i] = i;

	for (k = 0; k < length; k++) {
		p = 0;
		for (i = k; i < length; i++) {
			if (fabs(LU[i][k]) > p) {
				p = fabs(LU[i][k]);
				k_prime = i;
			}
		}

		//This is a singular matrix.
		if (p == 0) return -1;

		swap(int, P[k], P[k_prime]);
		for (i = 0; i < length; i++) {
			swap(double, LU[k][i], LU[k_prime][i]);
		}

		divisor = LU[k][k];
		for (i = k+1; i < length; i++) {
			LU[i][k] = LU[i][k]/divisor;
			for (j = k+1; j < length; j++) {
				LU[i][j] = LU[i][j] - LU[i][k]*LU[k][j];
			}
		}
	}
	return 0;
}

int LUP_decomposition__alloc(int length, double **A, double ***LU, int **p) {
	*p = (int*)malloc(sizeof(int) * length);
	init_matrix(LU, length, length);
	return LUP_decomposition(length, A, *LU, *p);
}

void LUP_solve(int length, double **LU, int *p, double *b, double *result) {
	int i, j;
	double *y;

	y = (double*)malloc(sizeof(double) * length);
	for (i = 0; i < length; i++) {
		y[i] = b[p[i]];
		for (j = 0; j < i; j++) {
			y[i] -= LU[i][j] * y[j];
		}
	}
	for (i = length - 1; i >= 0; i--) {
		result[i] = y[i];
		for (j = i+1; j < length; j++) {
			result[i] -= LU[i][j] * result[j];
		}
		result[i] /= LU[i][i];
	}
	free(y);
}

int matrix_invert(double** initial, int rows, int cols, double** result) {
	double **LU, **identity;
	int i, *p;
	int retval;

	retval= LUP_decomposition__alloc(rows, initial, &LU, &p);
	if (retval != 0) {
		free_matrix(&LU, rows, rows);
		free(p);
		return retval;
	}
	identity_matrix__alloc(&identity, rows, rows);

	for (i = 0; i < cols; i++) {
		LUP_solve(rows, LU, p, identity[i], result[i]);
	}
	matrix_transpose__inline(result, rows, cols);
	free_matrix(&LU, rows, rows);
	free_matrix(&identity, rows, rows);
	free(p);
	return 0;
}

int matrix_invert__alloc(double **initial, int rows, int cols, double ***result) {
	init_matrix(result, rows, rows);
	return matrix_invert(initial, rows, cols, *result);
}

#ifdef MATRIX_TEST
int main(int argc, char **argv) {
	double **matrix, **inverse, **result;
	double **LU, *result2, *b;
	int *p;

	matrix = (double**)malloc(sizeof(double*) * 4);
	matrix[0] = (double*)malloc(sizeof(double) * 4);
	matrix[1] = (double*)malloc(sizeof(double) * 4);
	matrix[2] = (double*)malloc(sizeof(double) * 4);
	matrix[3] = (double*)malloc(sizeof(double) * 4);
	matrix[0][0] = 2;
	matrix[0][1] = 0;
	matrix[0][2] = 2;
	matrix[0][3] = 0.6;
	matrix[1][0] = 3;
	matrix[1][1] = 3;
	matrix[1][2] = 4;
	matrix[1][3] = -2;
	matrix[2][0] = 5;
	matrix[2][1] = 5;
	matrix[2][2] = 4;
	matrix[2][3] = 2;
	matrix[3][0] = -1;
	matrix[3][1] = -2;
	matrix[3][2] = 3.4;
	matrix[3][3] = -1;

	fwrite_matrix(stdout, "matrix", matrix, 4, 4);

	matrix_invert__alloc(matrix, 4, 4, &inverse);
	fwrite_matrix(stdout, "inverse", inverse, 4, 4);

	matrix_multiply(matrix, 4, 4, inverse, 4, 4, &result);
	fwrite_matrix(stdout, "result", result, 4, 4);

	printf("\n\n");

	matrix[0][0] = 1;
	matrix[0][1] = 3;
	matrix[0][2] = 3;
	matrix[1][0] = 1;
	matrix[1][1] = 4;
	matrix[1][2] = 3;
	matrix[2][0] = 1;
	matrix[2][1] = 3;
	matrix[2][2] = 4;
	/*
	matrix[0][0] = 1;
	matrix[0][1] = 2;
	matrix[0][2] = 0;
	matrix[1][0] = 3;
	matrix[1][1] = 4;
	matrix[1][2] = 4;
	matrix[2][0] = 5;
	matrix[2][1] = 6;
	matrix[2][2] = 3;
	*/

	b = (double*)malloc(sizeof(double) * 3);
	result2 = (double*)malloc(sizeof(double) * 3);
	b[0] = 3;
	b[1] = 7;
	b[2] = 8;

	fwrite_matrix(stdout, "matrix", matrix, 3, 3);
	matrix_invert__alloc(matrix, 3, 3, &inverse);

	LUP_decomposition__alloc(3, matrix, &LU, &p);

	fwrite_matrix(stdout, "LU", LU, 3, 3);
	printf("p: %d %d %d\n", p[0], p[1], p[2]);

	LUP_solve(3, LU, p, b, result2);

	printf("result: %lf %lf %lf\n", result2[0], result2[1], result2[2]);
	
	fwrite_matrix(stdout, "inverse", inverse, 3, 3);
	matrix_multiply(matrix, 3, 3, inverse, 3, 3, &result);
	fwrite_matrix(stdout, "result", result, 3, 3);
}
#endif

