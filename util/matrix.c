#include "stdio.h"
#include "stdlib.h"

void matrix_print(FILE* file, const char* name, double** m, int r, int c) {
	int i, j;
	fprintf(file, "%s [%d x %d]:\n", name, r, c);
	for (i = 0; i < r; i++) {
		for (j = 0; j < c; j++) {
			fprintf(file, " %lf", m[i][j]);
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

void matrix_vector_multiply(double* m1, int c1, double** m2, int r2, int c2, double** result) {
	int j, k;

	(*result) = (double*)malloc(sizeof(double) * c2);
	for (j = 0; j < c2; j++) {
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

	printf("\n\ndid setup\n");
	matrix_print(stdout, "initial", initial, rows, cols);

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
		printf("working on row %d with value: %lf\n", i, working_value);

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
		matrix_print(stdout, "matrix", matrix, rows, cols);
		matrix_print(stdout, "inverse", (*inverse), rows, cols);
		printf("\n\n");
	}

	for (i = 0; i < rows; i++) {
		free(matrix[i]);
	}
	free(matrix);
}
