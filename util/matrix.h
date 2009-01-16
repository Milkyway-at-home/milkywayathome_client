#ifndef GEM_MATRIX_H
#define GEM_MATRIX_H

#include <stdio.h>

void init_matrix(double ***m, int r, int c);
void free_matrix(double ***m, int r, int c);

void matrix_print(FILE* file, const char* name, double**m, int r, int c);
void matrix_transpose(double** m, int r, int c, double*** result);
void matrix_multiply(double** m1, int r1, int c1, double** m2, int r2, int c2, double*** result);
void matrix_vector_multiply(double** m2, int r2, int c2, double* m1, int c1, double** result);
void vector_matrix_multiply(double* m1, int c1, double** m2, int c2, int r2, double** result);
void matrix_invert(double** initial, int rows, int cols, double*** inverse);

#endif
