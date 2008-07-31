#ifndef GEM_MATRIX_H
#define GEM_MATRIX_H

void matrix_multiply(double** m1, int r1, int c1, double** m2, int r2, int c2, double*** result);
void matrix_vector_multiply(double* m1, int c1, double** m2, int r2, int c2, double** result);
void matrix_invert(double** initial, int rows, int cols, double*** inverse);

#endif
