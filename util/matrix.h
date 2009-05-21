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

#ifndef FGDO_MATRIX_H
#define FGDO_MATRIX_H

#include <stdio.h>

extern char MATRIX__ERROR_MSG[1024];

#define MATRIX__MxM_INVALID_MATCH 1
#define MATRIX__READ_ERROR 2

void new_matrix(double ***m, int r, int c);
void free_matrix(double ***m, int r, int c);

void copy_matrix(double **m, int r, int c, double **result);
void copy_matrix__alloc(double **m, int r, int c, double ***result);

int fwrite_matrix(FILE* file, const char* name, double **m, int r, int c);
int fread_matrix(FILE* file, const char* name, double **m, int r, int c);
int fread_matrix__alloc(FILE* file, char** name, double ***m, int *r, int *c);

void matrix_transpose__alloc(double** m, int r, int c, double*** result);
void matrix_transpose__inline(double** m, int r, int c);

int matrix_multiply(double** m1, int r1, int c1, double** m2, int r2, int c2, double** result);
int matrix_multiply__alloc(double** m1, int r1, int c1, double** m2, int r2, int c2, double*** result);

int matrix_vector_multiply(double **m, int r1, int c1, double *v, int r2, double *result);
int matrix_vector_multiply__alloc(double **m, int r1, int c1, double *v, int r2, double **result);

int matrix_invert(double** initial, int rows, int cols, double** inverse);
int matrix_invert__alloc(double** initial, int rows, int cols, double*** inverse);
int matrix_invert__inline(double** initial, int rows, int cols);

#endif
