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

#ifndef GEM_MATRIX_H
#define GEM_MATRIX_H

#include <stdio.h>

void init_matrix(double ***m, int r, int c);
void free_matrix(double ***m, int r, int c);

void copy_matrix(double **m, int r, int c, double **result);
void copy_matrix__alloc(double **m, int r, int c, double ***result);

void fwrite_matrix(FILE* file, const char* name, double**m, int r, int c);

void matrix_transpose__alloc(double** m, int r, int c, double*** result);
void matrix_transpose__inline(double** m, int r, int c);

void matrix_multiply(double** m1, int r1, int c1, double** m2, int r2, int c2, double*** result);
void matrix_vector_multiply(double** m2, int r2, int c2, double* m1, int c1, double** result);
void vector_matrix_multiply(double* m1, int c1, double** m2, int c2, int r2, double** result);

void matrix_invert(double** initial, int rows, int cols, double** inverse);
void matrix_invert__alloc(double** initial, int rows, int cols, double*** inverse);

#endif
