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

#ifndef GEM_IO_UTIL_H
#define GEM_IO_UTIL_H

#include <stdio.h>
#include <stdlib.h>

void fwrite_double_array(FILE *file, const char *array_name, int size, double *array_t);
void fwrite_int_array(FILE *file, const char *array_name, int size, int *array_t);

int fread_double_array(FILE *file, const char *array_name, double** array_t);
int fread_double_array__no_alloc(FILE *file, const char *array_name, int number_parameters, double* array_t);
int fread_int_array(FILE *file, const char *array_name, int **array_t);
int fread_int_array__no_alloc(FILE *file, const char *array_name, int number_parameters, int *array_t);

#endif
