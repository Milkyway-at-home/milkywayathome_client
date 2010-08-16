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

#include "io_util.h"
#include "milkyway_util.h"

/****
	*	Functions for printing parameters to files/output
*****/
void fwrite_double_array(FILE *file, const char *array_name, double* array_t, size_t size)
{
	size_t i;

	fprintf(file, "%s[%u]: ", array_name, size);
	for (i = 0; i < size; i++)
    {
		fprintf(file, "%.20lf", array_t[i]);
		if (i < size - 1)
            fprintf(file, ", ");
	}
	fprintf(file, "\n");
}

void fwrite_int_array(FILE* file, const char* array_name, int* array_t, size_t size)
{
	size_t i;
	fprintf(file, "%s[%u]: ", array_name, size);
	for (i = 0; i < size; i++)
    {
		if (i == 0)
            fprintf(file, " %d", array_t[i]);
		else
            fprintf(file, ", %d", array_t[i]);
	}
	fprintf(file, "\n");
}

/****
	*	Functions for reading parameters from files
*****/
double* fread_double_array(FILE* file, const char* array_name, unsigned int* sizeOut)
{
    unsigned int i, size;
    double* arr;

	fscanf(file, array_name);
	fscanf(file, "[%u]: ", &size);

    arr = mallocSafe(sizeof(double) * size);

	for (i = 0; i < size; i++)
    {
		if (fscanf(file, "%lf", &arr[i]) != 1)
        {
			fprintf(stderr, "Error reading into %s\n", array_name);
			exit(-1);
		}
		if (i < size - 1)
            fscanf(file, ", ");
	}
	fscanf(file, "\n");

    if (sizeOut)
        *sizeOut = size;

	return arr;
}

int* fread_int_array(FILE *file, const char *array_name, unsigned int* sizeOut)
{
	unsigned int i, size;
	fscanf(file, array_name);
	fscanf(file, "[%u]: ", &size);
    int* arr;

	arr = mallocSafe(sizeof(int) * size);

	for (i = 0; i < size; i++)
    {
		if (fscanf(file, "%d", &arr[i]) != 1)
        {
			fprintf(stderr,"Error reading into %s\n", array_name);
			exit(-1);
		}

		if (i < size - 1)
            fscanf(file, ", ");
	}
	fscanf(file, "\n");

    if (sizeOut)
        *sizeOut = size;

	return arr;
}

