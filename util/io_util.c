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

/****
	*	Functions for printing parameters to files/output
*****/
void fwrite_double_array(FILE *file, const char *array_name, int size, double *array_t) {
	int i;

	fprintf(file, "%s[%d]: ", array_name, size);
	for (i = 0; i < size; i++) {
		fprintf(file, "%.20lf", array_t[i]);
		if (i < size-1) fprintf(file, ", ");
	}
	fprintf(file, "\n");
}

void fwrite_int_array(FILE *file, const char *array_name, int size, int *array_t) {
	int i;

	fprintf(file, "%s[%d]: ", array_name, size);
	for (i = 0; i < size; i++) {
		if (i == 0) fprintf(file, " %d", array_t[i]);
		else fprintf(file, ", %d", array_t[i]);
	}
	fprintf(file, "\n");
}

/****
	*	Functions for reading parameters from files
*****/
int fread_double_array(FILE *file, const char *array_name, double** array_t) {
	int i, size;
	fscanf(file, array_name);
	fscanf(file, "[%d]: ", &size);

	(*array_t) = (double*)malloc(sizeof(double) * size);

	for (i = 0; i < size; i++) {
		if (fscanf(file, "%lf", &(*array_t)[i]) != 1) {
			fprintf(stderr,"Error reading into %s\n",array_name);
			exit(-1);
		}
		if (i < size-1) fscanf(file, ", ");
	}
	fscanf(file, "\n");
	return size;
}

int fread_double_array__realloc(FILE *file, const char *array_name, int *size, double **array_t) {
	int i, read_size;
	int c;
	fscanf(file, array_name);
	c = fgetc(file);
	if (c == '[') {
		fscanf(file, "%d]: ", &read_size);
	} else {
		fscanf(file, "[%d]: ", &read_size);
	}

	if (read_size <= 0) {
		fprintf(stderr, "Error reading %s, invalid size: %d\n", array_name, read_size);
		return -1;
	}	
	if (read_size != *size) {
		*size = read_size;
		(*array_t) = (double*)realloc(*array_t, read_size * sizeof(double));
	}

	for (i = 0; i < read_size; i++) {
		if (fscanf(file, "%lf", &((*array_t)[i])) != 1) {
			fprintf(stderr, "Error reading into %s, invalid data\n", array_name);
			return -1;
		}
		if (i < read_size-1) fscanf(file, ", ");
	}
	fscanf(file, "\n");
	return read_size;
}

int fread_double_array__no_alloc(FILE *file, const char *array_name, int size, double *array_t) {
	int i, read_size;
	fscanf(file, array_name);
	fscanf(file, "[%d]: ", &read_size);

	if (read_size != size) {
		fprintf(stderr, "Error reading into %s, invalid length: %d, expected: %d\n", array_name, read_size, size);
		return -1;
	}

	for (i = 0; i < size; i++) {
		if (fscanf(file, "%lf", &(array_t[i])) != 1) {
			fprintf(stderr, "Error reading into %s, invalid data\n", array_name);
			return -1;
		}
		if (i < size-1) fscanf(file, ", ");
	}
	fscanf(file, "\n");
	return size;
}

int fread_int_array(FILE *file, const char *array_name, int **array_t) {
	int i, size;
	fscanf(file, array_name);
	fscanf(file, "[%d]: ", &size);

	(*array_t) = (int*)malloc(sizeof(int) * size);

	for (i = 0; i < size; i++) {
		if (fscanf(file, "%d", &(*array_t)[i]) != 1) {
			fprintf(stderr,"Error reading into %s\n",array_name);
			exit(-1);
		}
		if (i < size-1) fscanf(file, ", ");
	}
	fscanf(file, "\n");
	return size;
}

int fread_int_array__no_alloc(FILE *file, const char *array_name, int size, int *array_t) {
	int i, read_size;
	fscanf(file, array_name);
	fscanf(file, "[%d]: ", &read_size);

	if (read_size != size) {
		fprintf(stderr, "Error reading into %s, invalid length: %d, expected: %d\n", array_name, read_size, size);
		return -1;
	}

	for (i = 0; i < size; i++) {
		if (fscanf(file, "%d", &(array_t[i])) != 1) {
			fprintf(stderr, "Error reading into %s, invalid data\n", array_name);
			return -1;
		}
		if (i < size-1) fscanf(file, ", ");
	}
	fscanf(file, "\n");
	return size;
}
