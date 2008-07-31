#include <stdio.h>
#include <stdlib.h>

#include "io_util.h"

/****
	*	Functions for printing parameters to files/output
*****/
void print_double_array(FILE *file, const char *array_name, int size, double *array_t) {
	int i;

	fprintf(file, array_name);
	for (i = 0; i < size; i++) {
		fprintf(file, "%lf", array_t[i]);
		if (i < size-1) fprintf(file, ", ");
	}
	fprintf(file, "\n");
}

void print_int_array(FILE *file, const char *array_name, int size, int *array_t) {
	int i;

	fprintf(file, array_name);
	for (i = 0; i < size; i++) {
		if (i == 0) fprintf(file, " %d", array_t[i]);
		else fprintf(file, ", %d", array_t[i]);
	}
	fprintf(file, "\n");
}

/****
	*	Functions for reading parameters from files
*****/
void read_double_array(FILE *file, const char *array_name, int size, double** array_t) {
	int i;
	(*array_t) = (double*)malloc(sizeof(double) * size);

	fscanf(file, array_name);
	for (i = 0; i < size; i++) {
		if (fscanf(file, "%lf", &(*array_t)[i]) != 1) {
			fprintf(stderr,"Error reading into %s\n",array_name);
			exit(-1);
		}
		if (i < size-1) fscanf(file, ", ");
	}
	fscanf(file, "\n");
}

void read_int_array(FILE *file, const char *array_name, int size, int **array_t) {
	int i;
	(*array_t) = (int*)malloc(sizeof(int) * size);

	fscanf(file, array_name);
	for (i = 0; i < size; i++) {
		if (fscanf(file, "%d", &(*array_t)[i]) != 1) {
			fprintf(stderr,"Error reading into %s\n",array_name);
			exit(-1);
		}
		if (i < size-1) fscanf(file, ", ");
	}
	fscanf(file, "\n");
}
