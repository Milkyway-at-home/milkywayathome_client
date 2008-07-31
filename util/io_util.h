#ifndef GEM_IO_UTIL_H
#define GEM_IO_UTIL_H

#include <stdio.h>
#include <stdlib.h>

void print_double_array(FILE *file, const char *array_name, int size, double *array_t);
void print_int_array(FILE *file, const char *array_name, int size, int *array_t);

void read_double_array(FILE *file, const char *array_name, int size, double** array_t);
void read_int_array(FILE *file, const char *array_name, int size, int **array_t);

#endif
