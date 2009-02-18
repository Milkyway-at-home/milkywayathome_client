#ifndef FGDO_BOUNDS
#define FGDO_BOUNDS

#include "stdio.h"

typedef struct bounds {
	int number_parameters;
	double *min_bound, *max_bound;
	int *in_radians;
} BOUNDS;

void new_bounds(BOUNDS **bounds, int number_parameters, double *min_bound, double *max_bound, int *in_radians);
void free_bounds(BOUNDS **bounds);

void bound_parameters(double *parameters, BOUNDS *bounds);
void bound_velocity(double *parameters, double *velocity, BOUNDS *bounds);

void fwrite_bounds(FILE *file, BOUNDS *bounds);
void fread_bounds(FILE *file, BOUNDS **bounds);

#endif

