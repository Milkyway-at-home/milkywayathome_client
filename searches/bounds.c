#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "bounds.h"
#include "../util/io_util.h"

void new_bounds(BOUNDS **bounds, int number_parameters, double *min_bound, double *max_bound, int *in_radians) {
	int i;
	(*bounds) = (BOUNDS*)malloc(sizeof(BOUNDS));
	BOUNDS *b = (*bounds);

	b->number_parameters = number_parameters;
	b->min_bound = (double*)malloc(sizeof(double) * number_parameters);
	b->max_bound = (double*)malloc(sizeof(double) * number_parameters);
	b->in_radians = (int*)malloc(sizeof(int) * number_parameters);

	memcpy(b->min_bound, min_bound, sizeof(double) * number_parameters);
	memcpy(b->max_bound, max_bound, sizeof(double) * number_parameters);
	memcpy(b->in_radians, in_radians, sizeof(int) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		if (b->in_radians[i]) {
			b->min_bound[i] = -M_PI;
			b->max_bound[i] = M_PI;
		}
	}
}

void free_bounds(BOUNDS **bounds) {
	BOUNDS *b = (*bounds);
	free(b->min_bound);
	free(b->max_bound);
	free(b->in_radians);
	free((*bounds));
	(*bounds) = NULL;
}

void bound_parameters(double* parameters, BOUNDS *b) {
	int i;
	if (b->in_radians == NULL) {
		for (i = 0; i < b->number_parameters; i++) {
			if (parameters[i] > b->max_bound[i]) parameters[i] = b->max_bound[i];
			if (parameters[i] < b->min_bound[i]) parameters[i] = b->min_bound[i];
		}       
	} else {
		for (i = 0; i < b->number_parameters; i++) {
			if (b->in_radians[i]) {
				while (parameters[i] > M_PI || parameters[i] < -M_PI) {
					if (parameters[i] > M_PI) parameters[i] -= (M_PI + M_PI);
					if (parameters[i] < -M_PI) parameters[i] += (M_PI + M_PI);
				}
			} else {
				if (parameters[i] > b->max_bound[i]) parameters[i] = b->max_bound[i];
				if (parameters[i] < b->min_bound[i]) parameters[i] = b->min_bound[i];
			}
		}       
	}
}

void bound_velocity(double *parameters, double *velocity, BOUNDS *b) {
	int i;
	if (b->in_radians == NULL) {
		for (i = 0; i < b->number_parameters; i++) {
			if (parameters[i] + velocity[i] > b->max_bound[i]) velocity[i] = b->max_bound[i] - parameters[i];
			if (parameters[i] + velocity[i] < b->min_bound[i]) velocity[i] = b->min_bound[i] - parameters[i];
		}       
	} else {
		for (i = 0; i < b->number_parameters; i++) {
			if (b->in_radians[i]) {
				if (velocity[i] > M_PI) velocity[i] = M_PI;
				if (velocity[i] < -M_PI) velocity[i] = -M_PI;
			} else {
				if (parameters[i] + velocity[i] > b->max_bound[i]) velocity[i] = b->max_bound[i] - parameters[i];
				if (parameters[i] + velocity[i] < b->min_bound[i]) velocity[i] = b->min_bound[i] - parameters[i];
			}
		}       
	}
}

void bound_step(double *parameters, double *direction, double step, BOUNDS *b) {
}

void fwrite_bounds(FILE *file, BOUNDS *b) {
	print_double_array(file, "min_bound", b->number_parameters, b->min_bound);
	print_double_array(file, "max_bound", b->number_parameters, b->max_bound);
	print_int_array(file, "in_radians", b->number_parameters, b->in_radians);
}

void fread_bounds(FILE *file, BOUNDS **bounds) {
	(*bounds) = (BOUNDS*)malloc(sizeof(BOUNDS));
	BOUNDS *b = (*bounds);

	b->number_parameters = read_double_array(file, "min_bound", &(b->min_bound));
	read_double_array(file, "max_bound", &(b->max_bound));
	read_int_array(file, "in_radians", &(b->in_radians));
}

