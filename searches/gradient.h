#ifndef GEM_GRADIENT_H
#define GEM_GRADIENT_H

#include<stdio.h>

typedef struct gradient {
	int iteration;
	int set_values;
	int number_parameters;
	double* step;
	double* point;

	int** set_evaluations;
	double** evaluations;
	double* values;
} GRADIENT;


void get_gradient(int number_parameters, double *point, double *step, double *gradient);
int gradient_below_threshold(int number_parameters, double* gradient, double threshold);


#endif
