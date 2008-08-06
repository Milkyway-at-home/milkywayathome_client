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


void create_gradient(double* point, double* step, int iteration, int number_parameters, GRADIENT** gradi);
void free_gradient(GRADIENT* gradient);

int gradient__complete(GRADIENT* gradient);

void gradient__insert_individuals(GRADIENT* gradient, int number_individuals, double* fitness, char** metadata);
void gradient__get_individuals(GRADIENT* gradient, int number_individuals, double*** parameters, char*** metadata);

void synchronous_get_gradient(double* point, double* step, int number_parameters, GRADIENT** gradient);

void fprintf_gradient(FILE *file, GRADIENT* gradient);

#endif
