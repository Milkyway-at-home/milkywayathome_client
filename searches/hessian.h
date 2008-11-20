#ifndef GEM_HESSIAN_H
#define GEM_HESSIAN_H

#include "stdio.h"

typedef struct hessian {
	int iteration;
	int set_values;
	int number_parameters;
	double* step;
	double* point;

	int*** set_evaluations;
	double*** evaluations;
	double** values;
} HESSIAN;


void create_hessian(double* point, double* step, int iteration, int number_parameters, HESSIAN** hessian);
void free_hessian(HESSIAN* hessian);

int hessian__complete(HESSIAN* hessian);

void hessian__insert_individuals(HESSIAN* hessian, int number_individuals, double* fitness, char** metadata);
void hessian__get_individuals(HESSIAN* hessian, int number_individuals, double*** parameters, char*** metadata);

void randomized_hessian(double** actual_points, double* center, double* fitness, int number_points, int number_parameters, double*** hessian, double** gradient);

void fprintf_hessian(FILE* file, HESSIAN* hessian);

#endif
