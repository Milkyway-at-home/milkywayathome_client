#ifndef GEM_HESSIAN_H
#define GEM_HESSIAN_H

#include "stdio.h"

void get_hessian(int number_parameters, double *point, double *step, double **hessian);
void randomized_hessian(double** actual_points, double* center, double* fitness, int number_points, int number_parameters, double*** hessian, double** gradient);

#endif
