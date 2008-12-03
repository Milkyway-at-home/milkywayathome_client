#ifndef GEM_NEWTON_METHOD_H
#define GEM_NEWTON_METHOD_H

void newton_method(int number_parameters, double *point, double *step, int iterations);
void randomized_newton_method(int number_parameters, double *point, double *step, int evaluations_per_iteration, int iterations);
#endif
