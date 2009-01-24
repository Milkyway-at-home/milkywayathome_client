#ifndef GEM_NEWTON_METHOD_H
#define GEM_NEWTON_METHOD_H

void newton_method(int number_parameters, double *point, double *step, int iterations);
void randomized_newton_method(int number_parameters, double *point, double *step, int type, int evaluations_per_iteration, int iterations, int ls_evaluations, int ls_iterations);

#endif
