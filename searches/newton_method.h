#ifndef FGDO_NEWTON_METHOD_H
#define FGDO_NEWTON_METHOD_H

void newton_step_i(int number_parameters, double **hessian, double *gradient, double **step);
void newton_step(int number_parameters, double **hessian, double *gradient, double *step);

#endif
