#ifndef FGDO_GRADIENT_DESCENT_H
#define FGDO_GRADIENT_DESCENT_H

void synchronous_gradient_descent(int number_arguments, char **arguments, int number_parameters, double *point, double *step);
void synchronous_conjugate_gradient_descent(int number_arguments, char **arguments, int number_parameters, double *point, double *step);

#endif
