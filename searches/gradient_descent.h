#ifndef FGDO_GRADIENT_DESCENT_H
#define FGDO_GRADIENT_DESCENT_H

void synchronous_gradient_descent(char* search_path, char* search_parameters, double* point, double* step, int number_parameters);
void synchronous_conjugate_gradient_descent(char* search_path, char* search_parameters, double* point, double* step, int number_parameters);

#endif
