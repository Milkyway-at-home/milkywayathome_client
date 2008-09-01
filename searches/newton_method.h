#ifndef GEM_NEWTON_METHOD_H
#define GEM_NEWTON_METHOD_H

void synchronous_newton_method(char* search_path, char* search_parameters, double* point, double* step, int number_parameters);
void randomized_newton_method(char* search_path, char* search_parameters, double* point, double* step, int number_parameters);

#endif
