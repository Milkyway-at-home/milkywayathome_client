#ifndef GEM_NEWTON_METHOD_H
#define GEM_NEWTON_METHOD_H

#include "../evaluation/search_manager.h"

void synchronous_newton_method(char* search_path, char* search_parameters, double* point, double* step, int number_parameters);
void randomized_newton_method(char* search_path, char* search_parameters, double* point, double* step, int number_parameters);

int init_newton_method(char* search_name, SEARCH* search);

#endif
