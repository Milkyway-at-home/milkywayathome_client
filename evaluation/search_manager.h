#ifndef FGDO_SEARCH_MANAGER_H
#define FGDO_SEARCH_MANAGER_H

#include "../searches/search_parameters.h"

void initialize_search_manager();

int insert_workunit(char* search_name, double fitness, double* parameters, char* metadata);
int generate_workunits(int number_workunits, SEARCH_PARAMETERS ***parameters);

#endif
