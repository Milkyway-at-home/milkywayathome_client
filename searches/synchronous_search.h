#ifndef GEO_SYNCHRONOUS_SEARCH_H
#define GEO_SYNCHRONOUS_SEARCH_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../searches/population.h"

void synchronous_search(char *search_path, char *search_parameters, double *min_parameters, double *max_parameters, int number_parameters, 
				void (*init_search)(char*, char*, double*, double*, int, POPULATION**));
void synchronous_parallel_search(char *search_path, char *search_parameters, double *min_parameters, double *max_parameters, int number_parameters, 
				void (*init_search)(char*, char*, double*, double*, int, POPULATION**), int number_workers);

#endif
