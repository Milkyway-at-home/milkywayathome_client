#ifndef FGDO_ASYNCHRONOUS_PSO_H
#define FGDO_ASYNCHRONOUS_PSO_H

#include <stdio.h>

#include "asynchronous_search.h"
#include "bounds.h"
#include "population.h"
#include "redundancy.h"

#define NEWTON_ERROR_RANGE 1
#define NEWTON_UPDATE_RANGE 2
#define NEWTON_LINE_SEARCH 3

#define NMS_HESSIAN 1
#define NMS_LINE_SEARCH 2

typedef struct particle_swarm_optimization {
	int current_particle, size;
	int number_parameters;
	int remove_outliers;
	double w, c1, c2;
	long analyzed;

	BOUNDS *bounds;

	double global_best_fitness, *global_best;
	POPULATION *particles;
	POPULATION *velocities;
	POPULATION *local_best;

	REDUNDANCY *current_redundancy;
	REDUNDANCY **redundancies;
} PARTICLE_SWARM_OPTIMIZATION;

ASYNCHRONOUS_SEARCH* get_asynchronous_particle_swarm();

int create_particle_swarm(char* search_name, int number_arguments, char** arguments, int number_parameters, double *point, double *range, BOUNDS *bounds);
int read_particle_swarm(char* search_name, void** search_data);
int checkpoint_particle_swarm(char* search_name, void* search_data);
int pso_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);
int pso_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);
#endif
