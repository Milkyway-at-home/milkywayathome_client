#ifndef GEO_PARTICLE_SWARM_H
#define GEO_PARTICLE_SWARM_H

#include "population.h"

typedef struct particle_swarm {
	double velocity_weight;
	double local_weight;
	double global_weight;

	int next_generated;
	int local_best;

	double **velocity;
	double **local_best_parameters;
	double *local_best_fitness;
	double *global_best_parameters;
	double global_best_fitness;
} PARTICLE_SWARM;


void start_particle_swarm(char search_path[512], char search_parameters[512], double* min_parameters, double* max_parameters, int number_parameters, POPULATION **population);

void particle_swarm__insert_individual(POPULATION *population, double* parameters, double fitness, char *metadata);

void particle_swarm__get_individual(POPULATION *population, double** parameters, char **metadata);

#endif
