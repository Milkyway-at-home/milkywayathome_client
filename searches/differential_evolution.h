#ifndef GEO_DIFFERENTIAL_EVOLUTION_H
#define GEO_DIFFERENTIAL_EVOLUTION_H

#include "population.h"

typedef struct differential_evolution {
	int current_best;

	int parent_type;
	double parent_scale;
	int number_pairs, pair_type;
	double pair_scale;
	int recombination_type;
	double crossover_rate;

	int next_generated;
} DIFFERENTIAL_EVOLUTION;


void start_differential_evolution(char search_path[512], char search_parameters[512], double* min_parameters, double* max_parameters, int number_parameters, POPULATION **population);

void differential_evolution__insert_individual(POPULATION *population, double* parameters, double fitness, char *metadata);

void differential_evolution__get_individual(POPULATION *population, double** parameters, char **metadata);

#endif
