#ifndef GEM_GENETIC_SEARCH_H
#define GEM_GENETIC_SEARCH_H

#include "population.h"

typedef struct genetic_search {
	int recombination_type;

	double **parents;
	double *parent_fitness;
	int current_generated;

	double mutation_rate;
	int number_parents, number_children;
	double simplex_l1, simplex_l2;
	double crossover_rate, crossover_scale;
} GENETIC_SEARCH;

void parse_genetic_search(char search_parameters[512], GENETIC_SEARCH *gs, int *population_size, int *max_evaluations);

void genetic_search__insert_individual(POPULATION *population, double* parameters, double fitness, char *metadata);

void genetic_search__get_individual(POPULATION *population, double** parameters, char **metadata);

void start_genetic_search(char *search_path, char *search_parameters, double *min_parameters, double *max_parameters, int number_parameters, POPULATION **population);

#endif
