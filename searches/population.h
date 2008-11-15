#ifndef FGDO_POPULATION_H
#define FGDO_POPULATION_H

#include <stdio.h>
#include "../util/settings.h"

typedef struct population {
	int size;
	int max_size;

	int number_parameters;
	double** individuals;
	double* fitness;
} POPULATION;

int fread_population(FILE* file, POPULATION** p);
int read_population(char path[512], POPULATION** p);

int fwrite_population(FILE* file, POPULATION* population);
int write_population(char path[512], POPULATION* population);

int new_population(int max_size, int number_parameters, POPULATION** population);
void free_population(POPULATION* population);

void add(POPULATION* population, double* parameters, double fitness);
void insert_sorted(POPULATION* population, double* parameters, double fitness);
void replace(POPULATION* population, int position, double* parameters, double fitness);
void replace_if_better(POPULATION* population, int position, double* parameters, double fitness);

void get_n_distinct(POPULATION *population, int number_parents, POPULATION *n_distinct);

#endif
