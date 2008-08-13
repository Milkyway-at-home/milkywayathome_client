#ifndef GEO_POPULATION_H
#define GEO_POPULATION_H

#include <stdio.h>

#define metadata_size 2048

typedef struct population POPULATION;

typedef void (*get_function_type)(POPULATION*, double**, char**);

typedef void (*insert_function_type)(POPULATION*, double*, double, char*);

struct population {
	char* search_path;
	char* search_parameters;

	int current_size;
	int max_size;

	int current_evaluation;
	int max_evaluations;

	int number_parameters;
	double** individuals;
	double* fitness;

	double* min_parameters;
	double* max_parameters;

	get_function_type get_individual;
	insert_function_type insert_individual;

	void* parameters;
};

POPULATION* fread_population(FILE* file);
POPULATION* read_population(char path[512]);

void fwrite_population(FILE* file, POPULATION* population);
void write_population(char path[512], POPULATION* population);

void fwrite_population_statistics(FILE *file, POPULATION* population);

POPULATION* new_population(char search_path[512], char search_parameters[512], double *min_parameters, double *max_parameters, int number_parameters, int population_size, int max_evaluations);

void insert_sorted(POPULATION* population, double* parameters, double fitness);
void replace(POPULATION* population, int position, double* parameters, double fitness);
void replace_if_better(POPULATION* population, int position, double* parameters, double fitness);

void get_n_distinct(POPULATION *population, int number_parents, double ***parents, double **parent_fitness);
void bound_parameters(POPULATION *population, double *parameters);

#endif
