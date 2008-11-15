#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

#include "population.h"


int new_population(int max_size, int number_parameters, POPULATION** population) {
	(*population) = (POPULATION*)malloc(sizeof(POPULATION));
	(*population)->size = 0;
	(*population)->max_size = max_size;
	(*population)->number_parameters = number_parameters;
	(*population)->individuals = (double**)malloc(sizeof(double*) * max_size);
	(*population)->fitness = (double*)malloc(sizeof(double) * max_size);

	return 1;
}

void free_population(POPULATION* population) {
	int i;
	for (i = 0; i < population->max_size; i++) {
		if (population->individuals[i] != NULL) free(population->individuals[i]);
	}
	free(population->fitness);
}

void insert(POPULATION* population, double* parameters, double fitness) {
	if (population->individuals[population->size] == NULL) {
		population->individuals[population->size] = (double*)malloc(sizeof(double) * population->number_parameters);
	}
	memcpy(population->individuals[population->size], parameters, sizeof(double) * population->number_parameters);
	population->fitness[population->size] = fitness;
	population->size++;
}

void insert_sorted(POPULATION* population, double* parameters, double fitness) {
	int i, j;
        for (i = 0; i < population->size; i++) {
                if (fitness > population->fitness[i]) {
			if (population->size == population->max_size) free(population->individuals[population->size-1]);
			for (j = population->size - 1; j > i; j--) {
				population->individuals[j] = population->individuals[j-1];
				population->fitness[j] = population->fitness[j-1];
			}
			population->individuals[i] = (double*)malloc(sizeof(double) * population->number_parameters);
			memcpy(population->individuals[i], parameters, sizeof(double) * population->number_parameters);
			population->fitness[i] = fitness;
			break;
		}
        }
}

void replace(POPULATION* population, int position, double* parameters, double fitness) {
	if (population->individuals[position] == NULL) {
		population->individuals[position] = (double*)malloc(sizeof(double) * population->number_parameters);
		population->size++;
	}
	memcpy(population->individuals[position], parameters, sizeof(double) * population->number_parameters);
	population->fitness[position] = fitness;
}

void replace_if_better(POPULATION* population, int position, double* parameters, double fitness) {
	if (population->individuals[position] != NULL && population->fitness[position] <= fitness) return;
	if (population->individuals[position] == NULL) {
		population->individuals[position] = (double*)malloc(sizeof(double) * population->number_parameters);
		population->size++;
	}
	memcpy(population->individuals[position], parameters, sizeof(double) * population->number_parameters);
	population->fitness[position] = fitness;
}

void get_n_distinct(POPULATION *population, int number_parents, POPULATION **n_distinct) {
	int i, j, k, target;
	int *parent_positions;

	new_population(number_parents, population->number_parameters, n_distinct);

	parent_positions = (int*)malloc(sizeof(int) * number_parents);

	for (i = 0; i < number_parents; i++) {
		target = (int)(drand48() * (population->max_size - i));
		for (j = 0; j < i; j++) {
			for (k = 0; k < i; k++) {
				if (target == parent_positions[k]) target++;
			}
		}
		parent_positions[i] = target;
	}

	for (i = 0; i < number_parents; i++) {
		insert_sorted((*n_distinct), population->individuals[parent_positions[i]], population->fitness[parent_positions[i]]);
	}
}

int fread_population(FILE* file, POPULATION **population) {
	int i, j, size, max_size, number_parameters, position;
	double fitness;

	fscanf(file, "size: %d, max_size: %d\n", &size, &max_size);
	fscanf(file, "number_parameters: %d\n", &number_parameters);
	new_population(max_size, number_parameters, population);

	for (i = 0; i < size; i++) {
		fscanf(file, "[%d] %lf :", &position, &fitness);
		(*population)->fitness[position] = fitness;
		(*population)->individuals[position] = (double*)malloc(sizeof(double) * number_parameters);
		for (j = 0; j < number_parameters; j++) {
			fscanf(file, " %lf", &((*population)->individuals[position][j]));
		}
		fscanf(file, "\n");
	}
	return 1;
}

int read_population(char *filename, POPULATION **population) {
	int result;
	FILE *file;

	file = fopen(filename, "r");
	if (file == NULL) return -1;
	result = fread_population(file, population);
	fclose(file);
	return result;
}


int fwrite_population(FILE *file, POPULATION *population) {
	int i, j;

	fprintf(file, "size: %d, max_size: %d\n", population->size, population->max_size);
	fprintf(file, "number_parameters: %d\n", population->number_parameters);
	for (i = 0; i < population->max_size; i++) {
		if (population->individuals[i] == NULL) continue;

		fprintf(file, "[%d] %lf :", i, population->fitness[i]);
		for (j = 0; j < population->number_parameters; j++) {
			fprintf(file, " %lf", population->individuals[i][j]);
		}
		fprintf(file, "\n");
	}
	fflush(file);
	return 1;
}

int write_population(char *filename, POPULATION *population) {
	int result;
	FILE *file;

	file = fopen(filename, "w+");
	if (file == NULL) return -1;
	result = fwrite_population(file, population);
	fclose(file);
	return result;
}

void fwrite_population_statistics(FILE* file, POPULATION *population) {
	double min_fitness, avg_fitness, max_fitness;
	int i;

	min_fitness = population->fitness[0];
	max_fitness = population->fitness[0];
	avg_fitness = 0;
	for (i = 0; i < population->size; i++) {
		avg_fitness += population->fitness[i];
		if (min_fitness > population->fitness[i]) min_fitness = population->fitness[i];
		else if (max_fitness < population->fitness[i]) max_fitness = population->fitness[i];
	}
	avg_fitness /= population->size;

	fprintf(file, "min: %lf, avg: %lf, max: %lf", min_fitness, avg_fitness, max_fitness);
}
