#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

#include "population.h"
#include "regression.h"

#include "../util/matrix.h"

int new_population(int max_size, int number_parameters, POPULATION** population) {
	int i;
	(*population) = (POPULATION*)malloc(sizeof(POPULATION));
	(*population)->size = 0;
	(*population)->max_size = max_size;
	(*population)->number_parameters = number_parameters;
	(*population)->individuals = (double**)malloc(sizeof(double*) * max_size);
	for (i = 0; i < (*population)->max_size; i++) (*population)->individuals[i] = NULL;

	(*population)->os_names = (char**)malloc(sizeof(char*) * max_size);
	for (i = 0; i < (*population)->max_size; i++) (*population)->os_names[i] = NULL;

	(*population)->app_versions = (char**)malloc(sizeof(char*) * max_size);
	for (i = 0; i < (*population)->max_size; i++) (*population)->app_versions[i] = NULL;

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

void get_population_statistics(POPULATION *p, double *best_point, double *best_fitness, double *average_fitness, double *worst_fitness, double *standard_deviation) {
	double best, avg, worst, st_dev;
	int best_pos, i;

	best = p->fitness[0];
	avg = p->fitness[0];
	worst = p->fitness[0];
	best_pos = 0;
	printf("getting stats, size: %d, max_size: %d\n", p->size, p->max_size);
	for (i = 1; i < p->max_size; i++) {
		if (!individual_exists(p, i)) continue;
		avg += p->fitness[i];
		if (p->fitness[i] > best) {
			best = p->fitness[i];
			best_pos = i;
		}
		if (p->fitness[i] < worst) worst = p->fitness[i];
	}
	st_dev = 0;
	for (i = 0; i < p->size; i++) {
		if (!individual_exists(p, i)) continue;
		st_dev += (p->fitness[i] - best) * (p->fitness[i] - best);
	}
	st_dev /= p->size;
	st_dev = sqrt(st_dev);

	(*average_fitness) = avg / p->size;
	(*best_fitness) = best;
	(*worst_fitness) = worst;
	(*standard_deviation) = st_dev;
	for (i = 0; i < p->number_parameters; i++) best_point[i] = p->individuals[best_pos][i];
}

/********
	*	Functions to check for existence of individuals
 ********/
int population_contains(POPULATION *p, double fitness, double *point) {
	int i, j, count;

	for (i = 0; i < p->max_size; i++) {
		if (p->individuals[i] == NULL) continue;
		if (p->fitness[i] == fitness) return 1;

		count = 0;
		for (j = 0; j < p->number_parameters; j++) {
			if (fabs(point[j] - p->individuals[i][j]) > 10e-10) break;
			count++;
		}
		if (count == p->number_parameters) return 1;
	}
	return 0;
}

int individual_exists(POPULATION *p, int position) {
	return p->individuals[position] != NULL;
}

/********
	*	Functions for inserting individuals
 ********/

void insert_individual(POPULATION* population, int position, double* parameters, double fitness) {
	if (population->individuals[position] == NULL) {
		population->individuals[position] = (double*)malloc(sizeof(double) * population->number_parameters);
		population->size++;
	}
	memcpy(population->individuals[position], parameters, sizeof(double) * population->number_parameters);
	population->fitness[position] = fitness;
}

void insert_individual_info(POPULATION* population, int position, double* parameters, double fitness, char *os_name, char* app_version) {
	if (population->individuals[position] == NULL) {
		population->individuals[position] = (double*)malloc(sizeof(double) * population->number_parameters);
		population->size++;
	}
	memcpy(population->individuals[position], parameters, sizeof(double) * population->number_parameters);
	population->fitness[position] = fitness;

	if (population->app_versions[position] == NULL) population->app_versions[position] = (char*)malloc(sizeof(char) * 512);
	if (app_version != NULL) sprintf(population->app_versions[position], "%s", app_version);

	if (population->os_names[position] == NULL) population->os_names[position] = (char*)malloc(sizeof(char) * 512);
	if (os_name != NULL) sprintf(population->os_names[position], "%s", os_name);
}

void insert_incremental(POPULATION* population, double* parameters, double fitness) {
	if (population->size == population->max_size) return;

	if (population->individuals[population->size] == NULL) {
		population->individuals[population->size] = (double*)malloc(sizeof(double) * population->number_parameters);
	}
	memcpy(population->individuals[population->size], parameters, sizeof(double) * population->number_parameters);
	population->fitness[population->size] = fitness;
	population->size++;
}

void insert_incremental_info(POPULATION* population, double* parameters, double fitness, char* os_name, char* app_version) {
	if (population->size == population->max_size) return;

	if (population->individuals[population->size] == NULL) population->individuals[population->size] = (double*)malloc(sizeof(double) * population->number_parameters);
	memcpy(population->individuals[population->size], parameters, sizeof(double) * population->number_parameters);

	population->fitness[population->size] = fitness;
	if (population->app_versions[population->size] == NULL) population->app_versions[population->size] = (char*)malloc(sizeof(char) * 512);
	if (app_version != NULL) sprintf(population->app_versions[population->size], "%s", app_version);

	if (population->os_names[population->size] == NULL) population->os_names[population->size] = (char*)malloc(sizeof(char) * 512);
	if (os_name != NULL) sprintf(population->os_names[population->size], "%s", os_name);
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

void insert_sorted_info(POPULATION* population, double* parameters, double fitness, char *os_name, char* app_version) {
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

			if (population->app_versions[i] == NULL) population->app_versions[i] = (char*)malloc(sizeof(char) * 512);
			if (app_version != NULL) sprintf(population->app_versions[i], "%s", app_version);

			if (population->os_names[i] == NULL) population->os_names[i] = (char*)malloc(sizeof(char) * 512);
			if (os_name != NULL) sprintf(population->os_names[i], "%s", os_name);
			break;
		}
        }
}

/********
	*	Functions for removing individuals
 ********/

void remove_individual(POPULATION* population, int position) {
	free(population->individuals[position]);
	population->individuals[position] = NULL;
}

void remove_incremental(POPULATION* population, int position) {
	int i, j;
	for (i = position; i < population->size-1; i++) {
		population->fitness[i] = population->fitness[i+1];
		for (j = 0; j < population->number_parameters; j++) {
			population->individuals[i][j] = population->individuals[i+1][j];
			if (population->os_names[i] != NULL && population->os_names[i+1] != NULL) sprintf(population->os_names[i], population->os_names[i+1]);
			if (population->app_versions[i] != NULL && population->app_versions[i+1] != NULL) sprintf(population->app_versions[i], population->app_versions[i+1]);
		}
	}
	free(population->individuals[i]);
	population->individuals[i] = NULL;
	population->fitness[i] = -1;
	population->size--;
}

void remove_sorted(POPULATION *population, int position) {
	remove_incremental(population, position);
}

/********
	*	Get n distinct individuals from the population as a new population
 ********/
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

/********
	*	Functions for reading/writing populations
 ********/

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
		(*population)->app_versions[position] = (char*)malloc(sizeof(char) * 512);
		(*population)->os_names[position] = (char*)malloc(sizeof(char) * 512);

		fscanf(file, ", %s ", (*population)->os_names[position]);
		fgets((*population)->app_versions[position], 512, file);

		for (j = strlen((*population)->app_versions[position]); j >= 0; j--) {
			if ((*population)->app_versions[position][j] == '\0') {
			} else if ((*population)->app_versions[position][j] == ' ' || (*population)->app_versions[position][j] == 10 || (*population)->app_versions[position][j] == 13) {
				(*population)->app_versions[position][j] = '\0';
			} else {
				break;
			}
		}
	}
	(*population)->size = size;
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

int fwrite_individual(FILE *file, POPULATION *population, int position) {
	int j;
	fprintf(file, "[%d] %.20lf :", position, population->fitness[position]);
	for (j = 0; j < population->number_parameters; j++) {
		fprintf(file, " %.20lf", population->individuals[position][j]);
	}
	fprintf(file, ", %s %s\n", population->os_names[position], population->app_versions[position]);
	return 0;
}

int fwrite_population(FILE *file, POPULATION *population) {
	int i, count;

	fprintf(file, "size: %d, max_size: %d\n", population->size, population->max_size);
	fprintf(file, "number_parameters: %d\n", population->number_parameters);
	for (count = 0, i = 0; count < population->size && i < population->max_size; i++) {
		if (population->individuals[i] == NULL) continue;
		fwrite_individual(file, population, i);
		count++;
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
