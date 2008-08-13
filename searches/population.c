#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

#include "population.h"

POPULATION* new_population(char search_path[512], char search_parameters[512], double *min_parameters, double *max_parameters, int number_parameters, int population_size, int max_evaluations) {
	POPULATION* population;
	int i;

	population = (POPULATION*)malloc(sizeof(POPULATION));
	population->search_path = (char*)malloc(sizeof(char) * 512);
	population->search_parameters = (char*)malloc(sizeof(char) * 512);
	strcpy(population->search_path, search_path);
	strcpy(population->search_parameters, search_parameters);

	population->current_size = 0;
	population->max_size = population_size;
	population->current_evaluation = 0;
	population->max_evaluations = max_evaluations;

	population->number_parameters = number_parameters;
	population->min_parameters = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) population->min_parameters[i] = min_parameters[i];
	population->max_parameters = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) population->max_parameters[i] = max_parameters[i];

	population->fitness = (double*)malloc(sizeof(double) * population_size);
	population->individuals = (double**)malloc(sizeof(double*) * population_size);
	i = 0;

	write_population(search_path, population);

	return population;
}


POPULATION* fread_population(FILE* file) {
	POPULATION *population;
	int i, j;

	population = (POPULATION*)malloc(sizeof(POPULATION));

	fscanf(file, "%s\n", population->search_path);
	fscanf(file, "%s\n", population->search_parameters);

	fscanf(file, "current size: %d\n", &population->current_size);
	fscanf(file, "max size: %d\n", &population->max_size);
	fscanf(file, "current evaluation: %d\n", &population->current_evaluation);
	fscanf(file, "max evaluations: %d\n", &population->max_evaluations);
	fscanf(file, "number parameters: %d\n", &population->number_parameters);

	population->min_parameters = (double*)malloc(sizeof(double) * population->number_parameters);
	fscanf(file, "min_parameters:");
	for (i = 0; i < population->number_parameters; i++) {
		fscanf(file, " %lf", &population->min_parameters[i]);
	}
	fscanf(file, "\n");

	population->max_parameters = (double*)malloc(sizeof(double) * population->number_parameters);
	fscanf(file, "max_parameters:");
	for (i = 0; i < population->number_parameters; i++) {
		fscanf(file, " %lf", &population->max_parameters[i]);
	}
	fscanf(file, "\n");

	population->fitness = (double*)malloc(sizeof(double) * population->current_size);
	population->individuals = (double**)malloc(sizeof(double*) * population->current_size);
	fscanf(file, "population: <fitness> : <parameters>\n");
	for (i = 0; i < population->current_size; i++) {
		population->individuals[i] = (double*)malloc(sizeof(double) * population->number_parameters);
		fscanf(file, "%lf :", &population->fitness[i]);
		for (j = 0; j < population->number_parameters; j++) {
			fscanf(file, " %lf", &population->individuals[i][j]);
		}
		fscanf(file, "\n");
	}

	return population;
}


void fwrite_population(FILE *file, POPULATION *population) {
	int i, j;

	fprintf(file, "%s\n", population->search_path);
	fprintf(file, "%s\n", population->search_parameters);
	fprintf(file, "current size: %d\n", population->current_size);
	fprintf(file, "max size: %d\n", population->max_size);
	fprintf(file, "current evaluation: %d\n", population->current_evaluation);
	fprintf(file, "max evaluations: %d\n", population->max_evaluations);
	fprintf(file, "number parameters: %d\n", population->number_parameters);

	fprintf(file, "min_parameters:");
	for (i = 0; i < population->number_parameters; i++) {
		fprintf(file, " %lf", population->min_parameters[i]);
	}
	fprintf(file, "\n");

	fprintf(file, "max_parameters:");
	for (i = 0; i < population->number_parameters; i++) {
		fprintf(file, " %lf", population->max_parameters[i]);
	}
	fprintf(file, "\n");

	fprintf(file, "population: <fitness> : <parameters>\n");
	for (i = 0; i < population->current_size; i++) {
		fprintf(file, "%lf :", population->fitness[i]);
		for (j = 0; j < population->number_parameters; j++) {
			fprintf(file, " %lf", population->individuals[i][j]);
		}
		fprintf(file, "\n");
	}
	fflush(file);
}


void write_population(char path[512], POPULATION *population) {
	FILE *file;
	char search_path[1024];
	int result;

	result = mkdir(path, 0777);
	if (result == -1) {
		if (errno == EACCES) {
			fprintf(stderr, "ERROR CREATING DIRECTORY: %s\n\tWrite permission denied.\n", path);
			return;
		} else if (errno == EEXIST) {
			fprintf(stderr, "ERROR CREATING DIRECTORY: %s\n\tFile name already exists.\n", path);
		} else if (errno == EMLINK) {
			fprintf(stderr, "ERROR CREATING DIRECTORY: %s\n\tParent directory has too many entries.\n", path);
			return;
		} else if (errno == ENOSPC) {
			fprintf(stderr, "ERROR CREATING DIRECTORY: %s\n\tNot enough space on file system to create directory.\n", path);
			return;
		} else if (errno == EROFS) {
			fprintf(stderr, "ERROR CREATING DIRECTORY: %s\n\tParent directory is read-only.\n", path);
			return;
		} else {
			fprintf(stderr, "ERROR CREATING DIRECTORY: %s\n\tUnknown error!?!?\n", path);
			return;
		}
	}

	sprintf(search_path, "%s/population", path);
	file = fopen(search_path, "w");
	fwrite_population(file, population);
	fclose(file);
}

POPULATION* read_population(char path[512]) {
	FILE *file;
	POPULATION *population;
	char search_path[1024];

	sprintf(search_path, "%s/population", path);
	file = fopen(search_path, "r");
	population = fread_population(file);
	fclose(file);
	return population;
}

void fwrite_population_statistics(FILE* file, POPULATION *population) {
	double min_fitness, avg_fitness, max_fitness;
	int i;

	min_fitness = population->fitness[0];
	max_fitness = population->fitness[0];
	avg_fitness = 0;
	for (i = 0; i < population->current_size; i++) {
		avg_fitness += population->fitness[i];
		if (min_fitness > population->fitness[i]) min_fitness = population->fitness[i];
		else if (max_fitness < population->fitness[i]) max_fitness = population->fitness[i];
	}
	avg_fitness /= population->current_size;

	fprintf(file, "min: %lf, avg: %lf, max: %lf", min_fitness, avg_fitness, max_fitness);
}


void insert_sorted(POPULATION* population, double* parameters, double fitness) {
        int insert_position, i, j;
        insert_position = 0;
        for (i = 0; i < population->current_size; i++) {
                if (fitness < population->fitness[i]) break;
        }
        insert_position = i;

        population->current_evaluation++;
        if (insert_position == population->max_size) return;
                
        if (population->current_size == population->max_size) free(population->individuals[population->max_size-1]);
        else {
                population->current_size++;
                population->individuals = (double**)realloc(population->individuals, sizeof(double*) * population->current_size);
                population->fitness = (double*)realloc(population->fitness, sizeof(double) * population->current_size);
        } 
                
        for (i = population->current_size - 1; i > insert_position; i--) {
                population->individuals[i] = population->individuals[i-1];
                population->fitness[i] = population->fitness[i-1];
        }

        population->individuals[i] = (double*)malloc(sizeof(double) * population->number_parameters);
	for (j = 0; j < population->number_parameters; j++) population->individuals[i][j] = parameters[j];
        population->fitness[i] = fitness;
}


void replace(POPULATION* population, int position, double* parameters, double fitness) {
	int i;
	if (population->individuals[position] == NULL) {
		population->individuals[position] = (double*)malloc(sizeof(double) * population->number_parameters);
		population->current_size++;
	}

	population->fitness[position] = fitness;
	for (i = 0; i < population->number_parameters; i++) population->individuals[position][i] = parameters[i];
}


void replace_if_better(POPULATION* population, int position, double* parameters, double fitness) {
	int i;
	if (population->fitness[position] > fitness || population->individuals[position] == NULL) {
		if (population->individuals[position] == NULL) {
			population->individuals[position] = (double*)malloc(sizeof(double) * population->number_parameters);
			population->current_size++;
		}

		population->fitness[position] = fitness;
		for (i = 0; i < population->number_parameters; i++) population->individuals[position][i] = parameters[i];
	}
}

void get_n_distinct(POPULATION *population, int number_parents, double ***parents, double **parent_fitness) {
	int i, j, k, target;
	int *parent_positions;

	(*parent_fitness) = (double*)malloc(sizeof(double) * number_parents);
	(*parents) = (double**)malloc(sizeof(double*) * number_parents);
	parent_positions = (int*)malloc(sizeof(int) * number_parents);

	for (i = 0; i < number_parents; i++) {
		target = drand48() * (population->max_size - i);
		for (j = 0; j < i; j++) {
			for (k = 0; k < i; k++) {
				if (target == parent_positions[k]) target++;
			}
		}
		parent_positions[i] = target;
	}

	for (i = 0; i < number_parents; i++) {
		(*parents)[i] = (double*)malloc(sizeof(double) * population->number_parameters); 
		for (j = 0; j < population->number_parameters; j++) (*parents)[i][j] = population->individuals[parent_positions[i]][j];
		(*parent_fitness)[i] = population->fitness[parent_positions[i]];
	}
}

void bound_parameters(POPULATION *population, double* parameters) {
	int i;

	for (i = 0; i < population->number_parameters; i++) {
		if (parameters[i] > population->max_parameters[i]) parameters[i] = population->max_parameters[i];
		else if (parameters[i] < population->min_parameters[i]) parameters[i] = population->min_parameters[i];
	}
}
