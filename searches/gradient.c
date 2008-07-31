#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gradient.h"
#include "search.h"

void fprintf_gradient(FILE* file, GRADIENT* gradient) {
	int i;
	fprintf(file, "gradient:");

	for (i = 0; i < gradient->number_parameters; i++) {
		fprintf(file, " %lf", gradient->values[i]);
	}
	fprintf(file, "\n");
}

void create_gradient(double* point, double* step, int iteration, int number_parameters, GRADIENT** gradient) {
	int i, j;
	size_t param_size;

	param_size = sizeof(double) * number_parameters;

	(*gradient) = (GRADIENT*)malloc(sizeof(GRADIENT));
	(*gradient)->number_parameters = number_parameters;
	(*gradient)->set_values = 0;
	(*gradient)->iteration = iteration;


	(*gradient)->step = (double*)malloc(param_size);
	memcpy((*gradient)->step, step, param_size);

	(*gradient)->point = (double*)malloc(param_size);
	memcpy((*gradient)->point, point, param_size);

	(*gradient)->set_evaluations = (int**)malloc(sizeof(int*) * number_parameters);
	(*gradient)->evaluations = (double**)malloc(sizeof(double*) * number_parameters);
	(*gradient)->values = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		(*gradient)->set_evaluations[i] = (int*)malloc(sizeof(int*) * 2);
		(*gradient)->evaluations[i] = (double*)malloc(sizeof(double*) * 2);
		(*gradient)->values[i] = 0;
		for (j = 0; j < 2; j++) {
			(*gradient)->set_evaluations[i][j] = 0;
			(*gradient)->evaluations[i][j] = 0;
		}
	}

}

void free_gradient(GRADIENT* gradient) {
	int i;
	free(gradient->step);
	for (i = 0; i < gradient->number_parameters; i++) {
		free(gradient->set_evaluations[i]);
		free(gradient->evaluations[i]);
	}
	free(gradient->set_evaluations);
	free(gradient->evaluations);
	free(gradient->values);
}

int gradient__complete(GRADIENT* gradient) {
	return gradient->set_values == gradient->number_parameters;
}

void gradient__insert_individuals(GRADIENT* gradient, int number_individuals, double* fitness, char** metadata) {
	int i, iteration;
	int x, y;


	for (i = 0; i < number_individuals; i++) {
		sscanf(metadata[i], "iteration: %d, x: %d, y: %d", &iteration, &x, &y);
		if (gradient->set_evaluations[x][y] == 1 || gradient->iteration != iteration) continue;
		else {
			gradient->set_evaluations[x][y] = 1;
			gradient->evaluations[x][y] = fitness[i];

			if ((gradient->set_evaluations[x][0] + gradient->set_evaluations[x][1]) == 2) {
				gradient->values[x] = (gradient->evaluations[x][0] - gradient->evaluations[x][1])/(2*gradient->step[x]);
				gradient->set_values++;
			}
		}
	}
}

void gradient__get_individuals(GRADIENT* gradient, int number_individuals, double*** parameters, char*** metadata) {
	int i, j, current_individual;
	(*parameters) = (double**)malloc(sizeof(double*) * number_individuals);
	(*metadata) = (char**)malloc(sizeof(char*) * number_individuals);

	current_individual = 0;
	while (current_individual != number_individuals) {
		for (i = 0; i < gradient->number_parameters; i++) {
			for (j = 0; j < 2; j++) {
				if (!gradient->set_evaluations[i][j]) {
					(*parameters)[current_individual] = (double*)malloc(sizeof(double) * gradient->number_parameters);
					memcpy((*parameters)[current_individual], gradient->point, sizeof(double) * gradient->number_parameters);
					(*metadata)[current_individual] = (char*)malloc(sizeof(char) * METADATA_SIZE);
					sprintf((*metadata)[current_individual], "iteration: %d, x: %d, y: %d", gradient->iteration, i, j);
					switch (j) {
						case 0:
							(*parameters)[current_individual][i] += gradient->step[i];
						break;
						case 1:
							(*parameters)[current_individual][i] -= gradient->step[i];
						break;
					}
					current_individual++;
					if (current_individual == number_individuals) return;
				}
			}
		}
	}
}
