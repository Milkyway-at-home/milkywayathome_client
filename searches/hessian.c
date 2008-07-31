#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hessian.h"
#include "search.h"

void fprintf_hessian(FILE* file, HESSIAN* hessian) {
	int i, j;
	fprintf(file, "hessian:\n");
	for (i = 0; i < hessian->number_parameters; i++) {
		for (j = 0; j < hessian->number_parameters; j++) {
			fprintf(file, " %lf", hessian->values[i][j]);
		}
		fprintf(file, "\n");
	}
}


void create_hessian(double* point, double* step, int iteration, int number_parameters, HESSIAN** hessian) {
	int i, j, k;
	size_t param_size;

	param_size = sizeof(double) * number_parameters;

	(*hessian) = (HESSIAN*)malloc(sizeof(HESSIAN));
	(*hessian)->number_parameters = number_parameters;
	(*hessian)->set_values = 0;
	(*hessian)->iteration = iteration;


	(*hessian)->step = (double*)malloc(param_size);
	memcpy((*hessian)->step, step, param_size);

	(*hessian)->point = (double*)malloc(param_size);
	memcpy((*hessian)->point, point, param_size);

	(*hessian)->set_evaluations = (int***)malloc(sizeof(int**) * number_parameters);
	(*hessian)->evaluations = (double***)malloc(sizeof(double**) * number_parameters);
	(*hessian)->values = (double**)malloc(sizeof(double*) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		(*hessian)->set_evaluations[i] = (int**)malloc(sizeof(int*) * number_parameters);
		(*hessian)->evaluations[i] = (double**)malloc(sizeof(double*) * number_parameters);
		(*hessian)->values[i] = (double*)malloc(param_size);
		for (j = 0; j < number_parameters; j++) {
			(*hessian)->set_evaluations[i][j] = (int*)malloc(sizeof(int) * 4);
			(*hessian)->evaluations[i][j] = (double*)malloc(sizeof(double) * 4);
			for (k = 0; k < 4; k++) {
				(*hessian)->set_evaluations[i][j][k] = 0;
				(*hessian)->evaluations[i][j][k] = 0;
			}
		}
	}

}

void free_hessian(HESSIAN* hessian) {
	int i, j;
	free(hessian->step);
	for (i = 0; i < hessian->number_parameters; i++) {
		for (j = 0; j < hessian->number_parameters; j++) {
			free(hessian->set_evaluations[i][j]);
			free(hessian->evaluations[i][j]);
		}
		free(hessian->set_evaluations[i]);
		free(hessian->evaluations[i]);
		free(hessian->values[i]);
	}
	free(hessian->set_evaluations);
	free(hessian->evaluations);
	free(hessian->values);
}

int hessian__complete(HESSIAN* hessian) {
	return hessian->set_values == (hessian->number_parameters * hessian->number_parameters);
}

void hessian__insert_individuals(HESSIAN* hessian, int number_individuals, double* fitness, char** metadata) {
	int i, iteration;
	int x, y, z;
	double first, second;

	for (i = 0; i < number_individuals; i++) {
		sscanf(metadata[i], "iteration: %d, x: %d, y: %d, z: %d", &iteration, &x, &y, &z);
		if (hessian->set_evaluations[x][y][z] == 1 || hessian->iteration != iteration) continue;
		else {
			hessian->set_evaluations[x][y][z] = 1;
			hessian->evaluations[x][y][z] = fitness[i];
	
			if ((hessian->set_evaluations[x][y][0] + hessian->set_evaluations[x][y][1] + hessian->set_evaluations[x][y][2] + hessian->set_evaluations[x][y][3]) == 4) {
				first = (hessian->evaluations[x][y][0] - hessian->evaluations[x][y][1])/(2*hessian->step[y]);
				second = (hessian->evaluations[x][y][2] - hessian->evaluations[x][y][3])/(2*hessian->step[y]);
			
				hessian->values[x][y] = (first - second) / (2 * hessian->step[x]);
				hessian->set_values++;
			}
		}
	}
}

void hessian__get_individuals(HESSIAN* hessian, int number_individuals, double*** parameters, char*** metadata) {
	int i, j, k, current_individual;
	(*parameters) = (double**)malloc(sizeof(double*) * number_individuals);
	(*metadata) = (char**)malloc(sizeof(char*) * number_individuals);

	current_individual = 0;
	while (current_individual != number_individuals) {
		for (i = 0; i < hessian->number_parameters; i++) {
			for (j = 0; j < hessian->number_parameters; j++) {
				for (k = 0; k < 4; k++) {
					if (!hessian->set_evaluations[i][j][k]) {
						(*parameters)[current_individual] = (double*)malloc(sizeof(double) * hessian->number_parameters);
						memcpy((*parameters)[current_individual], hessian->point, sizeof(double) * hessian->number_parameters);
						(*metadata)[current_individual] = (char*)malloc(sizeof(char) * METADATA_SIZE);
						sprintf((*metadata)[current_individual], "iteration: %d, x: %d, y: %d, z: %d", hessian->iteration, i, j, k);
						switch (k) {
							case 0:
								(*parameters)[current_individual][i] += hessian->step[i];
								(*parameters)[current_individual][j] += hessian->step[j];
							break;
							case 1:
								(*parameters)[current_individual][i] += hessian->step[i];
								(*parameters)[current_individual][j] -= hessian->step[j];
							break;
							case 2:
								(*parameters)[current_individual][i] -= hessian->step[i];
								(*parameters)[current_individual][j] += hessian->step[j];
							break;
							case 3:
								(*parameters)[current_individual][i] -= hessian->step[i];
								(*parameters)[current_individual][j] -= hessian->step[j];
							break;
						}
						current_individual++;
						if (current_individual == number_individuals) return;
					}
				}
			}
		}
	}
}
