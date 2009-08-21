/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../mersenne_twister/dSFMT.h"

#include "recombination.h"

void mutate(double* parent, double* min_parameters, double* max_parameters, int number_parameters, double *parameters) {
	int i, target;

	target = (int)(drand48() * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		if (i == target) {
			parameters[i] = min_parameters[i] + (max_parameters[i] - min_parameters[i]) * drand48();
		} else {
			parameters[i] = parent[i];
		}
	}
}


void random_recombination(int number_parameters, double *min_parameters, double *max_parameters, double *parameters) {
	int i;
	for (i = 0; i < number_parameters; i++) parameters[i] = min_parameters[i] + (max_parameters[i] - min_parameters[i]) * drand48();
}

void range_recombination(int number_parameters, double* point, double* range, double *result) {
	int i;
	for (i = 0; i < number_parameters; i++) {
		result[i] = point[i] + (((drand48() * 2.0) - 1.0) * range[i]);
	}
}

double random_linear_recombination(int number_parameters, double min, double max, double* initial, double* step, double *parameters) {
	int i;
	double value = min + ((max - min) * drand48());
	for (i = 0; i < number_parameters; i++) {
		parameters[i] = initial[i] + (value * step[i]);
	}
	return value;
}


double* get_pair_sum(double **individuals, int number_individuals, int number_parameters, int number_pairs, double scale) {
	int pair1, pair2;
	double* parameters;
	int i, j;

	parameters = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) parameters[i] = 0;
	for (i = 0; i < number_pairs; i++) {
		pair1 = (int)(drand48() * number_individuals);
		pair2 = (int)(drand48() * (number_individuals - 1));
		if (pair2 == pair1) pair2++;

		for (j = 0; j < number_parameters; j++) {
			parameters[j] += individuals[pair1] - individuals[pair2];
		}
	}
	if (scale < 0) scale = drand48();
	for (i = 0; i < number_parameters; i++) parameters[i] /= (scale/number_pairs);
	return parameters;
}

double* get_dir_sum(double **individuals, double *fitness, int number_individuals, int number_parameters, int number_pairs, double scale) {
	int pair1, pair2;
	double* parameters;
	int i, j, temp;

	parameters = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) parameters[i] = 0;
	for (i = 0; i < number_pairs; i++) {
		pair1 = (int)(drand48() * number_individuals);
		pair2 = (int)(drand48() * (number_individuals - 1));
		if (pair2 == pair1) pair2++;

		if (fitness[pair1] > fitness[pair2]) {
			temp = pair2;
			pair2 = pair1;
			pair1 = temp;
		}

		for (j = 0; j < number_parameters; j++) {
			parameters[j] += individuals[pair1] - individuals[pair2];
		}
	}
	if (scale < 0) scale = drand48();
	for (i = 0; i < number_parameters; i++) parameters[i] /= (scale/number_pairs);
	return parameters;
}


double* binomial_recombination(double **parents, int number_parents, int number_parameters, double crossover_rate, double crossover_scale) {
	int i, selected;
	double *parameters;
	parameters = (double*)malloc(sizeof(double) * number_parameters);

	selected = (int)(drand48() * number_parameters);
	if (crossover_scale < 0) crossover_scale = drand48();

	for (i = 0; i < number_parameters; i++) {
		if (i == selected || drand48() < crossover_rate) {
			parameters[i] = parents[0][i] + (parents[0][i] - parents[1][i]) * crossover_scale;
		} else {
			parameters[i] = parents[0][i];
		}
	}
	return parameters;
}

double* exponential_recombination(double **parents, int number_parents, int number_parameters, double crossover_rate, double crossover_scale) {
	int i, selected;
	double *parameters;
	parameters = (double*)malloc(sizeof(double) * number_parameters);

	selected = (int)(drand48() * number_parameters);
	if (crossover_scale < 0) crossover_scale = drand48();
	for (i = 0; i < selected; i++) {
		if (drand48() < crossover_rate) break;
		parameters[i] = parents[0][i];
	}
	for (; i < number_parameters; i++) {
		parameters[i] = parents[0][i] + (parents[0][i] - parents[0][i]) * crossover_scale;
	}
	return parameters;
}

void average_recombination(double** parents, int number_parents, int number_parameters, double* parameters) {
	int i, j;
	for (i = 0; i < number_parameters; i++) {
		for (j = 0; j < number_parents; j++) {
			parameters[i] = parents[j][i];
		}
		parameters[i] /= number_parents;
	}
}

void lower_recombination(double** parents, int number_parents, int number_parameters, double* parameters) {
	int i;
	average_recombination(parents, number_parents, number_parameters, parameters);

	for (i = 0; i < number_parents; i++) {
		parameters[i] = 2.0 * parents[0][i] - parameters[i];
	}
}

void higher_recombination(double** parents, int number_parents, int number_parameters, double* parameters) {
	int i;
	average_recombination(parents, number_parents, number_parameters, parameters);

	for (i = 0; i < number_parents; i++) {
		parameters[i] = 2.0 * parents[1][i] - parameters[i];
	}
}


double simplex_recombination(double** parents, double* fitness, int number_parents, int number_parameters, double ls_center, double ls_outside, double* parameters) {
	int i, j, worst_parent;
	double distance;

	printf("finding worst parent\n");
	worst_parent = 0;
	for (i = 1; i < number_parents; i++) {
		if (fitness[i] > fitness[worst_parent]) worst_parent = i;
		printf("\tworst parent: %d\n", worst_parent);
	}

	/********
		*	Calcualte the centroid.
	 ********/
	printf("calculating centroid\n");
	printf("sizeof parameters: %d\n", sizeof(parameters));
	for (i = 0; i < number_parents; i++) {
		if (i == worst_parent) continue;
		printf("parameters:");
		for (j = 0; j < number_parameters; j++) {
			parameters[j] += parents[i][j] / (number_parents - 1);
			printf(" %.15lf", parameters[j]);
		}
		printf("\n");
	}
	printf("dividing [%d]:\n", number_parameters);
//	for (j = 0; j < number_parameters; j++) {
//		parameters[j] /= number_parents - 1;
//		printf(" %.15lf\n", parameters[j]);
//	}

	/********
		*	Generate points along the line from the worst parent to the centroid
	 ********/
	printf("generating point\n");
	distance = ls_center + dsfmt_gv_genrand_close_open() * (ls_outside - ls_center);
	printf("calculating point from distance: %.15lf\n", distance);
	for (j = 0; j < number_parameters; j++) {
		parameters[j] = parameters[j] + distance * (parameters[j] - parents[worst_parent][j]);
	}
	printf("returning distance\n");
	return distance;
}
