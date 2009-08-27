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



void get_pair_sum(double pair_weight, double *parent, double **pairs, int number_pairs, int number_parameters, double *result) {
	int i, j;

	for (i = 0; i < number_parameters; i++) {
		result[i] = 0;
		for (j = 0; j < number_pairs; j++) {
			result[i] += pairs[(j*2)][i] - pairs[(j*2)+1][i];
		}
		result[i] = parent[i] + pair_weight * result[i];
	}
}

void binomial_recombination(double crossover_rate, double *p1, double *p2, int number_parameters, double *result) {
	int i, selected;
	selected = (int)(dsfmt_gv_genrand_close_open() * (double)number_parameters);

	for (i = 0; i < number_parameters; i++) {
		if (i == selected || dsfmt_gv_genrand_close_open() < crossover_rate) {
			result[i] = p2[i];
		} else {
			result[i] = p1[i];
		}
	}
}

void exponential_recombination(double crossover_rate, double *p1, double *p2, int number_parameters, double *result) {
	int i, selected;
	selected = (int)(dsfmt_gv_genrand_close_open() * (double)number_parameters);

	for (i = 0; i < number_parameters; i++) {
		if (i == selected || dsfmt_gv_genrand_close_open() < crossover_rate) break;
		result[i] = p1[i];
	}

	for (; i < number_parameters; i++) result[i] = p2[i];
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

	worst_parent = 0;
	for (i = 1; i < number_parents; i++) {
		if (fitness[i] > fitness[worst_parent]) worst_parent = i;
	}

	/********
		*	Calcualte the centroid.
	 ********/
	for (i = 0; i < number_parents; i++) {
		if (i == worst_parent) continue;
		for (j = 0; j < number_parameters; j++) {
			parameters[j] += parents[i][j] / (number_parents - 1);
		}
	}

	/********
		*	Generate points along the line from the worst parent to the centroid
	 ********/
	distance = ls_center + dsfmt_gv_genrand_close_open() * (ls_outside - ls_center);
	for (j = 0; j < number_parameters; j++) {
		parameters[j] = parameters[j] + distance * (parameters[j] - parents[worst_parent][j]);
	}
	return distance;
}
