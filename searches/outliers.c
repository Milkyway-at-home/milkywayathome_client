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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "population.h"
#include "regression.h"
#include "../util/matrix.h"

int double_compare(const void *p1, const void *p2) {
	const double *d1 = (const double*)p1;
	const double *d2 = (const double*)p2;

	if (d1[0] < d2[0]) return -1;  
	else if (d1[0] == d2[0]) return 0;
	else return 1;
}

double get_distance(int number_parameters, double f1, double *p1, double f2, double *p2) {
	int i;
	double distance = 0;

	for (i = 0; i < number_parameters; i++) {
		distance += fabs(p1[i] - p2[i]);
	}
	return fabs(f1 - f2)/distance;
}

double distance_error_from(POPULATION *p, double fitness, double *parameters) {
	double distance;
	int i;

	distance = 0;
	for (i = 0; i < p->max_size; i++) {
		if (!individual_exists(p, i)) continue;
		distance += get_distance(p->number_parameters, p->fitness[i], p->individuals[i], fitness, parameters);
	}
	return distance / p->size;
}

int get_distance_errors(POPULATION *p, double **errors) {
	int i, j, current_i;
	double *e;
	(*errors) = (double*)malloc(sizeof(double) * p->size);
	e = (*errors);

	for (i = 0; i < p->size; i++) e[i] = 0.0;

	current_i = 0;
	for (i = 0; i < p->max_size; i++) {
		if (!individual_exists(p, i)) continue;

		for (j = 0; j < p->max_size; j++) {
			if (!individual_exists(p, j)) continue;
			if (i == j) continue;

			e[current_i] += get_distance(p->number_parameters, p->fitness[i], p->individuals[i], p->fitness[j], p->individuals[j]);

			if (isnan(e[current_i])) printf("ERROR: distance [%d] to [%d] resulted in nan\n", i, j);
		}
		current_i++;
	}
	for (i = 0; i < p->size; i++) e[i] /= current_i;
	return current_i;
}

void get_error_stats(double *errors, int error_size, double *min_error, double *max_error, double *median_error, double *average_error) {
	int i;
	double sum;

	sum = 0.0;
	(*min_error) = DBL_MAX;
	(*max_error) = DBL_MIN;
	for (i = 0; i < error_size; i++) {
		if (errors[i] < (*min_error)) (*min_error) = errors[i];
		if (errors[i] > (*max_error)) (*max_error) = errors[i];
		sum += errors[i];
	}
	(*average_error) = sum / error_size;
	(*median_error) = errors[error_size/2];
}

void population_error_stats(POPULATION *p, double *min_error, double *max_error, double *median_error, double *average_error) {
	double *errors;
	int error_size;

	error_size = get_distance_errors(p, &errors);
	get_error_stats(errors, error_size, min_error, max_error, median_error, average_error);
	free(errors);
}

#define REMOVE_OUTLIERS_INDIVIDUAL 1
#define REMOVE_OUTLIERS_INCREMENTAL 2
#define REMOVE_OUTLIERS_SORTED 3

void remove_outliers_helper(POPULATION *p, double range, int type) {
	double *errors, min_error, max_error, median_error, average_error;
	int i, error_size, current;

	error_size = get_distance_errors(p, &errors);
	get_error_stats(errors, error_size, &min_error, &max_error, &median_error, &average_error);

	printf("removing individuals, median_error: %.20lf, average_error: %.20lf\n", median_error, average_error);
	current = 0;
	for (i = 0; i < p->max_size; i++) {
		if (!individual_exists(p, i)) continue;
		printf("errors[%d]: %.20lf, fitness: %.20lf", i, errors[current], p->fitness[i]);
		if (errors[current] > range * average_error) {
			printf(" -- REMOVED");
			if (type == REMOVE_OUTLIERS_INDIVIDUAL) {
				printf(" INDIVIDUAL");
				remove_individual(p, i);
			} else if (type == REMOVE_OUTLIERS_INCREMENTAL || type == REMOVE_OUTLIERS_SORTED) {
				printf(" INCREMENTAL");
				remove_incremental(p, i);
				i--;
			}
		}
		printf("\n");
		current++;
	}
	free(errors);
}

void remove_outliers(POPULATION *p, double range) {
	remove_outliers_helper(p, range, REMOVE_OUTLIERS_INDIVIDUAL);
}

void remove_outliers_incremental(POPULATION *p, double range) {
	remove_outliers_helper(p, range, REMOVE_OUTLIERS_INCREMENTAL);
}

void remove_outliers_sorted(POPULATION *p, double range) {
	remove_outliers_helper(p, range, REMOVE_OUTLIERS_SORTED);
}

double get_distance2(POPULATION *p, int position, double *parameters) {
	double difference, distance, current_distance;
	int i, j;

	distance = 0;
	for (i = 0; i < p->max_size; i++) {
		if (!individual_exists(p, i) || i == position) continue;
		current_distance = 0;
		for (j = 0; j < p->number_parameters; j++) {
			difference = parameters[j] - p->individuals[i][j];
			current_distance += difference * difference;
		}
		distance += sqrt(current_distance);
	}
	return distance;
}

double get_error2(POPULATION *p, int position, double fitness) {
	double error, difference;
	int i;

	error = 0;
	for (i = 0; i < p->max_size; i++) {
		if (!individual_exists(p, i) || i == position) continue;
		difference = fitness - p->fitness[i];
		error += difference * difference;
	}
	return sqrt(error);
}

double distance_error2(POPULATION *p, int position, double fitness, double *parameters) {
	return get_error2(p, position, fitness) / get_distance2(p, position, parameters);
}

#ifdef OUTLIERS_MAIN
int main(int argc, char**argv) {
	POPULATION *p;
	double avg_distance, avg_error, distance, error;
	int i, j, count;

	read_population(argv[1], &p);

	for (i = 0; i < p->max_size; i++) {
		if (!individual_exists(p, i)) continue;

		count = 0;
		avg_distance = 0;
		avg_error = 0;
		for (j = 0; j < p->max_size; j++) {
			if (!individual_exists(p, j) || j == i) continue;
			distance = get_distance2(p, i, p->individuals[j]);
			error = get_error2(p, i, p->fitness[j]);

			avg_distance += distance;
			avg_error += error;
			count++;
		}
		avg_distance /= count;
		avg_error /= count;

		distance = get_distance2(p, i, p->individuals[i]);
		error = get_error2(p, i, p->fitness[i]);

		printf("distance: [%.15lf - %.15lf], error: [%.15lf - %.15lf], error/distance: [%.15lf - %.15lf], fitness: [%.15lf], - %.15lf\n", distance, avg_distance, error, avg_error, error/distance, avg_error/avg_distance, p->fitness[i], (error/distance)/(avg_error/avg_distance));
	}
}
#endif
