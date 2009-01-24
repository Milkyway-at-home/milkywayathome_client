#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int remove_outliers(POPULATION *p, double range) {
	double *median_diff, *median_diff2, *fitness;
	double median, median2;
	double **hessian, **hessian_error, *gradient, *gradient_error, c, c_error;
	double *errors, avg_error;
	POPULATION *sample_population;
	int i, j, current;

	median_diff2 = (double*)malloc(sizeof(double) * p->size);
	median_diff = (double*)malloc(sizeof(double) * p->size);
	fitness = (double*)malloc(sizeof(double) * p->size);

	memcpy(fitness, p->fitness, sizeof(double) * p->size);

	mergesort(fitness, p->size, sizeof(double), double_compare);
	median = fitness[p->size/2];
	for (i = 0; i < p->size; i++) {
		median_diff[i] = fabs(p->fitness[i] - median);
	}

	memcpy(median_diff2, median_diff, sizeof(double) * p->size);
	mergesort(median_diff2, p->size, sizeof(double), double_compare);
	median2 = median_diff2[p->size/2];
	for (i = 0; i < p->size; i++) {
		median_diff2[i] = fabs(median_diff[i] - median2);
	}

        init_matrix(&hessian, p->number_parameters, p->number_parameters);
        init_matrix(&hessian_error, p->number_parameters, p->number_parameters);
        gradient = (double*)malloc(sizeof(double) * p->number_parameters);
        gradient_error = (double*)malloc(sizeof(double) * p->number_parameters);
                
	new_population(p->size, p->number_parameters, &sample_population);
	sample_population->size = p->size - 1;
	for (i = 0; i < p->size-1; i++) {
		sample_population->fitness[i] = p->fitness[i+1];
		if (sample_population->individuals[i] == NULL) sample_population->individuals[i] = (double*)malloc(sizeof(double) * p->number_parameters);
		for (j = 0; j < p->number_parameters; j++) sample_population->individuals[i][j] = p->individuals[i+1][j];
	}

	errors = (double*)malloc(sizeof(double) * p->size);
	avg_error = 0;
	for (i = 0; i < p->size; i++) {
	        parabolic_2d_regression_error(sample_population->size, sample_population->number_parameters, sample_population->individuals, sample_population->fitness, hessian, hessian_error, gradient, gradient_error, &c, &c_error);
		errors[i] = fabs(p->fitness[i] - parabolic_2d_evaluate(p->number_parameters, p->individuals[i], hessian, gradient, c));

		if (i < p->size - 1) {
			sample_population->fitness[i] = p->fitness[i];
			for (j = 0; j < p->number_parameters; j++) sample_population->individuals[i][j] = p->individuals[i][j];
		}
		avg_error += errors[i];
	}
	avg_error /= p->size;

	printf("errors: median [%.20lf], median of median [%.20lf], error[avg: %.20lf], c: [%.20lf]\n", median, median2, avg_error, c);
	for (i = 0; i < p->size; i++) {
		printf("\t[%d] %.20lf, %.20lf, %.20lf", i, median_diff[i], median_diff2[i], errors[i]);
		if (errors[i] >= range * avg_error) printf(" -- REMOVED");
		printf("\n");
	}

	current = 0;
	for (i = 0; i < p->size; i++) {
		if (errors[i] < range * avg_error) {
			p->fitness[current] = p->fitness[i];
			for (j = 0; j < p->number_parameters; j++) p->individuals[current][j] = p->individuals[i][j];
			current++;
		}
	}
	p->size = current;

	free_population(sample_population);
	free(sample_population);

        free_matrix(&hessian, p->number_parameters, p->number_parameters);
        free_matrix(&hessian_error, p->number_parameters, p->number_parameters);
        free(gradient);
        free(gradient_error);

	free(median_diff);
	free(median_diff2);
	free(errors);

	return 0;
}
