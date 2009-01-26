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
	double *diff, *sorted_diff, median_diff;
	double distance;
	int i, j, k;

	/********
		*	calculate sqrt sum difference/(distance)^2 / (N-1)
	 ********/
	diff = (double*)malloc(sizeof(double) * p->size);
	sorted_diff = (double*)malloc(sizeof(double) * p->size);
	for (i = 0; i < p->size; i++) diff[i] = 0.0;
	for (i = 0; i < p->size; i++) {
		for (j = i+1; j < p->size; j++) {
			distance = 0;
			for (k = 0; k < p->number_parameters; k++) {
				distance += fabs(p->individuals[i][k] - p->individuals[j][k]);
			}
			diff[i] += fabs(p->fitness[i] - p->fitness[j]) / distance;
			diff[j] += fabs(p->fitness[i] - p->fitness[j]) / distance;
		}
		diff[i] /= p->size-1;
		sorted_diff[i] = diff[i];
	}
	qsort(sorted_diff, p->size, sizeof(double), double_compare);
	median_diff = sorted_diff[p->size/2];

	printf("removing individuals, median_diff: %.20lf\n", median_diff);
	for (i = 0; i < p->size; i++) {
		fwrite_individual(stdout, p, i);
		printf("diff[%d]: %.20lf", i, diff[i]);
		if (diff[i] > range * median_diff) {
			printf(" -- REMOVED");

			p->size--;
			p->fitness[i] = p->fitness[p->size];
			for (j = 0; j < p->number_parameters; j++) p->individuals[i][j] = p->individuals[p->size][j];

			if (p->os_names[p->size] == NULL) p->os_names[i] = NULL;
			else sprintf(p->os_names[i], "%s", p->os_names[p->size]);

			if (p->app_versions[p->size] == NULL) p->app_versions[i] = NULL;
			else sprintf(p->app_versions[i], "%s", p->app_versions[p->size]);
			i--;
		}
		printf("\n");
	}
	free(diff);
	free(sorted_diff);

	return 0;
}
