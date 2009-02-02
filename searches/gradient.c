#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gradient.h"
#include "../evaluation/evaluator.h"
#include "../util/settings.h"

void get_gradient(int number_parameters, double *point, double *step, double *gradient) {
	int j;
	double e1, e2, pj;
	for (j = 0; j < number_parameters; j++) {
		pj = point[j];
		point[j] = pj + step[j];
		e1 = evaluate(point);
		point[j] = pj - step[j];
		e2 = evaluate(point);
		point[j] = pj;

		gradient[j] = (e1 - e2)/(step[j] + step[j]);
		printf("\t\tgradient[%d]: %.20lf, (%.20lf - %.20lf)/(2 * %lf)\n", j, gradient[j], e1, e2, step[j]);
	}
}

int gradient_below_threshold(int number_parameters, double* gradient, double threshold) {
	int i;
	for (i = 0; i < number_parameters; i++) {
		if (gradient[i] > threshold || gradient[i] < -threshold) return 0;
	}
	return 1;
}
