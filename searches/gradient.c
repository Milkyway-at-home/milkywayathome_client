#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gradient.h"
#include "../evaluation/evaluator.h"
#include "../util/settings.h"

void get_gradient(int number_parameters, double *point, double *step, double *gradient) {
	int j;
	double e1, e2;
	for (j = 0; j < number_parameters; j++) {
		point[j] += step[j];
		e1 = evaluate(point);
		point[j] -= step[j] + step[j];
		e2 = evaluate(point);
		point[j] += step[j];

		gradient[j] = (e1 - e2)/(step[j] + step[j]);
	}
}

int gradient_below_threshold(int number_parameters, double* gradient, double threshold) {
	int i;
	for (i = 0; i < number_parameters; i++) {
		if (gradient[i] > threshold || gradient[i] < -threshold) return 0;
	}
	return 1;
}
