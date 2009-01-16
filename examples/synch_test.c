#include <stdio.h>
#include <stdlib.h>

#include "../searches/synchronous_newton_method.h"
#include "../evaluation/simple_evaluator.h"

double *min_parameters, *max_parameters;
int number_parameters;


double sum_of_squares(double* parameters) {
	int i;
	double sum;

	sum = 0.0;
	for (i = 0; i < number_parameters; i++) {
		sum += parameters[i] * parameters[i];
	}
	return -sum;
}



int main(int number_arguments, char** arguments) {
	int i;
	double* step;
	double* point;

	if (number_arguments != 3) {
		fprintf(stderr, "Invalid Arguments, proper usage:\n");
		fprintf(stderr, "\tsynch_test <number_parameters> <search_parameters>\n");

		return 0;
	}

	srand48(23);

	number_parameters = atoi(arguments[1]);
	min_parameters = (double*)malloc(sizeof(double) * number_parameters);
	max_parameters = (double*)malloc(sizeof(double) * number_parameters);
	step = (double*)malloc(sizeof(double) * number_parameters);
	point = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		min_parameters[i] = -100.0;
		max_parameters[i] = 100.0;
		step[i] = 2.0;
		point[i] = (drand48() * (double)(max_parameters[i] - min_parameters[i])) + (double)min_parameters[i];
	}

	init_simple_evaluator(sum_of_squares);
        printf("searching...\n");
	if (arguments[2][0] == 'n') {
                newton_method(number_parameters, point, step, 10);
	} else if (arguments[2][0] == 'r') {
                randomized_newton_method(number_parameters, point, step, 0, 200, 10, 0, 0);
        } else if (arguments[2][0] == 'l') {
		randomized_newton_method(number_parameters, point, step, 1, 200, 2, 200, 3);
	}

	return 0;
}
