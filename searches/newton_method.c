#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "newton_method.h"

#include "../util/matrix.h"
#include "../util/io_util.h"

void newton_step_i(int number_parameters, double **hessian, double *gradient, double **step) {
	double **inverse_hessian;
	int i;

	matrix_invert(hessian, number_parameters, number_parameters, &inverse_hessian);
	vector_matrix_multiply(gradient, number_parameters, inverse_hessian, number_parameters, number_parameters, step);

	for (i = 0; i < number_parameters; i++) free(inverse_hessian[i]);
	free(inverse_hessian);
}

void newton_step(int number_parameters, double **hessian, double *gradient, double *step) {
	double **inverse_hessian;
	double *s;
	int i;

	matrix_invert(hessian, number_parameters, number_parameters, &inverse_hessian);
	vector_matrix_multiply(gradient, number_parameters, inverse_hessian, number_parameters, number_parameters, &s);

	for (i = 0; i < number_parameters; i++) {
		step[i] = s[i];
		free(inverse_hessian[i]);
	}
	free(s);
	free(inverse_hessian);
}

