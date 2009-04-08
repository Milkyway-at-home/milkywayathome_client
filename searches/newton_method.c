#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "newton_method.h"

#include "../util/matrix.h"
#include "../util/io_util.h"

void newton_step__alloc(int number_parameters, double **hessian, double *gradient, double **step) {
	double **inverse_hessian;

	matrix_invert__alloc(hessian, number_parameters, number_parameters, &inverse_hessian);
	matrix_vector_multiply__alloc(inverse_hessian, number_parameters, number_parameters, gradient, number_parameters, step);

	free_matrix(&inverse_hessian, number_parameters, number_parameters);
}

void newton_step(int number_parameters, double **hessian, double *gradient, double *step) {
	double **inverse_hessian;

	matrix_invert__alloc(hessian, number_parameters, number_parameters, &inverse_hessian);
	matrix_vector_multiply(inverse_hessian, number_parameters, number_parameters, gradient, number_parameters, step);

	free_matrix(&inverse_hessian, number_parameters, number_parameters);
}

