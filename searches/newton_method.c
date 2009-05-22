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

