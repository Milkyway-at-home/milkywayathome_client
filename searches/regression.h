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

#ifndef FGDO_REGRESSION
#define FGDO_REGRESSION

double parabolic_2d_evaluate(int number_points, double *x, double **hessian, double *gradient, double c);

void parabolic_2d_regression(int number_points, int x_length, double **x, double *y, double **hessian, double *gradient, double *c);
void parabolic_2d_regression_error(int number_points, int x_length, double **x, double *y, double **hessian, double **hessian_error, double *gradient, double *gradient_error, double *c, double *c_error);

double parabolic_center(double a, double b, double c);
void parabolic_regression(int number_points, double *x, double *y, double *a, double *b, double *c, double *a_error, double *b_error, double *c_error);

#endif
