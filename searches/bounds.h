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

#ifndef FGDO_BOUNDS
#define FGDO_BOUNDS

#include "stdio.h"

typedef struct bounds {
	int number_parameters;
	double *min_bound, *max_bound;
	int *in_radians;
} BOUNDS;

void new_bounds(BOUNDS **bounds, int number_parameters, double *min_bound, double *max_bound, int *in_radians);
void free_bounds(BOUNDS **bounds);

void bound_parameters(double *parameters, BOUNDS *bounds);
void bound_velocity(double *parameters, double *velocity, BOUNDS *bounds);

void fwrite_bounds(FILE *file, BOUNDS *bounds);
void fread_bounds(FILE *file, BOUNDS **bounds);

#endif

