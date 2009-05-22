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

#ifndef GEM_GRADIENT_H
#define GEM_GRADIENT_H

#include <stdio.h>

typedef struct gradient {
	int iteration;
	int set_values;
	int number_parameters;
	double* step;
	double* point;

	int** set_evaluations;
	double** evaluations;
	double* values;
} GRADIENT;


void get_gradient__checkpointed(int number_parameters, double *point, double *step, double *gradient, char *checkpoint_file);
void get_gradient(int number_parameters, double *point, double *step, double *gradient);
int gradient_below_threshold(int number_parameters, double* gradient, double threshold);


#endif
