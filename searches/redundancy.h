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

#ifndef FGDO_REDUNDANCY_H
#define FGDO_REDUNDANCY_H

#include "stdio.h"

typedef struct redundancy REDUNDANCY;

struct redundancy {
	double *parameters;
	double *velocity;
	double fitness;

	REDUNDANCY *next;
};

int parameters_match(int number_parameters, double *p1, double *p2);
int fitness_match(double f1, double f2);

void free_redundancy(REDUNDANCY **r);
void new_redundancy(REDUNDANCY **r, double fitness, int number_parameters, double *parameters, double *velocity);

int fread_redundancy(FILE *file, REDUNDANCY **r, int number_parameters, int *particle);
int fwrite_redundancy(FILE *file, REDUNDANCY *r, int number_parameters, int particle);
#endif
