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

#ifndef FGDO_OUTLIERS_H
#define FGDO_OUTLIERS_H

#include "population.h"

double distance_error_from(POPULATION *p, double fitness, double *parameters);

void get_distance_errors(POPULATION *p, double **errors);
void get_error_stats(double *errors, int error_size, double *min_error, double *max_error, double *median_error, double *average_error);
void population_error_stats(POPULATION *p, double *min_error, double *max_error, double *median_error, double *average_error);

double get_error2(POPULATION *p, int position, double fitness);
double get_distance2(POPULATION *p, int position, double *parameters);
double distance_error2(POPULATION *p, int position, double fitness, double *parameters);

int remove_outliers(POPULATION *p, double range);
int remove_outliers_incremental(POPULATION *p, double range);
int remove_outliers_sorted(POPULATION *p, double range);
#endif
