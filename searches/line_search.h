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

#ifndef FGDO_LINE_SEARCH_H
#define FGDO_LINE_SEARCH_H

extern const char *LS_STR[];
#define LS_SUCCESS	0
#define LS_LOOP1_NAN	1
#define LS_LOOP2_NAN	2
#define LS_LOOP3_NAN	3
#define LS_LOOP1_TOL	4
#define LS_LOOP1_MAX	5
#define LS_LOOP2_MAX	6


int line_search(double* point, double initial_fitness, double* direction, int number_parameters, double* minimum, double* fitness, int *evaluations_done);

void randomized_line_search(int number_parameters, double *point, double *step, int ls_evaluations, int ls_iterations, double **new_point, double *fitness);

#endif
