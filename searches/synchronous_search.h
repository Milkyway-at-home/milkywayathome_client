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

#ifndef GEO_SYNCHRONOUS_SEARCH_H
#define GEO_SYNCHRONOUS_SEARCH_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../searches/population.h"

void synchronous_search(char *search_path, char *search_parameters, double *min_parameters, double *max_parameters, int number_parameters, 
				void (*init_search)(char*, char*, double*, double*, int, POPULATION**));
void synchronous_parallel_search(char *search_path, char *search_parameters, double *min_parameters, double *max_parameters, int number_parameters, 
				void (*init_search)(char*, char*, double*, double*, int, POPULATION**), int number_workers);

#endif
