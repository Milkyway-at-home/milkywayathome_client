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

#ifndef FGDO_ASYNCHRONOUS_GENETIC_SEARCH_H
#define FGDO_ASYNCHRONOUS_GENETIC_SEARCH_H

#include <stdio.h>

#include "asynchronous_search.h"
#include "bounds.h"
#include "population.h"
#include "redundancies.h"

#define GENETIC_AVERAGE 1
#define GENETIC_SIMPLEX 2

typedef struct genetic_search {
	int current_evaluation;
	int number_parameters;

	double mutation_rate;

	int type;
	int number_parents;
	double ls_center, ls_outside;

	BOUNDS *bounds;
	POPULATION *population;
	REDUNDANCIES *redundancies;
} GENETIC_SEARCH;

ASYNCHRONOUS_SEARCH* get_asynchronous_genetic_search();

int create_genetic_search(char* search_name, int number_arguments, char** arguments, int number_parameters, double *point, double *range, BOUNDS *bounds);
int read_genetic_search(char* search_name, void** search_data);
int write_genetic_search(char* search_name, void* search_data);
int gs_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);
int gs_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS *sp);
#endif
