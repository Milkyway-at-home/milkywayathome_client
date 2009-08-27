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

#ifndef FGDO_ASYNCHRONOUS_DE_H
#define FGDO_ASYNCHRONOUS_DE_H

#include <stdio.h>

#include "asynchronous_search.h"
#include "bounds.h"
#include "population.h"
#include "redundancies.h"

#define DE_PARENT_BEST 0
#define DE_PARENT_CURRENT 1
#define DE_PARENT_RANDOM 2

#define DE_RECOMBINATION_NONE 0
#define DE_RECOMBINATION_BINOMIAL 1
#define DE_RECOMBINATION_EXPONENTIAL 2

typedef struct differential_evolution {
	double crossover_rate, pair_weight;
	int recombination_type, recombination_pairs;
	int parent_type;

	int current_individual, analyzed;
	int number_parameters;

	BOUNDS *bounds;

	int best_individual_position;
	double best_individual_fitness, *best_individual;

	int population_size;
	POPULATION *population;

	REDUNDANCIES *redundancies;
} DIFFERENTIAL_EVOLUTION;

ASYNCHRONOUS_SEARCH* get_asynchronous_differential_evolution();

int create_differential_evolution(char *search_name, int number_arguments, char** arguments, int number_parameters, double *point, double *range, BOUNDS *bounds);
int write_differential_evolution(char *search_name, void *search_data);
int read_differential_evolution(char *search_name, void **search_data);
int de_insert_parameters(char *search_name, void *search_data, SEARCH_PARAMETERS *sp);
int de_generate_parameters(char *search_name, void *search_data, SEARCH_PARAMETERS *sp);
#endif
