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

#ifndef GEO_DIFFERENTIAL_EVOLUTION_H
#define GEO_DIFFERENTIAL_EVOLUTION_H

#include "population.h"

typedef struct differential_evolution {
	int current_best;

	int parent_type;
	double parent_scale;
	int number_pairs, pair_type;
	double pair_scale;
	int recombination_type;
	double crossover_rate;

	int next_generated;
} DIFFERENTIAL_EVOLUTION;

void parse_differential_evolution(char search_parameters[512], DIFFERENTIAL_EVOLUTION *de, int *population_size, int *max_evaluations);

void start_differential_evolution(char search_path[512], char search_parameters[512], double* min_parameters, double* max_parameters, int number_parameters, POPULATION **population);

void differential_evolution__insert_individual(POPULATION *population, double* parameters, double fitness, char *metadata);

void differential_evolution__get_individual(POPULATION *population, double** parameters, char **metadata);

#endif
