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

#ifndef GEM_GENETIC_SEARCH_H
#define GEM_GENETIC_SEARCH_H

#include "population.h"

typedef struct genetic_search {
	int recombination_type;

	double **parents;
	double *parent_fitness;
	int current_generated;

	double mutation_rate;
	int number_parents, number_children;
	double simplex_l1, simplex_l2;
	double crossover_rate, crossover_scale;
} GENETIC_SEARCH;

void parse_genetic_search(char search_parameters[512], GENETIC_SEARCH *gs, int *population_size, int *max_evaluations);

void genetic_search__insert_individual(POPULATION *population, double* parameters, double fitness, char *metadata);

void genetic_search__get_individual(POPULATION *population, double** parameters, char **metadata);

void start_genetic_search(char *search_path, char *search_parameters, double *min_parameters, double *max_parameters, int number_parameters, POPULATION **population);

#endif
