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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "genetic_search.h"
#include "recombination.h"
#include "population.h"
#include "../settings.h"

#define GS_AVERAGE 0
#define GS_DOUBLE_SHOT 1
#define GS_SIMPLEX 2
#define GS_BINOMIAL 3
#define GS_EXPONENTIAL 4

char* get_population_file(char search_path[512]) {
	char *population_file = (char*)malloc(sizeof(char) * 2048);

	strcat(population_file, search_path);
	strcat(population_file, "/population");

	return population_file;
}

void invalid_genetic_search() {
	fprintf(stderr, "Genetic Search was illdefined:\n");
	fprintf(stderr, "Search naming conventions are:\n");
	fprintf(stderr, "  Genetic Search:\n");
	fprintf(stderr, "    gs/average/<mutation_rate>/<population_size>/<max_evaluations>\n");
	fprintf(stderr, "    gs/doubleshot/<mutation_rate>/<population_size>/<max_evaluations>\n");
	fprintf(stderr, "    gs/simplex(parents=<number_parents>,children=<number_children>,l1=<l1>,l2=<l2>)/<mutation_rate>/<population_size>/<max_evaluations>\n");
	fprintf(stderr, "    gs/binomial(crossover_rate=<crossover_rate>,crossover_scale=<crossover_scale|random>])/<mutation_rate>/<population_size>/<max_evaluations>\n");
	fprintf(stderr, "    gs/exponential(crossover_rate=<crossover_rate>,crossover_scale=<crossover_scale|random>])/<mutation_rate>/<population_size>/<max_evaluations>\n");
}

void parse_genetic_search(char search_parameters[512], GENETIC_SEARCH *gs, int *population_size, int *max_evaluations) {
	gs->number_parents = 2;
	gs->mutation_rate = -1;
	gs->recombination_type = -1;
	gs->number_children = -1;
	gs->simplex_l1 = -1;
	gs->simplex_l2 = -1;
	gs->crossover_rate = -1;
	gs->crossover_scale = -1;
	gs->current_generated = 0;

	if (search_parameters[3] == 'a' || search_parameters[3] == 'A') {
		gs->recombination_type = GS_AVERAGE;
		if (sscanf(search_parameters, "gs/average/%lf/%d/%d", &gs->mutation_rate, population_size, max_evaluations) != 3) {
			invalid_genetic_search();
			return;
		}

	} else if (search_parameters[3] == 'd' || search_parameters[3] == 'D') {
		gs->recombination_type = GS_DOUBLE_SHOT;
		gs->number_children = 3;
		if (sscanf(search_parameters, "gs/doubleshot/%lf/%d/%d", &gs->mutation_rate, population_size, max_evaluations) != 3) {
			invalid_genetic_search();
			return;
		}

	} else if (search_parameters[3] == 's' || search_parameters[3] == 'S') {
		gs->recombination_type = GS_SIMPLEX;
		if (sscanf(search_parameters, "gs/simplex(parents=%d,children=%d,l1=%lf,l2=%lf)/%lf/%d/%d", &gs->number_parents, &gs->number_children, &gs->simplex_l1, &gs->simplex_l2, &gs->mutation_rate, population_size, max_evaluations) != 7) {
			invalid_genetic_search();
			return;
		}

	} else if (search_parameters[3] == 'b' || search_parameters[3] == 'B') {
		gs->recombination_type = GS_BINOMIAL;
		if (sscanf(search_parameters, "gs/binomial(crossover_rate=%lf,crossover_scale=%lf)/%lf/%d/%d", &gs->crossover_rate, &gs->crossover_scale, &gs->mutation_rate, population_size, max_evaluations) == 5) {
		} else if (sscanf(search_parameters, "gs/binomial(crossover_rate=%lf,crossover_scale=random)/%lf/%d/%d", &gs->crossover_rate, &gs->mutation_rate, population_size, max_evaluations) == 4) {
		} else {
			invalid_genetic_search();
			return;
		}

	} else if (search_parameters[3] == 'e' || search_parameters[3] == 'E') {
		gs->recombination_type = GS_EXPONENTIAL;
		if (sscanf(search_parameters, "gs/exponential(crossover_rate=%lf,crossover_scale=%lf)/%lf/%d/%d", &gs->crossover_rate, &gs->crossover_scale, &gs->mutation_rate, population_size, max_evaluations) == 5) {
		} else if (sscanf(search_parameters, "gs/exponential(crossover_rate=%lf,crossover_scale=random)/%lf/%d/%d", &gs->crossover_rate, &gs->mutation_rate, population_size, max_evaluations) == 4) {
		} else {
			invalid_genetic_search();
			return;
		}
	} else {
		(*max_evaluations) = 0;
		(*population_size) = 0;
		invalid_genetic_search();
		return;
	}
}

void genetic_search__insert_individual(POPULATION *population, double* parameters, double fitness, char *metadata) {
	insert_sorted(population, parameters, fitness);
}

void genetic_search__get_individual(POPULATION *population, double **parameters, char **metadata) {
	int target, i;
	GENETIC_SEARCH *gs;

	gs = (GENETIC_SEARCH*) population->parameters;

	(*metadata) = (char*)malloc(sizeof(char) * METADATA_SIZE);
	if (population->current_size <population->max_size) {
		(*parameters) = random_recombination(population->min_parameters, population->max_parameters, population->number_parameters);

	} else if (drand48() < gs->mutation_rate) {
		(*parameters) = (double*)malloc(sizeof(double));
		target = (int)(drand48() * population->current_size);
		(*parameters) = mutate(population->individuals[target], population->min_parameters, population->max_parameters, population->number_parameters);

	} else {
		gs = (GENETIC_SEARCH*)population->parameters;
		(*parameters) = (double*)malloc(sizeof(double));

		if (gs->recombination_type == GS_SIMPLEX || gs->recombination_type == GS_DOUBLE_SHOT) {
			if (gs->current_generated == 0) {
				get_n_distinct(population, gs->number_parents, &gs->parents, &gs->parent_fitness);
			}

			if (gs->recombination_type == GS_SIMPLEX) {
				(*parameters) = simplex_recombination(gs->parents, gs->parent_fitness, gs->number_parents, population->number_parameters, gs->simplex_l1, gs->simplex_l2);
				bound_parameters(population, (*parameters));
			} else if (gs->recombination_type == GS_DOUBLE_SHOT) {
				if (gs->current_generated == 0) {
					(*parameters) = lower_recombination(gs->parents, gs->number_parents, population->number_parameters);
				} else if (gs->current_generated == 1) {
					(*parameters) = average_recombination(gs->parents, gs->number_parents, population->number_parameters);
				} else if (gs->current_generated == 2) {
					(*parameters) = higher_recombination(gs->parents, gs->number_parents, population->number_parameters);
				}
				bound_parameters(population, (*parameters));
			}

			gs->current_generated++;
			if (gs->current_generated == gs->number_children) {
				gs->current_generated = 0;
				for (i = 0; i < gs->number_parents; i++) free(gs->parents[i]);
				free(gs->parents);
				free(gs->parent_fitness);
			}
		} else {
			get_n_distinct(population, gs->number_parents, &gs->parents, &gs->parent_fitness);

			if (gs->recombination_type == GS_AVERAGE) {
				(*parameters) = average_recombination(gs->parents, gs->number_parents, population->number_parameters);
			} else if (gs->recombination_type == GS_BINOMIAL) {
				(*parameters) = binomial_recombination(gs->parents, gs->number_parents, population->number_parameters, gs->crossover_rate, gs->crossover_scale);
			} else if (gs->recombination_type == GS_EXPONENTIAL) {
				(*parameters) = exponential_recombination(gs->parents, gs->number_parents, population->number_parameters, gs->crossover_rate, gs->crossover_scale);
			}
			bound_parameters(population, (*parameters));

			for (i = 0; i < gs->number_parents; i++) free(gs->parents[i]);
			free(gs->parents);
			free(gs->parent_fitness);
		}
	}
}


void start_genetic_search(char search_path[512], char search_parameters[512], double *min_parameters, double *max_parameters, int number_parameters, POPULATION **population) {
	GENETIC_SEARCH *gs;
        int population_size;
        int max_evaluations;
        
	gs = (GENETIC_SEARCH*)malloc(sizeof(GENETIC_SEARCH));
        parse_genetic_search(search_parameters, gs, &population_size, &max_evaluations);
        (*population) = new_population(search_path, search_parameters, min_parameters, max_parameters, number_parameters, population_size, max_evaluations);

        (*population)->parameters = gs;
	(*population)->get_individual = genetic_search__get_individual;
	(*population)->insert_individual = genetic_search__insert_individual;
}
