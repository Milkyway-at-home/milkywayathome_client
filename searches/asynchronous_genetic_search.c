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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "asynchronous_search.h"
#include "asynchronous_genetic_search.h"
#include "regression.h"
#include "newton_method.h"
#include "outliers.h"
#include "population.h"
#include "recombination.h"
#include "search_log.h"
#include "search_parameters.h"

#include "../evaluation/search_manager.h"
#include "../util/settings.h"
#include "../util/matrix.h"
#include "../util/io_util.h"

#define DSFMT_MEXP 19937
#define DSFMT_DO_NOT_USE_OLD_NAMES
#include "../../mersenne_twister/dSFMT.h"


ASYNCHRONOUS_SEARCH* get_asynchronous_genetic_search() {
	ASYNCHRONOUS_SEARCH *as = (ASYNCHRONOUS_SEARCH*)malloc(sizeof(ASYNCHRONOUS_SEARCH));
	as->search_qualifier = (char*)malloc(sizeof(char) * SEARCH_QUALIFIER_SIZE);
	strcpy(as->search_qualifier, "gs");
	as->create_search = create_genetic_search;
	as->read_search = read_genetic_search;
	as->checkpoint_search = checkpoint_genetic_search;
	as->generate_parameters = gs_generate_parameters;
	as->insert_parameters = gs_insert_parameters;
	return as;
}

int create_genetic_search(char* search_name, int number_arguments, char** arguments, int number_parameters, double *point, double *range, BOUNDS *bounds) {
	char search_directory[FILENAME_SIZE];
	GENETIC_SEARCH* gs;
	int i, population_size;

	gs = (GENETIC_SEARCH*)malloc(sizeof(GENETIC_SEARCH));

	gs->type = GENETIC_SIMPLEX;
	gs->number_parameters = number_parameters;
	gs->current_evaluation = 0;
	gs->mutation_rate = 0.3;
	gs->redundancy_rate = 0.3;
	gs->number_parents = 2;
	gs->ls_center = -1.5;
	gs->ls_outside = 1.5;
	gs->bounds = bounds;
	gs->no_redundancy = 0;

	population_size = 100;

	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-gs_parents")) {
			gs->number_parents = atoi(arguments[++i]);
		} else if (!strcmp(arguments[i], "-gs_mutation_rate")) {
			gs->mutation_rate = atof(arguments[++i]);
		} else if (!strcmp(arguments[i], "-gs_redundancy_rate")) {
			gs->redundancy_rate = atof(arguments[++i]);
		} else if (!strcmp(arguments[i], "-gs_ls_range")) {
			gs->ls_center = atof(arguments[++i]);
			gs->ls_outside = atof(arguments[++i]);
		} else if (!strcmp(arguments[i], "-gs_type")) {
			i++;
			if (!strcmp(arguments[i], "average")) {
				gs->type = GENETIC_AVERAGE;
			} else if (!strcmp(arguments[i], "simplex")) {
				gs->type = GENETIC_SIMPLEX;
			} else {
				sprintf(AS_MSG, "unknown search type: %s", arguments[i]);
				return AS_CREATE_FAIL;
			}
		} else if (!strcmp(arguments[i], "-gs_population_size")) {
			population_size = atoi(arguments[++i]);
		} else if (!strcmp(arguments[i], "-no_redundancy")) {
			gs->no_redundancy = 1;
		}
	}

	sprintf(search_directory, "%s/%s", get_working_directory(), search_name);
	printf("making directory: %s\n", search_directory);
	mkdir(search_directory, 0777);

	printf("population size: %d, number_parameters: %d\n", population_size, gs->number_parameters);
	new_population(population_size, gs->number_parameters, &(gs->population));
	initialize_redundancies(&(gs->redundancies));

	return checkpoint_genetic_search(search_name, gs);
}

int read_genetic_search(char *search_name, void** search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *file;
	GENETIC_SEARCH *gs;

	*search_data = (GENETIC_SEARCH*)malloc(sizeof(GENETIC_SEARCH));
	gs = (GENETIC_SEARCH*)(*search_data);

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	file = fopen(search_filename, "r");

	if (file == NULL) {
		fprintf(stderr, "Error reading genetic search: [%s]\n", search_name);
		return AS_READ_ERROR;
	}

	printf("reading genetic search\n");
	fscanf(file, "type: %d\n", &(gs->type));
	fscanf(file, "current_evaluation: %d\n", &(gs->current_evaluation));
	fscanf(file, "no_redundancy: %d, redundancy_rate: %lf\n", &(gs->no_redundancy), &(gs->redundancy_rate));
	fscanf(file, "mutation_rate: %lf\n", &(gs->mutation_rate));
	fscanf(file, "number_parents: %d, ls_center: %lf, ls_outside: %lf\n", &(gs->number_parents), &(gs->ls_center), &(gs->ls_outside));
	fread_bounds(file, &(gs->bounds));
	fclose(file);

	printf("read file\n");

	sprintf(population_filename, "%s/%s/redundancies", get_working_directory(), search_name);
	if (0 > read_redundancies(population_filename, &(gs->redundancies))) return AS_READ_ERROR;

	printf("read redundancies\n");

	sprintf(population_filename, "%s/%s/population", get_working_directory(), search_name);
	if (0 > read_population(population_filename, &(gs->population))) return AS_READ_ERROR;

	printf("read population\n");

	gs->number_parameters = gs->population->number_parameters;

	printf("init mersenne twister\n");
	dsfmt_gv_init_gen_rand((int)time(NULL));
	printf("initialized\n");

	return AS_READ_SUCCESS; 
}


int checkpoint_genetic_search(char* search_name, void* search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *file;
	GENETIC_SEARCH *gs = (GENETIC_SEARCH*)search_data;

	printf("search_name: %s\n", search_name);

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	file = fopen(search_filename, "w");

	if (file == NULL) {
		fprintf(stderr, "Error writing genetic search: [%s]\n", search_name);
		return AS_CP_ERROR;
	}
	
	fprintf(file, "type: %d\n", gs->type);
	fprintf(file, "current_evaluation: %d\n", gs->current_evaluation);
	fprintf(file, "no_redundancy: %d, redundancy_rate: %.15lf\n", gs->no_redundancy, gs->redundancy_rate);
	fprintf(file, "mutation_rate: %.15lf\n", gs->mutation_rate);
	fprintf(file, "number_parents: %d, ls_center: %.15lf, ls_outside: %.15lf\n", gs->number_parents, gs->ls_center, gs->ls_outside);
	fwrite_bounds(file, gs->bounds);
	fclose(file);

	sprintf(population_filename, "%s/%s/redundancies", get_working_directory(), search_name);
	if (0 > write_redundancies(population_filename, gs->redundancies)) return AS_CP_ERROR;

	sprintf(population_filename, "%s/%s/population", get_working_directory(), search_name);
	if (0 > write_population(population_filename, gs->population)) return AS_CP_ERROR;

	sprintf(AS_MSG, "evaluation: %d", gs->current_evaluation);

	printf("search_name: %s\n", search_name);

	return AS_CP_SUCCESS; 
}

int gs_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	POPULATION *p, *parents;
	GENETIC_SEARCH *gs;

	gs = (GENETIC_SEARCH*)(search_data);
	p = gs->population;

	printf("redundancy rate: %.25lf, mutation rate: %.25lf, random: %.25lf\n", gs->redundancy_rate, gs->mutation_rate, dsfmt_gv_genrand_close_open());
	printf("search parameters size: %d\n", sizeof(sp->parameters));

	if (gs->no_redundancy == 0 && gs->redundancies->redundancy_list != NULL && dsfmt_gv_genrand_close_open() < gs->redundancy_rate) {
		printf("generating redundancy\n");
		generate_redundancy(gs->redundancies, gs->number_parameters, sp->parameters, sp->metadata);
		strcat(sp->metadata, ", redundancy");
	} else if (p->size < p->max_size) {
		printf("generating random parameters\n");
		random_recombination(p->number_parameters, gs->bounds->min_bound, gs->bounds->max_bound, sp->parameters);
		sprintf(sp->metadata, "ev: %d, random", gs->current_evaluation);
	} else if (dsfmt_gv_genrand_close_open() < gs->mutation_rate) {
		printf("generating mutation\n");
		mutate(p->individuals[(int)(dsfmt_gv_genrand_close_open() * p->size)], gs->bounds->min_bound, gs->bounds->max_bound, gs->number_parameters, sp->parameters);
		sprintf(sp->metadata, "ev: %d, mutation", gs->current_evaluation);
	} else {
		printf("generating recombination\n");
		get_n_distinct(p, gs->number_parents, &parents);
		printf("got n distinct: %d\n", gs->number_parents);

		if (gs->type == GENETIC_AVERAGE) {
			printf("average recombination\n");
			average_recombination(parents->individuals, gs->number_parents, gs->number_parameters, sp->parameters);
			sprintf(sp->metadata, "ev: %d, average", gs->current_evaluation);
		} else if (gs->type == GENETIC_SIMPLEX) {
			printf("simplex recombination\n");
			free(sp->parameters);
			sp->parameters = (double*)malloc(sizeof(double) * gs->number_parameters);

			double point = simplex_recombination(parents->individuals, parents->fitness, gs->number_parents, gs->number_parameters, gs->ls_center, gs->ls_outside, sp->parameters);
			sprintf(sp->metadata, "ev: %d, simplex: %.20lf", gs->current_evaluation, point);
		} else {
			printf("UNKNOWN TYPE\n");
		}

		bound_parameters(sp->parameters, gs->bounds);

		free_population(parents);
		free(parents);
	}

	return AS_GEN_SUCCESS;
}


int gs_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	double best, average, median, worst, deviation;
	int verify_result;
	GENETIC_SEARCH *gs = (GENETIC_SEARCH*)search_data;
	gs->current_evaluation++;

	if (sp->fitness > -2.0) {
		sprintf(AS_MSG, "f: %.15lf, w: %.15lf, g: %.15lf, invalid", sp->fitness, gs->population->fitness[gs->population->size-1], gs->population->fitness[0]);
		return AS_INSERT_FITNESS_INVALID;
	} else if (!population_contains(gs->population, sp->fitness, sp->parameters)) {
		if (gs->population->size < gs->population->max_size || sp->fitness > gs->population->fitness[gs->population->size-1]) {
			if (gs->no_redundancy == 0) {
				verify_result = verify_with_insert(gs->redundancies, gs->number_parameters, sp->fitness, sp->parameters, sp->metadata, sp->hostid);
				if (verify_result == VERIFY_VALID) {
					int position = insert_sorted(gs->population, sp->parameters, sp->fitness);
					sprintf(AS_MSG, "f: %.15lf, w: %.15lf, g: %.15lf, inserted: %d", sp->fitness, gs->population->fitness[gs->population->size-1], gs->population->fitness[0], position);

					get_population_statistics(gs->population, &best, &average, &median, &worst, &deviation);
					log_printf(search_name, "%ld -- b: %.15lf, a: %.15lf, m: %.15lf, w: %.15lf, d: %.15lf, insert at: %ld\n", gs->current_evaluation, best, average, median, worst, deviation, position);
				} else {
					if (gs->population->size > 0) {
						sprintf(AS_MSG, "f: %.15lf, w: %.15lf, g: %.15lf, verifying", sp->fitness, gs->population->fitness[gs->population->size-1], gs->population->fitness[0]);
					} else {
						sprintf(AS_MSG, "f: %.15lf, w: %.15lf, g: %.15lf, verifying", sp->fitness, 0.0, 0.0);
					}
				}
			} else {
				int position = insert_sorted(gs->population, sp->parameters, sp->fitness);
				sprintf(AS_MSG, "f: %.15lf, w: %.15lf, g: %.15lf, inserted: %d", sp->fitness, gs->population->fitness[gs->population->size-1], gs->population->fitness[0], position);

				get_population_statistics(gs->population, &best, &average, &median, &worst, &deviation);
				log_printf(search_name, "%ld -- b: %.15lf, a: %.15lf, m: %.15lf, w: %.15lf, d: %.15lf, insert at: %ld\n", gs->current_evaluation, best, average, median, worst, deviation, position);
			}
		} else {
			verify_without_insert(gs->redundancies, gs->number_parameters, sp->fitness, sp->parameters, sp->metadata, sp->hostid);
			sprintf(AS_MSG, "f: %.15lf, w: %.15lf, g: %.15lf, not inserted", sp->fitness, gs->population->fitness[gs->population->size-1], gs->population->fitness[0]);
		}
	} else {
		sprintf(AS_MSG, "f: %.15lf, w: %.15lf, g: %.15lf, duplicate", sp->fitness, gs->population->fitness[gs->population->size-1], gs->population->fitness[0]);
	}

	return AS_INSERT_SUCCESS;
}

