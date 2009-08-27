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
#include <float.h>
#include <time.h>

#include "asynchronous_search.h"
#include "asynchronous_differential_evolution.h"
#include "bounds.h"
#include "outliers.h"
#include "population.h"
#include "recombination.h"
#include "redundancies.h"
#include "search_log.h"
#include "search_parameters.h"

#include "../evaluation/search_manager.h"
#include "../util/settings.h"
#include "../util/matrix.h"
#include "../util/io_util.h"

/********
 * 	Mersenne Twister Includes
 ********/
#define DSFMT_MEXP 19937
#define DSFMT_DO_NOT_USE_OLD_NAMES
#include "../../mersenne_twister/dSFMT.h"

ASYNCHRONOUS_SEARCH* get_asynchronous_differential_evolution() {
	ASYNCHRONOUS_SEARCH *as = (ASYNCHRONOUS_SEARCH*)malloc(sizeof(ASYNCHRONOUS_SEARCH));
	as->search_qualifier = (char*)malloc(sizeof(char) * SEARCH_QUALIFIER_SIZE);
	strcpy(as->search_qualifier, "de");
	as->create_search = create_differential_evolution;
	as->read_search = read_differential_evolution;
	as->generate_parameters = de_generate_parameters;
	as->insert_parameters = de_insert_parameters;
	return as;
}

int create_differential_evolution(char* search_name, int number_arguments, char** arguments, int number_parameters, double* point, double *range, BOUNDS* bounds) {
	char search_directory[FILENAME_SIZE];
	DIFFERENTIAL_EVOLUTION *de;
	int i;
	double redundancy_rate;

	de = (DIFFERENTIAL_EVOLUTION*)malloc(sizeof(DIFFERENTIAL_EVOLUTION));
	de->recombination_type = DE_RECOMBINATION_NONE;
	de->parent_type = DE_PARENT_BEST;
	de->recombination_pairs = 1;
	de->pair_weight = 0.5;
	de->crossover_rate = 0.5;
	redundancy_rate = 1;

	de->current_individual = 0;
	de->analyzed = 0;

	de->bounds = bounds;
	de->number_parameters = number_parameters;

	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-population_size")) de->population_size = atoi(arguments[++i]);
		else if (!strcmp(arguments[i], "-redundancy_rate")) redundancy_rate = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-pair_weight")) de->pair_weight = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-crossover_rate")) de->crossover_rate = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-recombination_pairs")) de->recombination_pairs = atoi(arguments[++i]);
		else if (!strcmp(arguments[i], "-recombination_type")) {
			i++;
			if (!strcmp(arguments[i], "none")) de->recombination_type = DE_RECOMBINATION_NONE;
			else if (!strcmp(arguments[i], "binomial")) de->recombination_type = DE_RECOMBINATION_BINOMIAL;
			else if (!strcmp(arguments[i], "exponential")) de->recombination_type = DE_RECOMBINATION_EXPONENTIAL;
			else {
				printf("unknown recombination type: %s\n", arguments[i]);
				exit(0);
			}
		} else if (!strcmp(arguments[i], "-parent_type")) {
			i++;
			if (!strcmp(arguments[i], "best")) de->parent_type = DE_PARENT_BEST;
			else if (!strcmp(arguments[i], "current")) de->parent_type = DE_PARENT_CURRENT;
			else if (!strcmp(arguments[i], "random")) de->parent_type = DE_PARENT_RANDOM;
			else {
				printf("unknown parent type: %s\n", arguments[i]);
				exit(0);
			}
		}
	}

	printf("initializing best individual\n");
	de->best_individual = (double*)malloc(sizeof(double) * de->number_parameters);
	de->best_individual_fitness = -DBL_MAX;
	de->best_individual_position = -1;

	printf("initializing population, size: %d, number_parameters: %d\n", de->population_size, de->number_parameters);
	new_population(de->population_size, de->number_parameters, &(de->population));
	for (i = 0; i < de->population_size; i++) de->population->fitness[i] = 0;

	printf("initializing redundancies, rate: %lf (1.0 == dynamic rate)\n", redundancy_rate);
	initialize_redundancies(de->population_size, redundancy_rate, &(de->redundancies));

	sprintf(search_directory, "%s/%s", get_working_directory(), search_name);
	printf("making directory: %s\n", search_directory);
	mkdir(search_directory, 0777);

	return write_differential_evolution(search_name, de);
}

int write_differential_evolution(char* search_name, void* search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	DIFFERENTIAL_EVOLUTION *de = (DIFFERENTIAL_EVOLUTION*)search_data;

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "w");
	if (search_file == NULL) return AS_CP_ERROR;

	fprintf(search_file, "population_size: %d, number_parameters: %d\n", de->population_size, de->number_parameters);
	fprintf(search_file, "analyzed: %d, current_individual: %d\n", de->analyzed, de->current_individual);
	fprintf(search_file, "parent_type: %d\n", de->parent_type);
	fprintf(search_file, "pair_weight: %.15lf, crossover_rate: %.15lf, recombination_pairs: %d, recombination_type: %d\n", de->pair_weight, de->crossover_rate, de->recombination_pairs, de->recombination_type);

	fprintf(search_file, "best_individual_position: %d, best_individual_fitness: %.20lf\n", de->best_individual_position, de->best_individual_fitness);
	fwrite_double_array(search_file, "best_individual", de->number_parameters, de->best_individual);

	fwrite_bounds(search_file, de->bounds);

	fclose(search_file);

	sprintf(population_filename, "%s/%s/redundancies", get_working_directory(), search_name);
	if (0 > write_redundancies(population_filename, de->redundancies)) return AS_CP_ERROR;

	sprintf(population_filename, "%s/%s/population", get_working_directory(), search_name);
	if (0 > write_population(population_filename, de->population)) return AS_CP_ERROR;

	return AS_CP_SUCCESS;
}

int read_differential_evolution(char* search_name, void** search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	DIFFERENTIAL_EVOLUTION *de;

	(*search_data) = (DIFFERENTIAL_EVOLUTION*)malloc(sizeof(DIFFERENTIAL_EVOLUTION));
	de = (DIFFERENTIAL_EVOLUTION*)(*search_data);

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "r");
	if (search_file == NULL) return AS_READ_ERROR;

	fscanf(search_file, "population_size: %d, number_parameters: %d\n", &(de->population_size), &(de->number_parameters));
	fscanf(search_file, "analyzed: %d, current_individual: %d\n", &(de->analyzed), &(de->current_individual));
	fscanf(search_file, "parent_type: %d\n", &(de->parent_type));
	fscanf(search_file, "pair_weight: %lf, crossover_rate: %lf, recombination_pairs: %d, recombination_type: %d\n", &(de->pair_weight), &(de->crossover_rate), &(de->recombination_pairs), &(de->recombination_type));

	fscanf(search_file, "best_individual_position: %d, best_individual_fitness: %lf\n", &(de->best_individual_position), &(de->best_individual_fitness));
	fread_double_array(search_file, "best_individual", &(de->best_individual));

	fread_bounds(search_file, &(de->bounds));

	fclose(search_file);

	sprintf(population_filename, "%s/%s/redundancies", get_working_directory(), search_name);
	if (0 > read_redundancies(population_filename, &(de->redundancies)) ) return AS_READ_ERROR;

	sprintf(population_filename, "%s/%s/population", get_working_directory(), search_name);
	if (0 > read_population(population_filename, &(de->population)) ) return AS_READ_ERROR;

	dsfmt_gv_init_gen_rand((int)time(NULL));

	return AS_READ_SUCCESS;
}

int de_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	DIFFERENTIAL_EVOLUTION *de = (DIFFERENTIAL_EVOLUTION*)search_data;

	if (generate_redundancy(de->redundancies, sp->number_parameters, sp->parameters, sp->metadata)) {
		sprintf(strchr(sp->metadata, 0), ", redundancy");
	} else if (de->population->size < de->population_size) {
		random_recombination(de->number_parameters, de->bounds->min_bound, de->bounds->max_bound, sp->parameters);
		sprintf(sp->metadata, "i: %d, random", de->current_individual);

		de->current_individual++;
		if (de->current_individual >= de->population_size) de->current_individual = 0;
	} else {
		double *parent, *pair_sum;
		parent = NULL;
		POPULATION *recombination_pairs;

		if (de->parent_type == DE_PARENT_BEST) {
			parent = de->best_individual;
			get_n_distinct_exclude(de->population, 2 * de->recombination_pairs, &recombination_pairs, de->best_individual_position);
		} else if (de->parent_type == DE_PARENT_CURRENT) {
			parent = de->population->individuals[de->current_individual];
			get_n_distinct_exclude(de->population, 2 * de->recombination_pairs, &recombination_pairs, de->current_individual);
		} else if (de->parent_type == DE_PARENT_RANDOM) {
			int position = dsfmt_gv_genrand_close_open() * de->population_size;
			parent = de->population->individuals[position];
			get_n_distinct_exclude(de->population, 2 * de->recombination_pairs, &recombination_pairs, position);
		}


		pair_sum = (double*)malloc(sizeof(double) * de->number_parameters);
		get_pair_sum(de->pair_weight, parent, recombination_pairs->individuals, de->recombination_pairs, recombination_pairs->number_parameters, pair_sum);

		if (de->recombination_type == DE_RECOMBINATION_NONE) {
			memcpy(sp->parameters, pair_sum, sizeof(double) * de->number_parameters);
		} else if (de->recombination_type == DE_RECOMBINATION_BINOMIAL) {
			binomial_recombination(de->crossover_rate, pair_sum, de->population->individuals[de->current_individual], de->number_parameters, sp->parameters);
		} else if (de->recombination_type == DE_RECOMBINATION_EXPONENTIAL) {
			exponential_recombination(de->crossover_rate, pair_sum, de->population->individuals[de->current_individual], de->number_parameters, sp->parameters);
		}

		sprintf(sp->metadata, "i: %d", de->current_individual);

		de->current_individual++;
		if (de->current_individual >= de->population_size) de->current_individual = 0;

		free(pair_sum);
		free_population(recombination_pairs);
		free(recombination_pairs);
	}

	return AS_GEN_SUCCESS;
}

int de_parse(DIFFERENTIAL_EVOLUTION *de, SEARCH_PARAMETERS *sp, int *position) {
	int i;
	char *metadata;

	metadata = (char*)malloc(sizeof(char) * METADATA_SIZE);
	memcpy(metadata, sp->metadata, sizeof(char) * METADATA_SIZE);

	if (isnan(sp->fitness)) return AS_INSERT_FITNESS_NAN;

	for (i = 0; i < de->number_parameters; i++) {
		if (isnan(sp->parameters[i])) return AS_INSERT_PARAMETERS_NAN;
		if (sp->parameters[i] < de->bounds->min_bound[i] || sp->parameters[i] > de->bounds->max_bound[i]) return AS_INSERT_OUT_OF_BOUNDS;
	}

	if (1 != sscanf(metadata, "i: %d", position)) return AS_INSERT_INVALID_METADATA;
	free(metadata);

	return 0;
}

int de_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	DIFFERENTIAL_EVOLUTION *de = (DIFFERENTIAL_EVOLUTION*)search_data;
	int result, position, verify_result;

	result = de_parse(de, sp, &position);
	if (result != 0) {
		sprintf(AS_MSG, "parse error");
		return result;
	}

	de->analyzed++;

	if (!population_contains(de->population, sp->fitness, sp->parameters)) {
		if (!individual_exists(de->population, position) || sp->fitness > de->population->fitness[position]) {
			verify_result = verify_with_insert(de->redundancies, sp->number_parameters, sp->fitness, sp->parameters, sp->metadata, sp->hostid);
			if (verify_result == VERIFY_VALID) {
			        double previous, best, average, median, worst, deviation;

				if (!individual_exists(de->population, position)) previous = 0;
				else previous = de->population->fitness[position];

				insert_individual_info(de->population, position, sp->parameters, sp->fitness, sp->host_os, sp->app_version);

				get_population_statistics(de->population, &best, &average, &median, &worst, &deviation);
				if (sp->fitness > de->best_individual_fitness) {
					de->best_individual_fitness = sp->fitness;
					de->best_individual_position = position;
					memcpy(de->best_individual, sp->parameters, sizeof(double) * de->number_parameters);
					sprintf(AS_MSG, "i[%d]: %.15lf, l: %.15lf, g: %.15lf, global best", position, sp->fitness, previous, de->best_individual_fitness);
					log_printf(search_name, "%ld -- b: %.15lf, a: %.15lf, m: %.15lf, w: %.15lf, d: %.15lf, global best\n", de->analyzed, best, average, median, worst, deviation);
				} else {
					sprintf(AS_MSG, "i[%d]: %.15lf, l: %.15lf, g: %.15lf, local best", position, sp->fitness, previous, de->best_individual_fitness);
					log_printf(search_name, "%ld -- b: %.15lf, a: %.15lf, m: %.15lf, w: %.15lf, d: %.15lf, local best\n", de->analyzed, best, average, median, worst, deviation);
				}
				write_differential_evolution(search_name, de);
			} else {
				if (!individual_exists(de->population, position)) { 
					if (de->population->size == 0) {
						sprintf(AS_MSG, "i[%d]: %.15lf, l: %.15lf, g: %.15lf, verifying", position, sp->fitness, 0.0, 0.0);
					} else {
						sprintf(AS_MSG, "i[%d]: %.15lf, l: %.15lf, g: %.15lf, verifying", position, sp->fitness, 0.0, de->best_individual_fitness);
					}
				} else {
					sprintf(AS_MSG, "i[%d]: %.15lf, l: %.15lf, g: %.15lf, verifying", position, sp->fitness, de->population->fitness[position], de->best_individual_fitness);
				}
			}
		} else {
			verify_without_insert(de->redundancies, de->number_parameters, sp->fitness, sp->parameters, sp->metadata, sp->hostid);
			sprintf(AS_MSG, "i[%d]: %.15lf, l: %.15lf, g: %.15lf, low fitness", position, sp->fitness, de->population->fitness[position], de->best_individual_fitness);
		}
	} else {
		sprintf(AS_MSG, "i[%d]: %.15lf, l: %.15lf, g: %.15lf, duplicate", position, sp->fitness, de->population->fitness[position], de->best_individual_fitness);
	}
	return AS_INSERT_SUCCESS;
}
