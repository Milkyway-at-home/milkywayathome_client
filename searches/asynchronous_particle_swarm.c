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
#include "asynchronous_particle_swarm.h"
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

ASYNCHRONOUS_SEARCH* get_asynchronous_particle_swarm() {
	ASYNCHRONOUS_SEARCH *as = (ASYNCHRONOUS_SEARCH*)malloc(sizeof(ASYNCHRONOUS_SEARCH));
	as->search_qualifier = (char*)malloc(sizeof(char) * SEARCH_QUALIFIER_SIZE);
	strcpy(as->search_qualifier, "ps");
	as->create_search = create_particle_swarm;
	as->read_search = read_particle_swarm;
	as->checkpoint_search = checkpoint_particle_swarm;
	as->generate_parameters = pso_generate_parameters;
	as->insert_parameters = pso_insert_parameters;
	return as;
}

int create_particle_swarm(char* search_name, int number_arguments, char** arguments, int number_parameters, double* point, double *range, BOUNDS* bounds) {
	char search_directory[FILENAME_SIZE];
	PARTICLE_SWARM_OPTIMIZATION *pso;
	int i;

	pso = (PARTICLE_SWARM_OPTIMIZATION*)malloc(sizeof(PARTICLE_SWARM_OPTIMIZATION));
	pso->w = 0.5;
	pso->c0 = 1.0;
	pso->c1 = 2.0;
	pso->c2 = 2.0;
	pso->size = 50;

	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-pso_size")) pso->size = atoi(arguments[++i]);
		else if (!strcmp(arguments[i], "-w")) pso->w = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-c0")) pso->c0 = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-c1")) pso->c1 = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-c2")) pso->c2 = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-redundancy_rate")) pso->redundancy_rate = atof(arguments[++i]);
	}

	sprintf(search_directory, "%s/%s", get_working_directory(), search_name);
	printf("making directory: %s\n", search_directory);
	mkdir(search_directory, 0777);

	pso->current_particle = 0;
	pso->analyzed = 0;

	pso->number_parameters = number_parameters;
	pso->bounds = bounds;

	pso->global_best_fitness = -DBL_MAX;
	pso->global_best = (double*)malloc(sizeof(double) * pso->number_parameters);
	for (i = 0; i < pso->size; i++) pso->global_best[i] = 0.0;

	new_population(pso->size, pso->number_parameters, &(pso->local_best));
	for (i = 0; i < pso->size; i++) pso->local_best->fitness[i] = 0;
	new_population(pso->size, pso->number_parameters, &(pso->particles));
	new_population(pso->size, pso->number_parameters, &(pso->velocities));

	initialize_redundancies(&(pso->redundancies));

	return checkpoint_particle_swarm(search_name, pso);
}

int write_particle_swarm(char* search_name, void* search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	PARTICLE_SWARM_OPTIMIZATION *pso = (PARTICLE_SWARM_OPTIMIZATION*)search_data;

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "w");
	if (search_file == NULL) return AS_CP_ERROR;

	fprintf(search_file, "current_particle: %d, size: %d\n", pso->current_particle, pso->size);
	fprintf(search_file, "w: %lf, c0: %lf, c1: %lf, c2: %lf\n", pso->w, pso->c0, pso->c1, pso->c2);
	fprintf(search_file, "redundancy_rate: %lf\n", pso->redundancy_rate);
	fprintf(search_file, "analyzed: %ld\n", pso->analyzed);
	fprintf(search_file, "number_parameters: %d\n", pso->number_parameters); 

	fwrite_bounds(search_file, pso->bounds);

	fprintf(search_file, "global_best_fitness: %.20lf\n", pso->global_best_fitness);
	fwrite_double_array(search_file, "global_best", pso->number_parameters, pso->global_best);
	fclose(search_file);

	sprintf(population_filename, "%s/%s/redundancies", get_working_directory(), search_name);
	if (0 > write_redundancies(population_filename, pso->redundancies)) return AS_CP_ERROR;

	sprintf(population_filename, "%s/%s/particles", get_working_directory(), search_name);
	if (0 > write_population(population_filename, pso->particles)) return AS_CP_ERROR;

	sprintf(population_filename, "%s/%s/velocities", get_working_directory(), search_name);
	if (0 > write_population(population_filename, pso->velocities)) return AS_CP_ERROR;

	sprintf(population_filename, "%s/%s/local_best", get_working_directory(), search_name);
	if (0 > write_population(population_filename, pso->local_best)) return AS_CP_ERROR;

	return AS_CP_SUCCESS;
}

int read_particle_swarm(char* search_name, void** search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	PARTICLE_SWARM_OPTIMIZATION *pso;

	(*search_data) = (PARTICLE_SWARM_OPTIMIZATION*)malloc(sizeof(PARTICLE_SWARM_OPTIMIZATION));
	pso = (PARTICLE_SWARM_OPTIMIZATION*)(*search_data);

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "r");
	if (search_file == NULL) return AS_READ_ERROR;

	fscanf(search_file, "current_particle: %d, size: %d\n", &(pso->current_particle), &(pso->size));
	fscanf(search_file, "w: %lf, c0: %lf, c1: %lf, c2: %lf\n", &(pso->w), &(pso->c0), &(pso->c1), &(pso->c2));
	fscanf(search_file, "redundancy_rate: %lf\n", &(pso->redundancy_rate));
	fscanf(search_file, "analyzed: %ld\n", &(pso->analyzed));
	fscanf(search_file, "number_parameters: %d\n", &(pso->number_parameters));

	fread_bounds(search_file, &(pso->bounds));

	fscanf(search_file, "global_best_fitness: %lf\n", &(pso->global_best_fitness));
	fread_double_array(search_file, "global_best", &(pso->global_best));
	fclose(search_file);

	sprintf(population_filename, "%s/%s/redundancies", get_working_directory(), search_name);
	if (0 > read_redundancies(population_filename, &(pso->redundancies)) ) return AS_READ_ERROR;

	sprintf(population_filename, "%s/%s/particles", get_working_directory(), search_name);
	if (0 > read_population(population_filename, &(pso->particles)) ) return AS_READ_ERROR;
	
	sprintf(population_filename, "%s/%s/velocities", get_working_directory(), search_name);
	if (0 > read_population(population_filename, &(pso->velocities)) ) return AS_READ_ERROR;

	sprintf(population_filename, "%s/%s/local_best", get_working_directory(), search_name);
	if (0 > read_population(population_filename, &(pso->local_best)) ) return AS_READ_ERROR;

	dsfmt_gv_init_gen_rand((int)time(NULL));

	return AS_READ_SUCCESS;
}


int checkpoint_particle_swarm(char* search_name, void* search_data) {
	PARTICLE_SWARM_OPTIMIZATION *pso = (PARTICLE_SWARM_OPTIMIZATION*)search_data;

	sprintf(AS_MSG, "current_particle: %d", pso->current_particle);
	return write_particle_swarm(search_name, search_data);
}

int pso_generate_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	PARTICLE_SWARM_OPTIMIZATION *pso = (PARTICLE_SWARM_OPTIMIZATION*)search_data;
	int i;

	if (pso->redundancies->redundancy_list != NULL && dsfmt_gv_genrand_close_open() < pso->redundancy_rate) {
		generate_redundancy(pso->redundancies, sp->number_parameters, sp->parameters, sp->metadata);
		sprintf(strchr(sp->metadata, 0), ", redundancy");
	} else if (!individual_exists(pso->local_best, pso->current_particle)) {
		/********
			*	This particle hasn't yet been created.
		 ********/
		random_recombination(pso->number_parameters, pso->bounds->min_bound, pso->bounds->max_bound, sp->parameters);
		sprintf(sp->metadata, "p: %d, v:", pso->current_particle);
		for (i = 0; i < pso->number_parameters; i++) sprintf(strchr(sp->metadata, 0), " %lf", 0.0);
		pso->current_particle++;
		if (pso->current_particle >= pso->size) pso->current_particle = 0;
	} else {
		double *local_best, *velocity, *particle;
		sprintf(sp->metadata, "p: %d, v:", pso->current_particle);

		local_best = pso->local_best->individuals[pso->current_particle];
		velocity = pso->velocities->individuals[pso->current_particle];
		particle = pso->particles->individuals[pso->current_particle];

		for (i = 0; i < pso->number_parameters; i++) {
			velocity[i] = (pso->w * velocity[i]) + pso->c0 * ((pso->c1 * dsfmt_gv_genrand_close_open() * (local_best[i] - particle[i])) + (pso->c2 * dsfmt_gv_genrand_close_open() * (pso->global_best[i] - particle[i])));
		}
		bound_velocity(particle, velocity, pso->bounds);

		for (i = 0; i < pso->number_parameters; i++) {
			particle[i] = particle[i] + velocity[i];
			if (isnan(particle[i])) return AS_GEN_FAIL;
		}
		bound_parameters(particle, pso->bounds);
		
		pso->current_particle++;
		if (pso->current_particle >= pso->size) pso->current_particle = 0;

		for (i = 0; i < pso->number_parameters; i++) {
			sp->parameters[i] = particle[i];
			sprintf(strchr(sp->metadata, 0), " %.20lf", velocity[i]);
		}
	}

	return AS_GEN_SUCCESS;
}

int parse(PARTICLE_SWARM_OPTIMIZATION *pso, SEARCH_PARAMETERS *sp, int *particle, double *velocity) {
	int i;
	char *current_token;
	char *metadata;

	metadata = (char*)malloc(sizeof(char) * METADATA_SIZE);
	memcpy(metadata, sp->metadata, sizeof(char) * METADATA_SIZE);

	if (isnan(sp->fitness)) return AS_INSERT_FITNESS_NAN;
//	if (sp->fitness > -2.0) return AS_INSERT_FITNESS_INVALID;

	for (i = 0; i < pso->number_parameters; i++) {
		if (isnan(sp->parameters[i])) return AS_INSERT_PARAMETERS_NAN;
		if (sp->parameters[i] < pso->bounds->min_bound[i] || sp->parameters[i] > pso->bounds->max_bound[i]) return AS_INSERT_OUT_OF_BOUNDS;
	}

	if (1 != sscanf(metadata, "p: %d, v:", particle)) return AS_INSERT_INVALID_METADATA;
	current_token = strtok(strchr(metadata, 'v'), ", :");
	for (i = 0; i < pso->number_parameters; i++) {
		if (current_token == NULL) return AS_INSERT_INVALID_METADATA;
		velocity[i] = atof(current_token);
		current_token = strtok(NULL, " ");
	}
	free(metadata);

	return 0;
}

void insert_particle(char* search_name, PARTICLE_SWARM_OPTIMIZATION *pso, int particle, double *velocity, SEARCH_PARAMETERS *sp) {
	double previous;
	double best, average, median, worst, deviation;

	if (!individual_exists(pso->local_best, particle)) previous = 0;
	else previous = pso->local_best->fitness[particle];
	/********
		*	Actually insert the particle.
	 ********/
	insert_individual_info(pso->particles, particle, sp->parameters, sp->fitness, sp->host_os, sp->app_version);
	insert_individual_info(pso->velocities, particle, velocity, sp->fitness, sp->host_os, sp->app_version);
	insert_individual_info(pso->local_best, particle, sp->parameters, sp->fitness, sp->host_os, sp->app_version);

	get_population_statistics(pso->local_best, &best, &average, &median, &worst, &deviation);
	if (sp->fitness > pso->global_best_fitness) {
		pso->global_best_fitness = sp->fitness;
		memcpy(pso->global_best, sp->parameters, sizeof(double) * pso->number_parameters);
		sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, global best", particle, sp->fitness, previous, pso->global_best_fitness);
		log_printf(search_name, "%ld -- b: %.15lf, a: %.15lf, m: %.15lf, w: %.15lf, d: %.15lf, global best\n", pso->analyzed, best, average, median, worst, deviation);
	} else {
		sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, local best", particle, sp->fitness, previous, pso->global_best_fitness);
		log_printf(search_name, "%ld -- b: %.15lf, a: %.15lf, m: %.15lf, w: %.15lf, d: %.15lf, local best\n", pso->analyzed, best, average, median, worst, deviation);
	}
	write_particle_swarm(search_name, pso);
}

int pso_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	PARTICLE_SWARM_OPTIMIZATION *pso = (PARTICLE_SWARM_OPTIMIZATION*)search_data;
	int result, particle, verify_result;
	double *velocity = (double*)malloc(sizeof(double) * pso->number_parameters);

	result = parse(pso, sp, &particle, velocity);
	if (result != 0) {
		free(velocity);
		sprintf(AS_MSG, "parse error");
		return result;
	}

	pso->analyzed++;

	if (!population_contains(pso->local_best, sp->fitness, sp->parameters)) {
		if (!individual_exists(pso->local_best, particle) || sp->fitness > pso->local_best->fitness[particle]) {
			verify_result = verify_with_insert(pso->redundancies, sp->number_parameters, sp->fitness, sp->parameters, sp->metadata, sp->hostid);
			if (verify_result == VERIFY_VALID) {
				insert_particle(search_name, pso, particle, velocity, sp);
			} else {
				sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, verifying", particle, sp->fitness, pso->local_best->fitness[particle], pso->global_best_fitness);
			}
		} else {
			sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, low fitness", particle, sp->fitness, pso->local_best->fitness[particle], pso->global_best_fitness);
		}
	} else {
		sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, duplicate", particle, sp->fitness, pso->local_best->fitness[particle], pso->global_best_fitness);
	}
	return AS_INSERT_SUCCESS;
}
