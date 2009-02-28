#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <float.h>

#include "asynchronous_search.h"
#include "asynchronous_particle_swarm.h"
#include "bounds.h"
#include "outliers.h"
#include "population.h"
#include "recombination.h"
#include "redundancy.h"
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
	int i, remove_outliers, size;
	double w, c1, c2;

	remove_outliers = 1;
	size = 400;
	w = 0.95;
	c1 = 2.0;
	c2 = 2.0;
	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-pso_size")) size = atoi(arguments[++i]);
		else if (!strcmp(arguments[i], "-w")) w = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-c1")) c1 = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-c2")) c2 = atof(arguments[++i]);
		else if (!strcmp(arguments[i], "-remove_outliers")) remove_outliers = 1;
	}

	sprintf(search_directory, "%s/%s", get_working_directory(), search_name);
	printf("making directory: %s\n", search_directory);
	mkdir(search_directory, 0777);

	pso = (PARTICLE_SWARM_OPTIMIZATION*)malloc(sizeof(PARTICLE_SWARM_OPTIMIZATION));
	pso->remove_outliers = remove_outliers;
	pso->current_particle = 0;
	pso->size = size;
	pso->number_parameters = number_parameters;
	pso->w = w;
	pso->c1 = c1;
	pso->c2 = c2;
	pso->analyzed = 0;

	pso->global_best_fitness = -DBL_MAX;
	pso->global_best = (double*)malloc(sizeof(double) * pso->number_parameters);
	for (i = 0; i < pso->size; i++) pso->global_best[i] = 0.0;

	pso->bounds = bounds;

	new_population(pso->size, pso->number_parameters, &(pso->local_best));
	for (i = 0; i < pso->size; i++) pso->local_best->fitness[i] = 0;
	new_population(pso->size, pso->number_parameters, &(pso->particles));
	new_population(pso->size, pso->number_parameters, &(pso->velocities));

	pso->redundancies = (REDUNDANCY**)malloc(sizeof(REDUNDANCY*) * pso->size);
	for (i = 0; i < pso->size; i++) pso->redundancies[i] = NULL;
	pso->current_redundancy = NULL;

	return checkpoint_particle_swarm(search_name, pso);
}

int write_particle_swarm(char* search_name, void* search_data) {
	char search_filename[FILENAME_SIZE], population_filename[FILENAME_SIZE];
	FILE *search_file;
	int i;
	PARTICLE_SWARM_OPTIMIZATION *pso = (PARTICLE_SWARM_OPTIMIZATION*)search_data;

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "w+");
	if (search_file == NULL) return AS_CP_ERROR;

	fprintf(search_file, "current_particle: %d, size: %d\n", pso->current_particle, pso->size);
	fprintf(search_file, "w: %lf, c1: %lf, c2: %lf\n", pso->w, pso->c1, pso->c2);
	fprintf(search_file, "analyzed: %ld\n", pso->analyzed);
	fprintf(search_file, "number_parameters: %d\n", pso->number_parameters); 
	fwrite_bounds(search_file, pso->bounds);
	fprintf(search_file, "global_best_fitness: %.20lf\n", pso->global_best_fitness);
	print_double_array(search_file, "global_best", pso->number_parameters, pso->global_best);

	fprintf(search_file, "redundancies:\n");
	for (i = 0; i < pso->size; i++) {
		if (pso->redundancies[i] != NULL) {
			REDUNDANCY *r = pso->redundancies[i];
			while (r != NULL) {
				fwrite_redundancy(search_file, r, pso->number_parameters, i);
				r = r->next;
			}
		}
	}
	fclose(search_file);

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
	int particle, current, i;
	int count;
	REDUNDANCY *r, *current_r;

	(*search_data) = (PARTICLE_SWARM_OPTIMIZATION*)malloc(sizeof(PARTICLE_SWARM_OPTIMIZATION));
	pso = (PARTICLE_SWARM_OPTIMIZATION*)(*search_data);

	sprintf(search_filename, "%s/%s/search", get_working_directory(), search_name);
	search_file = fopen(search_filename, "r");
	if (search_file == NULL) return AS_READ_ERROR;

	fscanf(search_file, "current_particle: %d, size: %d\n", &(pso->current_particle), &(pso->size));
	fscanf(search_file, "w: %lf, c1: %lf, c2: %lf\n", &(pso->w), &(pso->c1), &(pso->c2));
	fscanf(search_file, "analyzed: %ld\n", &(pso->analyzed));
	fscanf(search_file, "number_parameters: %d\n", &(pso->number_parameters));

	fread_bounds(search_file, &(pso->bounds));

	fscanf(search_file, "global_best_fitness: %lf\n", &(pso->global_best_fitness));
	read_double_array(search_file, "global_best", &(pso->global_best));

	pso->redundancies = (REDUNDANCY**)malloc(sizeof(REDUNDANCY*) * pso->size);
	for (i = 0; i < pso->size; i++) pso->redundancies[i] = NULL;
	fscanf(search_file, "redundancies:\n");
	current = -1;
	current_r = NULL;
	count = 0;
	while (fread_redundancy(search_file, &r, pso->number_parameters, &particle)) {
		count++;
		if (particle != current) {
			current = particle;
			current_r = pso->redundancies[current];
		}

		if (current_r == NULL) {
			pso->redundancies[current] = r;
			current_r = r;
		} else {
			current_r->next = r;
			current_r = r;
		}
	}
	fclose(search_file);
	printf("read %d redundancies\n", count);

	sprintf(population_filename, "%s/%s/particles", get_working_directory(), search_name);
	if (0 > read_population(population_filename, &(pso->particles)) ) return AS_READ_ERROR;
	
	sprintf(population_filename, "%s/%s/velocities", get_working_directory(), search_name);
	if (0 > read_population(population_filename, &(pso->velocities)) ) return AS_READ_ERROR;

	sprintf(population_filename, "%s/%s/local_best", get_working_directory(), search_name);
	if (0 > read_population(population_filename, &(pso->local_best)) ) return AS_READ_ERROR;

	pso->remove_outliers = 1;
	pso->current_redundancy = NULL;
//	remove_outliers(pso->local_best, 0.15);

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

	if (!individual_exists(pso->local_best, pso->current_particle)) {
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
		if (pso->current_redundancy != NULL) {
			printf("current redundancy not null\n");
			/**
			 * Generate parameters from this redundancy and move the current redundancy
			 */
			particle = pso->current_redundancy->parameters;
			velocity = pso->current_redundancy->velocity;
			pso->current_redundancy = pso->current_redundancy->next;
		} else {
			printf("current redundancy is null\n");
			local_best = pso->local_best->individuals[pso->current_particle];
			velocity = pso->velocities->individuals[pso->current_particle];
			particle = pso->particles->individuals[pso->current_particle];

			for (i = 0; i < pso->number_parameters; i++) {
				velocity[i] = (pso->w * velocity[i]) + (pso->c1 * dsfmt_gv_genrand_close_open() * (local_best[i] - particle[i])) + (pso->c2 * dsfmt_gv_genrand_close_open() * (pso->global_best[i] - particle[i]));
			}
			bound_velocity(particle, velocity, pso->bounds);

			for (i = 0; i < pso->number_parameters; i++) {
				particle[i] = particle[i] + velocity[i];
				if (isnan(particle[i])) return AS_GEN_FAIL;
			}
			bound_parameters(particle, pso->bounds);
			
			pso->current_particle++;
			if (pso->current_particle >= pso->size) pso->current_particle = 0;
			pso->current_redundancy = pso->redundancies[pso->current_particle];
		}

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

	if (isnan(sp->fitness)) return AS_INSERT_FITNESS_NAN;

	for (i = 0; i < pso->number_parameters; i++) {
		if (isnan(sp->parameters[i])) return AS_INSERT_PARAMETERS_NAN;
		if (sp->parameters[i] < pso->bounds->min_bound[i] || sp->parameters[i] > pso->bounds->max_bound[i]) return AS_INSERT_OUT_OF_BOUNDS;
	}

	if (1 != sscanf(sp->metadata, "p: %d, v:", particle)) return AS_INSERT_INVALID_METADATA;
	current_token = strtok(strchr(sp->metadata, 'v'), ", :");
	for (i = 0; i < pso->number_parameters; i++) {
		if (current_token == NULL) return AS_INSERT_INVALID_METADATA;
		velocity[i] = atof(current_token);
		current_token = strtok(NULL, " ");
	}

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

int check_particle(PARTICLE_SWARM_OPTIMIZATION *pso, int particle, double *velocity, SEARCH_PARAMETERS *sp) {
	/**
	 * See if this particle has any saved local best values
	 * If this matches a local best value update the local best
	 * 	remove the match from the queue
	 * 	remove all queued matches with lower fitness
	 * If this matches parameters but not fitness, remove the match from the queue
	 */
	REDUNDANCY *r, *r_prev, *n, *n_prev;
	int match;
	r = pso->redundancies[particle];
	r_prev = NULL;
	while (r != NULL) {
		printf("checking parameters\n");
		if (parameters_match(pso->number_parameters, r->parameters, sp->parameters)) {
			match = fitness_match(r->fitness, sp->fitness);
			printf ("parameter match for fitness: %.20lf, fitness match %.20lf? %d\n", r->fitness, sp->fitness, match);

			if (match) {
				n_prev = NULL;
				n = pso->redundancies[particle];
				while (n != NULL) {
					if (n->fitness <= sp->fitness) {
						REDUNDANCY *temp;
						printf("removing redundancy with fitness: %.20lf\n", n->fitness);
						temp = n->next;
						if (n_prev == NULL) {
							printf("n_prev is null!\n");
							pso->redundancies[particle] = n->next;
						} else {
							printf("n_prev is not null!\n");
							n_prev->next = n->next;
						}
						printf("freeing!\n");
						free_redundancy(&n);
						printf("updating n\n");
						n = temp;
					} else {
						printf("not removing fitness: %.20lf\n", n->fitness);
						n_prev = n;
						n = n->next;
					}
				}
			} else {
				printf("removing redundancy with fitness: %.20lf\n", r->fitness);
				if (r_prev == NULL) {
					pso->redundancies[particle] = r->next;
					free_redundancy(&r);
				} else {
					r_prev->next = r->next;
					free_redundancy(&r);
				}
			}
			pso->current_redundancy = NULL;

			return match;
		}
		r_prev = r;
		r = r->next;
	}
	return 0;
}

int verify_particle(PARTICLE_SWARM_OPTIMIZATION *pso, int particle, double *velocity, SEARCH_PARAMETERS *sp) {
	REDUNDANCY *r;
	r = pso->redundancies[particle];
	if (r == NULL) {
		new_redundancy(&(pso->redundancies[particle]), sp->fitness, pso->number_parameters, sp->parameters, velocity);
	} else {
		while (r->next != NULL) r = r->next;
		new_redundancy(&(r->next), sp->fitness, pso->number_parameters, sp->parameters, velocity);
	}
	return 0;
}

int pso_insert_parameters(char* search_name, void* search_data, SEARCH_PARAMETERS* sp) {
	PARTICLE_SWARM_OPTIMIZATION *pso = (PARTICLE_SWARM_OPTIMIZATION*)search_data;
	int result, particle;
	double *velocity = (double*)malloc(sizeof(double) * pso->number_parameters);

	result = parse(pso, sp, &particle, velocity);
	if (result != 0) {
		free(velocity);
		sprintf(AS_MSG, "parse error");
		return result;
	}

	pso->analyzed++;

	if (!population_contains(pso->local_best, sp->fitness, sp->parameters)) {
		/**
		 * Verify the particle before comparing fitness for insert, so we can remove bad redundancies.
		 */

		if (pso->local_best->size >= 40) {
			double min_error, max_error, median_error, average_error;
			double particle_error;
			population_error_stats(pso->local_best, &min_error, &max_error, &median_error, &average_error);
			particle_error = distance_error_from(pso->local_best, sp->fitness, sp->parameters);
			if (particle_error > average_error * 30.0) {
				sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, outlier", particle, sp->fitness, pso->local_best->fitness[particle], pso->global_best_fitness);
				return AS_INSERT_OUTLIER;
			}
		}

		if (!individual_exists(pso->local_best, particle)) {
			insert_particle(search_name, pso, particle, velocity, sp);
		} else {
			if (check_particle(pso, particle, velocity, sp)) {
				insert_particle(search_name, pso, particle, velocity, sp);
			} else {
				if (sp->fitness > pso->local_best->fitness[particle]) {
					verify_particle(pso, particle, velocity, sp);
					sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, verifying", particle, sp->fitness, pso->local_best->fitness[particle], pso->global_best_fitness);
				} else {
					sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, low fitness", particle, sp->fitness, pso->local_best->fitness[particle], pso->global_best_fitness);
				}
			}
		}	
	} else {
		sprintf(AS_MSG, "p[%d]: %.15lf, l: %.15lf, g: %.15lf, duplicate", particle, sp->fitness, pso->local_best->fitness[particle], pso->global_best_fitness);
	}
	return AS_INSERT_SUCCESS;
}
