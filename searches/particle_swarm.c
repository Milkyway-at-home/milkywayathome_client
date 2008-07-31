#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "particle_swarm.h"
#include "recombination.h"
#include "population.h"

#include "../evaluation/evaluator.h"

void invalid_particle_swarm() {
	fprintf(stderr, "Particle swarm was illdefined:\n");
	fprintf(stderr, "Search naming conventions are:\n");
	fprintf(stderr, "  ps/<velocity_weight>/<local_weight>/<global_weight>/<population_size>/<max_evaluations>\n");
}

void parse_particle_swarm(char search_parameters[512], PARTICLE_SWARM *ps, int *population_size, int *max_evaluations) {
	printf("search_parameters: %s\n", search_parameters);

	if (sscanf(search_parameters, "ps/%lf/%lf/%lf/%d/%d", &ps->velocity_weight, &ps->local_weight, &ps->global_weight, population_size, max_evaluations) != 5) {
		(*population_size) = 0;
		(*max_evaluations) = 0;
		invalid_particle_swarm();
		return;
	}
	ps->next_generated = 0;
}

void start_particle_swarm(char search_path[512], char search_parameters[512], double* min_parameters, double* max_parameters, int number_parameters, POPULATION **population) {
	PARTICLE_SWARM *ps;
	FILE* ps_file;
	char ps_path[1024];
	int population_size;
	int max_evaluations;
	int i, j;

	ps = (PARTICLE_SWARM*)malloc(sizeof(PARTICLE_SWARM));
	parse_particle_swarm(search_parameters, ps, &population_size, &max_evaluations);
	(*population) = new_population(search_path, search_parameters, min_parameters, max_parameters, number_parameters, population_size, max_evaluations);
	(*population)->parameters = ps;
	(*population)->get_individual = particle_swarm__get_individual;
	(*population)->insert_individual = particle_swarm__insert_individual;

	ps->velocity = (double**)malloc(sizeof(double*) * population_size);
	ps->local_best_parameters = (double**)malloc(sizeof(double*) * population_size);
	ps->local_best_fitness = (double*)malloc(sizeof(double) * population_size);

	sprintf(ps_path, "%s/particle_swarm", search_path);
	ps_file = fopen(ps_path, "r");
	if (ps_file) {
		//Search exists.
		fscanf(ps_file, "global best fitnss: %lf\n", &ps->global_best_fitness);
		fscanf(ps_file, "global best:");
		for (i = 0; i < (*population)->number_parameters; i++) fscanf(ps_file, " %lf", &ps->global_best_parameters[i]);
		fscanf(ps_file, "\n");

		fscanf(ps_file, "local best fitness : local best parameters : velocity\n");
		for (i = 0; i < (*population)->current_size; i++) {
			fscanf(ps_file, "%lf :", &ps->local_best_fitness[i]);
			ps->local_best_parameters[i] = (double*)malloc(sizeof(double) * (*population)->number_parameters);
			for (j = 0; j < (*population)->number_parameters; j++) fscanf(ps_file, " %lf", &ps->local_best_parameters[i][j]);
			fscanf(ps_file, " :");
			ps->velocity[i] = (double*)malloc(sizeof(double) * (*population)->number_parameters);
			for (j = 0; j < (*population)->number_parameters; j++) fscanf(ps_file, " %lf", &ps->velocity[i][j]);
			fscanf(ps_file, "\n");
		}
	}
}


void particle_swarm__process_metadata(char *metadata, int *position) {
	sscanf(metadata, "position:%d", position);
}

void particle_swarm__insert_individual(POPULATION* population, double* parameters, double fitness, char *metadata) {
	PARTICLE_SWARM *ps;
	int position, i;

	ps = (PARTICLE_SWARM*)population->parameters;
	if (ps->global_best_parameters == NULL) {
		ps->global_best_parameters = (double*)malloc(sizeof(double) * population->number_parameters);
		ps->global_best_fitness = fitness;
		for (i = 0; i < population->number_parameters; i++) ps->global_best_parameters[i] = parameters[i];
		printf("updated global best[%d] to: %lf\n", position, fitness);
	}

	population->current_evaluation++;
	particle_swarm__process_metadata(metadata, &position);

	if (ps->local_best_parameters[position] == NULL) {
		ps->local_best_parameters[position] = (double*)malloc(sizeof(double) * population->number_parameters);
		ps->local_best_fitness[position] = fitness;
		for (i = 0; i < population->number_parameters; i++) ps->local_best_parameters[position][i] = parameters[i];
		printf("updated local best[%d] to: %lf\n", position, fitness);
	} else if (ps->local_best_fitness[position] > fitness) {
		ps->local_best_fitness[position] = fitness;
		for (i = 0; i < population->number_parameters; i++) ps->local_best_parameters[position][i] = parameters[i];
		printf("updated local best[%d] to: %lf\n", position, fitness);
	}

	replace(population, position, parameters, fitness);
	if (ps->global_best_fitness > fitness) {
		ps->global_best_fitness = fitness;
		for (i = 0; i < population->number_parameters; i++) ps->global_best_parameters[i] = parameters[i];
		printf("updated global best[%d] to: %lf\n", position, fitness);
	}
}

void particle_swarm__get_individual(POPULATION *population, double **parameters, char **metadata) {
	PARTICLE_SWARM *ps;
	double* velocity;
	double* parent;
	double* global_best;
	double* local_best;
	double rand1, rand2;
	int i;

	ps = (PARTICLE_SWARM*)population->parameters;
	(*metadata) = (char*)malloc(sizeof(char) * metadata_size);
	if (population->current_size < population->max_size) {
		(*parameters) = random_recombination(population->min_parameters, population->max_parameters, population->number_parameters);
		sprintf((*metadata), "position:%d,random", ps->next_generated);
		ps->next_generated++;
		if (ps->next_generated >= population->max_size) ps->next_generated = 0;
		return;
	}

	if (ps->velocity[ps->next_generated] == NULL) {
		ps->velocity[ps->next_generated] = random_recombination(population->min_parameters, population->max_parameters, population->number_parameters);
		for (i = 0; i < population->number_parameters; i++) {
			ps->velocity[ps->next_generated][i] = population->individuals[ps->next_generated][i] - ps->velocity[ps->next_generated][i];
		}
	}

	parent = population->individuals[ps->next_generated];
	velocity = ps->velocity[ps->next_generated];
	global_best = ps->global_best_parameters;
	local_best = ps->local_best_parameters[ps->next_generated];

	rand1 = drand48();
	rand2 = drand48();
	for (i = 0; i < population->number_parameters; i++) {
		velocity[i] = (ps->velocity_weight * velocity[i]) + (ps->local_weight * rand1 * (local_best[i] - parent[i])) + (ps->global_weight * rand2 * (global_best[i] - parent[i]));
	}

	(*parameters) = (double*)malloc(sizeof(double) * population->number_parameters);
	for (i = 0; i < population->number_parameters; i++) {
		(*parameters)[i] = velocity[i] + parent[i];
	}
	bound_parameters(population, (*parameters));
	sprintf((*metadata), "position:%d,pso", ps->next_generated);

	ps->next_generated++;
	if (ps->next_generated >= population->max_size) ps->next_generated = 0;
}
