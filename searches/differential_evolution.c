#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "differential_evolution.h"
#include "recombination.h"
#include "population.h"

#include "../evaluation/evaluator.h"

#define DE_RAND 0
#define DE_BEST 1
#define DE_CURRENT_TO_RAND 2
#define DE_CURRENT_TO_BEST 3

#define DE_PAIR 0
#define DE_DIR 1

#define DE_BINOMIAL 0
#define DE_EXPONENTIAL 1
#define DE_DIRECTIONAL 2
#define DE_NONE 3

#define metadata_size 2048

void invalid_differential_evolution() {
	fprintf(stderr, "Differential Evolution was illdefined:\n");
	fprintf(stderr, "Search naming conventions are:\n");
	fprintf(stderr, "  de/<parent_selection>/<number_pairs>/<pair_recombination>/<population_size>/<max_evaluations>\n");
	fprintf(stderr, "    parent_selection:\n");
	fprintf(stderr, "      best\n");
	fprintf(stderr, "      rand\n");
	fprintf(stderr, "      current_to_best(parent_scale=<parent_scale>)\n");
	fprintf(stderr, "      current_to_rand(parent_scale=<parent_scale>)\n");
	fprintf(stderr, "    pair_recombination:\n");
	fprintf(stderr, "      binomial(crossover_rate=<crossover_rate>,crossover_scale=<crossover_scale|random>)\n");
	fprintf(stderr, "      exponential(crossover_rate=<crossover_rate>,crossover_scale=<crossover_scale|random>)\n");
	fprintf(stderr, "      directional(crossover_rate=<crossover_rate>,crossover_scale=<crossover_scale|random>)\n");
	fprintf(stderr, "      none\n");
}

void parse_differential_evolution(char search_parameters[512], DIFFERENTIAL_EVOLUTION *de, int *population_size, int *max_evaluations) {
	char rest[512];
	char rest2[512];

	printf("search_parameters: %s\n", search_parameters);

	if (sscanf(search_parameters, "de/best/%d/%s", &de->number_pairs, rest) == 2) de->parent_type = DE_BEST;
	else if (sscanf(search_parameters, "de/rand/%d/%s", &de->number_pairs, rest) == 2) de->parent_type = DE_RAND;
	else if (sscanf(search_parameters, "de/current_to_best(parent_scale=%lf)/%d/%s", &de->parent_scale, &de->number_pairs, rest) == 3) de->parent_type = DE_CURRENT_TO_BEST;
	else if (sscanf(search_parameters, "de/current_to_rand(parent_scale=%lf)/%d/%s", &de->parent_scale, &de->number_pairs, rest) == 3) de->parent_type = DE_CURRENT_TO_RAND;
	else {
		invalid_differential_evolution();
		return;
	}
	printf("rest: %s, pairs: %d\n", rest, de->number_pairs);

	de->pair_scale = -1;
	if (sscanf(rest, "none/%s", rest2) == 1) de->recombination_type = DE_NONE;
	else if (sscanf(rest, "binomial(crossover_rate=%lf,crossover_scale=random)/%s", &de->crossover_rate, rest2) == 2) de->recombination_type = DE_BINOMIAL;
	else if (sscanf(rest, "binomial(crossover_rate=%lf,crossover_scale=%lf)/%s", &de->crossover_rate, &de->pair_scale, rest2) == 3) de->recombination_type = DE_BINOMIAL;
	else if (sscanf(rest, "exponential(crossover_rate=%lf,crossover_scale=random)/%s", &de->crossover_rate, rest2) == 2) de->recombination_type = DE_EXPONENTIAL;
	else if (sscanf(rest, "exponential(crossover_rate=%lf,crossover_scale=%lf)/%s", &de->crossover_rate, &de->pair_scale, rest2) == 3) de->recombination_type = DE_EXPONENTIAL;
	else if (sscanf(rest, "directional(crossover_rate=%lf,crossover_scale=random)/%s", &de->crossover_rate, rest2) == 2) de->recombination_type = DE_DIRECTIONAL;
	else if (sscanf(rest, "directional(crossover_rate=%lf,crossover_scale=%lf)/%s", &de->crossover_rate, &de->pair_scale, rest2) == 3) de->recombination_type = DE_DIRECTIONAL;
	else {
		invalid_differential_evolution();
		return;
	}
	printf("rest2: %s\n", rest2);

	if (sscanf(rest2, "%d/%d", population_size, max_evaluations) != 2) {
		invalid_differential_evolution();
		return;
	}

	de->current_best = -1;
}

void print_differential_evolution(DIFFERENTIAL_EVOLUTION *de) {
	printf("current best: %d\n", de->current_best);
	printf("parent_type: %d, parent_scale: %lf\n", de->parent_type, de->parent_scale);
	printf("number_pairs: %d, pair_type: %d, pair_scale: %lf\n", de->number_pairs, de->pair_type, de->pair_scale);
	printf("recombination_type: %d, crossover_rate: %lf\n", de->recombination_type, de->crossover_rate);
	printf("next_generated: %d\n", de->next_generated);
}


void start_differential_evolution(char search_path[512], char search_parameters[512], double* min_parameters, double* max_parameters, int number_parameters, POPULATION **population) {
	DIFFERENTIAL_EVOLUTION* de;
	int population_size;
	int max_evaluations;

	de = (DIFFERENTIAL_EVOLUTION*)malloc(sizeof(DIFFERENTIAL_EVOLUTION));
	parse_differential_evolution(search_parameters, de, &population_size, &max_evaluations);
	(*population) = new_population(search_path, search_parameters, min_parameters, max_parameters, number_parameters, population_size, max_evaluations);
	(*population)->parameters = de;
	(*population)->get_individual = differential_evolution__get_individual;
	(*population)->insert_individual = differential_evolution__insert_individual;
}


void differential_evolution__process_metadata(char *metadata, int *position) {
	sscanf(metadata, "position:%d", position);
}

void differential_evolution__insert_individual(POPULATION* population, double* parameters, double fitness, char *metadata) {
	DIFFERENTIAL_EVOLUTION *de;
	int position, i;

	de = (DIFFERENTIAL_EVOLUTION*)population->parameters;
	if (de->current_best < 0) {
		de->current_best = 0;
		for (i = 0; i < population->current_size; i++) {
			if (population->fitness[de->current_best] > population->fitness[i]) de->current_best = i;
		}
	}

	population->current_evaluation++;
	differential_evolution__process_metadata(metadata, &position);

	replace_if_better(population, position, parameters, fitness);
	if (population->fitness[de->current_best] > fitness) de->current_best = position;
}

void differential_evolution__get_individual(POPULATION *population, double **parameters, char **metadata) {
	DIFFERENTIAL_EVOLUTION *de;
	double *parent;
	double *pair_sum;
	int i, target;

	de = (DIFFERENTIAL_EVOLUTION*)population->parameters;
	(*metadata) = (char*)malloc(sizeof(char) * metadata_size);
	if (population->current_size < population->max_size) {
		(*parameters) = random_recombination(population->min_parameters, population->max_parameters, population->number_parameters);
		sprintf((*metadata), "position:%d,random", de->next_generated);
		de->next_generated++;
		if (de->next_generated >= population->max_size) de->next_generated = 0;
		return;
	}

	parent = (double*)malloc(sizeof(double) * population->number_parameters);
	if (de->parent_type == DE_RAND) {
		target = drand48() * population->current_size;
		for (i = 0; i < population->number_parameters; i++) parent[i] = population->individuals[target][i];
	} else if (de->parent_type == DE_BEST) {
		for (i = 0; i < population->number_parameters; i++) parent[i] = population->individuals[de->current_best][i];
	} else if (de->parent_type == DE_CURRENT_TO_RAND) {
		target = drand48() * population->current_size;
		for (i = 0; i < population->number_parameters; i++) {
			parent[i] = population->individuals[de->next_generated][i] + de->parent_scale * (population->individuals[target][i] - population->individuals[de->next_generated][i]);
		}
	} else if (de->parent_type == DE_CURRENT_TO_BEST) {
		for (i = 0; i < population->number_parameters; i++) {
			parent[i] = population->individuals[de->next_generated][i] + de->parent_scale * (population->individuals[de->current_best][i] - population->individuals[de->next_generated][i]);
		}
	} else {
		fprintf(stderr, "ERROR: unknown parent type for differential evolution.\n");
		return;
	}

	if (de->pair_type == DE_PAIR) {
		pair_sum = get_pair_sum(population->individuals, population->current_size, population->number_parameters, de->number_pairs, de->pair_scale);
	} else if (de->pair_type == DE_DIR) {
		pair_sum = get_dir_sum(population->individuals, population->fitness, population->current_size, population->number_parameters, de->number_pairs, de->pair_scale);
	} else {
		fprintf(stderr, "ERROR: unknown pair type for differential evolution.\n");
		return;
	}

	(*parameters) = (double*)malloc(sizeof(double) * population->number_parameters);
	if (de->recombination_type == DE_BINOMIAL) {
		target = drand48() * population->number_parameters;
		for (i = 0; i < population->number_parameters; i++) {
			if (i == target || drand48() < de->crossover_rate) (*parameters)[i] = parent[i] + pair_sum[i];
			else (*parameters)[i] = population->individuals[de->next_generated][i];
		}
	} else if (de->recombination_type == DE_EXPONENTIAL) {
		target = drand48() * population->number_parameters;
		for (i = 0; i < population->number_parameters; i++) {
			if (i == target || drand48() < de->crossover_rate) break;
			else (*parameters)[i] = population->individuals[de->next_generated][i];
		}
		for (; i < population->number_parameters; i++) (*parameters)[i] = parent[i] + pair_sum[i];
	} else if (de->recombination_type == DE_NONE || de->recombination_type == DE_DIRECTIONAL) {
		for (i = 0; i < population->number_parameters; i++) (*parameters)[i] = parent[i] + pair_sum[i];
	} else {
		fprintf(stderr, "ERRROR: unknown recombination type for differential evolution.\n");
		return;
	}
	sprintf((*metadata), "position:%d,parent_type:%d,pair_type:%d,recombination_type:%d", de->next_generated, de->parent_type, de->pair_type, de->recombination_type);
	de->next_generated++;

	if (de->next_generated >= population->max_size) de->next_generated = 0;

	free(parent);
	free(pair_sum);
}
