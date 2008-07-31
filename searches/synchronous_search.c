#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "population.h"
#include "../evaluation/evaluator.h"

void synchronous_search(char *search_path, char *search_parameters, double *min_parameters, double *max_parameters, int number_parameters, void (*init_search)(char*, char*, double*, double*, int, POPULATION**)) {
	POPULATION *population;
	double fitness;
	double *individual;
	char *metadata;

	init_search(search_path, search_parameters, min_parameters, max_parameters, number_parameters, &population);
	while (population->current_evaluation < population->max_evaluations) {
		(population->get_individual)(population, &individual, &metadata);
		fitness = evaluate(individual);
		(population->insert_individual)(population, individual, fitness, metadata);

		printf("iteration %d: current: %lf, ", population->current_evaluation, fitness);
		fwrite_population_statistics(stdout, population);
		printf("\n");

		free(individual);
		free(metadata);
	}
}

/*
void synchronous_parallel_search(POPULATION* population, double (*likelihood_function)(double*), double (*likelihood_compose)(double*, int), int number_workers) {
	double fitness;
	double **individual;
	char **metadata;

	evaluator_multiple(likelihood_function, likelihood_compose, number_workers);
	while (population->current_evaluation < population->max_evaluations) {
		for (i = 0; i < number_workers; i++) {
			(population->get_individual)(population, &individual[i], &metadata[i]);
		}

		fitness = evaluate(individual);

		for (i = 0; i < number_workers; i++) {
			(population->insert_individual)(population, individuals[i], fitness[i], metadata[i]);
			printf("iteration %d: current: %lf, ", population->current_evaluation, fitness[i]);
			fwrite_population_statistics(stdout, population);
			printf("\n");
		}
		for (i = 0; i < number_workers; i++) {
			free(individuals[i]);
			free(metadata[i]);
		}
		free(individual);
		free(metadata);
		free(fitness);
	}
}
*/
