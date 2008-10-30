#include "boinc_db.h"
#include "error_numbers.h"
#include "backend_lib.h"
#include "parse.h"
#include "util.h"
          
#include "sched_config.h"
#include "sched_util.h"
#include "sched_msgs.h"

#include "../settings.h"
#include "../searches/search_parameters.h"
#include "../searches/genetic_search.h"
#include "../searches/differential_evolution.h"
#include "../searches/particle_swarm.h"

#include <string.h>
#include <iostream>
#include <sstream>

int number_searches = 0;
POPULATION** searches;


void init_search_manager(int argc, char** argv, 


POPULATION* get_search(char* search_name) {
	char search_type[512];
	char population_path[2048], population_file_name[2048];
	FILE* population_file;
	int i, population_size, max_evaluations;

	sscanf(search_name, "%s_", search_type);
	sprintf(population_file_name, "%s%s/population", BOINC_SEARCH_PATH, search_name);
	sprintf(population_path, "%s%s", BOINC_SEARCH_PATH, search_name);
	for (i = 0; i < number_searches; i++) {
		if (!strcmp(searches[i]->search_path, population_path)) {
			return searches[i];
		}
	}

	population_file = fopen(population_file_name, "r");
	if (population_file != NULL) {
		searches = (POPULATION**)realloc(searches, sizeof(POPULATION*) * (number_searches + 1));
		searches[number_searches] = fread_population(population_file);

		if (!strcmp(search_type, "gs")) {
		        searches[number_searches]->parameters = (GENETIC_SEARCH*)malloc(sizeof(GENETIC_SEARCH));
			parse_genetic_search(searches[number_searches]->search_parameters, (GENETIC_SEARCH*)searches[number_searches]->parameters, &population_size, &max_evaluations);
			searches[number_searches]->get_individual = genetic_search__get_individual;
			searches[number_searches]->insert_individual = genetic_search__insert_individual;
		} else if (!strcmp(search_type, "de")) {
		        searches[number_searches]->parameters = (DIFFERENTIAL_EVOLUTION*)malloc(sizeof(DIFFERENTIAL_EVOLUTION));
			parse_differential_evolution(searches[number_searches]->search_parameters, (DIFFERENTIAL_EVOLUTION*)searches[number_searches]->parameters, &population_size, &max_evaluations);
			searches[number_searches]->get_individual = differential_evolution__get_individual;
			searches[number_searches]->insert_individual = differential_evolution__insert_individual;
		} else if (!strcmp(search_type, "pso")) {
		        searches[number_searches]->parameters = (PARTICLE_SWARM*)malloc(sizeof(PARTICLE_SWARM));
			parse_particle_swarm(searches[number_searches]->search_parameters, (PARTICLE_SWARM*)searches[number_searches]->parameters, &population_size, &max_evaluations);
			searches[number_searches]->get_individual = particle_swarm__get_individual;
			searches[number_searches]->insert_individual = particle_swarm__insert_individual;
		} else if (!strcmp(search_type, "newton")) {
			return NULL;
		}

		number_searches++;
		return searches[number_searches-1];
	} else {
		fprintf(stderr, "ERROR: unknown search: %s\n", population_file_name);
		return NULL;
	}
}


int insert_workunit(char* search_name, double fitness, double* parameters, char* metadata) {
	POPULATION* target;

	target = get_search(search_name);
	if (target == NULL) return 0;
	target->insert_individual(target, parameters, fitness, metadata);
	return 1;
}

int generate_workunits(int number_workunits, SEARCH_PARAMETERS*** parameters) {
	int generated, current, i, j;

	(*parameters) = (SEARCH_PARAMETERS**)malloc(sizeof(SEARCH_PARAMETERS*) * number_workunits);

	generated = 0;
	current = 0;
	for (i = 0; i < number_searches; i++) {
		generated = (number_workunits - current) / (number_searches - i);
		for (j = 0; j < generated; j++) {
			strcpy((*parameters)[current]->search_path, searches[i]->search_path);

			(*parameters)[current]->number_parameters = searches[i]->number_parameters;
			searches[i]->get_individual(searches[i], &((*parameters)[current]->parameters), &((*parameters)[current]->metadata));
			current++;
		}
	}
	return current;
}
