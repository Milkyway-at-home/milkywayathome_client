#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/********
	*	FGDO includes
 ********/
#include "../evaluation/search_manager.h"
#include "../searches/genetic_search.h"
#include "../searches/newton_method.h"
#include "../searches/differential_evolution.h"
#include "../searches/particle_swarm.h"
#include "../searches/search_parameters.h"

/********
	*	Astronomy includes
 ********/
#include "add_workunit.h"


void print_arguments() {
	printf("Usage:\n");
	printf("\t-d <working_directory>, default: ./\n");
	printf("\t-s <search_name>, required.\n");
	printf("\t-wus <wus_to_generate>, required.\n");
	exit(1);
}

int main(int argc, char** argv) {
	int wus_to_generate, i;

	/********
		*	Arguments
		*		-d <working_directory>
		*		-s <search_name>
		*		-wus <wus_to_generate>
	 ********/
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-wus") == 0) wus_to_generate = atoi(argv[++i]);
		if (strcmp(argv[i], "-s") == 0) strcpy(search_name, argv[++i]);
	}

	register_search("nm", start_newton_method);
	register_search("gs", start_genetic_search);
	register_search("de", start_differential_evolution);
	register_search("pso", start_particle_swarm);

	init_add_workunit();

	init_search_manager(argc, argv, add_workunit);
	manage_search(search_name);
	generate_workunits(search_name, wus_to_generate);
}
