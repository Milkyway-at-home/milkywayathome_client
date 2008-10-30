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
	register_search("nm", start_newton_method);
	register_search("gs", start_genetic_search);
	register_search("de", start_differential_evolution);
	register_search("pso", start_particle_swarm);

	init_add_workunit();

	init_search_manager(argc, argv, add_workunit);
	start_search_manager();
}
