#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/********
	*	FGDO includes
 ********/
#include "../evaluation/search_manager.h"
#include "../evaluation/boinc_search_manager.h"
#include "../searches/asynchronous_newton_method.h"
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
	int i;
	char search_name[512];

	/********
		*	Arguments
		*		-s <search_name>
		*		-wus <wus_to_generate>
	 ********/
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-s") == 0) strcpy(search_name, argv[++i]);
	}

	register_search("nm", init_newton_method);
//	register_search("gs", init_genetic_search);
//	register_search("de", init_differential_evolution);
//	register_search("pso", init_particle_swarm);

	init_add_workunit();
	init_boinc_search_manager(argc, argv, add_workunit);
	generate_workunits();
}
