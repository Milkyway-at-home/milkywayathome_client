#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/********
	*	FGDO includes
 ********/
#include "../evaluation/boinc_search_manager.h"
#include "../evaluation/search_manager.h"
#include "../searches/asynchronous_newton_method.h"

void print_arguments() {
	printf("Usage:\n");
	printf("\t-cwd <working_directory>, default: ./\n");
	printf("\t-wus <wus_to_generate>, required.\n");
	exit(1);
}

int main(int argc, char** argv) {
	register_search(asynchronous_newton_method);
//	register_search("gs", start_genetic_search);
//	register_search("de", start_differential_evolution);
//	register_search("pso", start_particle_swarm);

	init_boinc_search_manager(argc, argv);
	start_search_manager();
}
