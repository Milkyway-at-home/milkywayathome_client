#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/********
	*	FGDO includes
 ********/
#include "../evaluation/boinc_search_manager.h"
#include "../evaluation/search_manager.h"
#include "../searches/asynchronous_newton_method.h"
#include "../searches/asynchronous_differential_evolution.h"
#include "../searches/asynchronous_particle_swarm.h"
#include "../searches/asynchronous_genetic_search.h"

void print_arguments() {
	printf("Usage:\n");
	printf("\t-cwd <working_directory>, default: ./\n");
	printf("\t-wus <wus_to_generate>, required.\n");
	exit(1);
}

int main(int argc, char** argv) {
	register_search(get_asynchronous_newton_method());
	register_search(get_asynchronous_differential_evolution());
	register_search(get_asynchronous_particle_swarm());
	register_search(get_asynchronous_genetic_search());

	init_boinc_search_manager(argc, argv);
	start_search_manager();
}
