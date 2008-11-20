#include <stdlib.h>

#include "evaluator.h"
#include "mpi_search_manager.h"
#include "search_manager.h"
#include "../searches/search_parameters.h"


void start_mpi_search_manager(int argc, char** argv) {
	SEARCH_PARAMETERS **sp;
	int generated, generation_rate, i;

	printf("[mpi_search_manager]: init.\n");
	init_search_manager(argc, argv);
	printf("[mpi_search_manager]: successful.\n");

	generated = 1;
	generation_rate = get_generation_rate();
	sp = (SEARCH_PARAMETERS**)malloc(sizeof(SEARCH_PARAMETERS*) * generation_rate);

	while (generated > 0) {
		generated = generate_search_parameters(sp);
		for (i = 0; i < generated; i++) {
			sp[i]->fitness = evaluate(sp[i]->parameters);
			insert_search_parameters(sp[i]);
			free_search_parameters(sp[i]);
		}
	}
}
