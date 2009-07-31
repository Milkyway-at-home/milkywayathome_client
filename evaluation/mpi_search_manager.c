/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "evaluator.h"
#include "mpi_search_manager.h"
#include "search_manager.h"
#include "../searches/asynchronous_search.h"
#include "../searches/bounds.h"
#include "../searches/search_parameters.h"


void start_mpi_search_manager(int argc, char** argv, int number_parameters) {
	SEARCH_PARAMETERS **sp;
	int generated, generation_rate, i, total;

	printf("[mpi_search_manager]: init.\n");
	init_search_manager(argc, argv);
	printf("[mpi_search_manager]: successful.\n");

	generated = 1;
	generation_rate = get_generation_rate();
	printf("[mpi_search_manager]: generation_rate: %d\n", generation_rate);

	sp = (SEARCH_PARAMETERS**)malloc(sizeof(SEARCH_PARAMETERS*) * generation_rate);
	for (i = 0; i < generation_rate; i++) {
		init_search_parameters(&(sp[i]), number_parameters);
	}

	total = 0;
	while (generated > 0) {
//		printf("[mpi_search_manager]: generating %d search parameters.\n", generation_rate);
		generated = generate_search_parameters(sp);
//		printf("[mpi_search_manager]: %d generated.\n", generated);
//		printf("[mpi_search_manager]: evaluating %d search parameters.\n", generated);
		for (i = 0; i < generated; i++) {
			sp[i]->fitness = evaluate(sp[i]->parameters);
			insert_search_parameters(sp[i]);
			total++;
			printf("[%d]: %.20lf -- %s\n", total, sp[i]->fitness, AS_MSG);
		}
//		printf("[mpi_search_manager]: evaluated %d search parameters.\n", total);
	}
}


void mpi_asynchronous_search(int number_arguments, char** arguments, int number_parameters, double *point, double *range, BOUNDS *bounds) {
	int i, generate_result, insert_result, retval;
	MANAGED_SEARCH *ms = NULL;
	SEARCH_PARAMETERS *sp;

	for (i = 0; i < number_arguments; i++) {
		printf("checking argument: %s\n", arguments[i]);

		if (!strcmp(arguments[i], "-s")) {
			char *qualifier;
			i++;
			get_qualifier_from_name(arguments[i], &qualifier);
			if (!search_exists(arguments[i])) {
				ASYNCHRONOUS_SEARCH *as;
				printf("getting search from qualifier: %s\n", qualifier);
				as = get_registered_search(qualifier);
				printf("creating search\n");
				retval = as->create_search(arguments[i], number_arguments, arguments, number_parameters, point, range, bounds);
				if (retval) {
					printf("ERROR creating search: %s -- [%s]\n", arguments[i], AS_MSG);
                                        return;
				}
			}
                        manage_search(arguments[i]);
                        free(qualifier);
			ms = get_search(arguments[i]);
			break;
		}
	}
	if (ms == NULL) return;

	printf("Searching...\n");

	init_search_parameters(&sp, number_parameters);
	sprintf(sp->search_name, "%s", ms->search_name);

	i = 0;
	generate_result = AS_GEN_SUCCESS;
	while (generate_result == AS_GEN_SUCCESS) {
		generate_result = ms->search->generate_parameters(ms->search_name, ms->search_data, sp);
		sp->fitness = evaluate(sp->parameters);
		insert_result = ms->search->insert_parameters(ms->search_name, ms->search_data, sp);
		printf("[%d] %s, %s, %s\n", i, AS_MSG, AS_INSERT_STR[insert_result], sp->metadata);
		i++;
		if (i % 100 == 0) ms->search->checkpoint_search(ms->search_name, ms->search_data);
		free(sp->parameters);
		sp->number_parameters = number_parameters;
		sp->parameters = malloc(sizeof(double) * number_parameters);
	}
}
