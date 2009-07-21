/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

/********
	*	FGDO includes
 ********/
#include "asynchronous_search.h"
#include "../evaluation/search_manager.h"
#include "../evaluation/evaluator.h"

const char *AS_GEN_STR[] = { "success", "search completed", "ERROR", "failure" };
const char *AS_INSERT_STR[] = { "success", "search completed", "fitness is NAN", "fitness is invalid", "parameters contain NAN", "parameters out of bounds", "ERROR", "out of range", "out of iteration", "bad metadata", "not unique", "outlier", "invalid metadata"};
const char *AS_CP_STR[] = { "success", "search completed", "ERROR" };

char AS_MSG[1024] = "";
char AS_VERIFY_MSG[1024] = "";

void asynchronous_search__init(int number_arguments, char **arguments, int number_parameters, double *point, double *range, BOUNDS* bounds) {
	int i, retval;

	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-s")) {
			char *qualifier;
			i++;
			if (!search_exists(arguments[i])) {
				get_qualifier_from_name(arguments[i], &qualifier);
				retval = (get_registered_search(qualifier))->create_search(arguments[i], number_arguments, arguments, number_parameters, point, range, bounds);
				if (retval) {
					printf("ERROR creating search: %s -- [%s]\n", arguments[i], AS_MSG);
					continue;
				}
			}
			manage_search(arguments[i]);
			free(qualifier);
		}
	}
}

void asynchronous_search(int number_arguments, char** arguments, int number_parameters, double *point, double *range, BOUNDS* bounds) {
        int i, generate_result, insert_result, retval, position;
        MANAGED_SEARCH *ms = NULL;
        SEARCH_PARAMETERS *sp;

	init_search_manager(number_arguments, arguments);
        for (i = 0; i < number_arguments; i++) {
                if (!strcmp(arguments[i], "-s")) {
                        char *qualifier;
                        i++;

                        get_qualifier_from_name(arguments[i], &qualifier);
                        if (!search_exists(arguments[i])) {
                                ASYNCHRONOUS_SEARCH *as = get_registered_search(qualifier);
                                retval = as->create_search(arguments[i], number_arguments, arguments, number_parameters, point, range, bounds);
                                if (retval) {
                                        printf("ERROR creating search: %s -- [%s]\n", arguments[i], AS_MSG);
                                        return;
                                }
                        }
			position = manage_search(arguments[i]);
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
		sp->hostid = (int)(10 * drand48());
                insert_result = ms->search->insert_parameters(ms->search_name, ms->search_data, sp);
                printf("[%d] [%-150s] [%-20s] [%s] [%d]\n", i, AS_MSG, AS_INSERT_STR[insert_result], sp->metadata, sp->hostid);
                i++;
                if (i % 100 == 0) ms->search->checkpoint_search(ms->search_name, ms->search_data);
                free(sp->parameters);
                sp->number_parameters = number_parameters;
                sp->parameters = (double*)malloc(sizeof(double) * number_parameters);
        }
}
