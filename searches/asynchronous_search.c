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
        int i, generate_result, insert_result, retval;
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
                sp->parameters = (double*)malloc(sizeof(double) * number_parameters);
        }
}
