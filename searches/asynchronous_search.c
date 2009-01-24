#include "stdlib.h"
#include "stdio.h"
#include "string.h"

/********
	*	FGDO includes
 ********/
#include "asynchronous_search.h"
#include "../evaluation/search_manager.h"

const char *AS_GEN_STR[] = { "success", "search completed", "ERROR", "failure" };
const char *AS_INSERT_STR[] = { "success", "search completed", "fitness is NAN", "fitness is invalid", "parameters contain NAN", "parameters out of bounds", "ERROR", "out of range", "out of iteration", "bad metadata", "not unique"};
const char *AS_CP_STR[] = { "success", "search completed", "ERROR" };

char AS_MSG[1024] = "";
char AS_VERIFY_MSG[1024] = "";

void asynchronous_search__init(int number_arguments, char **arguments, int number_parameters, double *point, double *range, double *min_bound, double *max_bound) {
	int i, retval;

	for (i = 0; i < number_arguments; i++) {
		if (!strcmp(arguments[i], "-s")) {
			char *qualifier;
			i++;
			if (!search_exists(arguments[i])) {
				get_qualifier_from_name(arguments[i], &qualifier);
				retval = (get_registered_search(qualifier))->create_search(arguments[i], number_arguments, arguments, number_parameters, point, range, min_bound, max_bound);
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
