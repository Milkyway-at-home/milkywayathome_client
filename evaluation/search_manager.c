#include <string.h>
#include <stdlib.h>

#include "search_manager.h"

#include "../util/settings.h"
#include "../searches/search_parameters.h"
#include "../searches/asynchronous_search.h"


/********
	*	Initialization
 ********/

int generation_rate = 10;

int get_generation_rate() { return generation_rate; }

void init_search_manager(int argc, char **argv) {
	int i;

	for (i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-cwd")) {
			set_working_directory(argv[++i]);
			printf("set working directory: %s\n", get_working_directory());
		} else if (!strcmp(argv[i], "-gen")) {
			generation_rate = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "-s")) {
			if (manage_search(argv[++i]) < 0) {
				fprintf(stderr, "ERROR: unknown search %s specified.\n", argv[i]);
			}
		}
	}
}


/********
	*	Handle registered searches.
 ********/
int number_registered_searches = 0;
ASYNCHRONOUS_SEARCH **registered_searches;

int get_registered_search_pos(char* search_qualifier) {
	int cmp, i;
	for (i = 0; i < number_registered_searches; i++) {
		cmp = strcmp(search_qualifier, registered_searches[i]->search_qualifier);
		if (cmp == 0)		return i;
		else if (cmp < 0)	return -1 - i;
	}
	return -1 - i;
}

ASYNCHRONOUS_SEARCH* get_registered_search(char* search_qualifier) {
	int pos = get_registered_search_pos(search_qualifier);
	if (pos >= 0) return registered_searches[pos];
	else return NULL;
}

void register_search(ASYNCHRONOUS_SEARCH as) {
	int i;
	int pos = get_registered_search_pos(as.search_qualifier);
	if (pos < 0) {
		/********
			*	Search is not known. Inorder position to put the search is -(pos + 1).
		 ********/
		pos = -(pos + 1);
		number_registered_searches++;
		registered_searches = (ASYNCHRONOUS_SEARCH**)realloc(registered_searches, sizeof(ASYNCHRONOUS_SEARCH*) * number_registered_searches);
		for (i = number_registered_searches-1; i > pos; i--) {
			registered_searches[i] = registered_searches[i-1];
		}
		registered_searches[pos] = &as;
	} else {
		fprintf(stderr, "ERROR registering search %s, already known.\n", as.search_qualifier);
	}
}


/********
	*	Manage searches.
 ********/

int number_searches = 0;
MANAGED_SEARCH** searches;  

int get_search_pos(char* search_name) {
	int cmp, i;
	for (i = 0; i < number_searches; i++) {
		cmp = strcmp(search_name, searches[i]->search_name);
		if (cmp == 0)		return i;
		else if (cmp < 0)	return -1 - i;
	}
	return -1 - i;
}

MANAGED_SEARCH* get_search(char* search_name) {
	int pos = get_search_pos(search_name);
	if (pos >= 0) return searches[pos];
	else return NULL;
}

int get_qualifier_from_name(char* search_name, char** search_qualifier) {
	int i;
	(*search_qualifier) = (char*)malloc(sizeof(char) * SEARCH_QUALIFIER_SIZE);

	for (i = 0; i < (int)strlen(search_name); i++) {
		if (search_name[i] == '_') break;
	}
	if (i >= (int)strlen(search_name)) return -1;
	strncpy((*search_qualifier), search_name, i);
	(*search_qualifier)[i] = '\0';
	return i;
}

int manage_search(char* search_name) {
	char *search_qualifier;
	MANAGED_SEARCH *ms;
	ASYNCHRONOUS_SEARCH *as;
	void *search_data;
	int search_pos, i, success;
	/********
		*	Check to see if the search is already being managed.
	 ********/
	printf("getting search pos\n");
	search_pos = get_search_pos(search_name);
	printf("search pos: %d\n", search_pos);
	if (search_pos >= 0) {
		fprintf(stderr, "ERROR manging search %s, already being managed.\n", search_name);
		return -1;
	}
	search_pos = -(search_pos + 1);
	printf("new search pos: %d\n", search_pos);

	printf("getting qualifier\n");
	success = get_qualifier_from_name(search_name, &search_qualifier);
	printf("qualifier: %s\n", search_qualifier);
	if (success < 0) {
		fprintf(stderr, "ERROR managing search %s, has invalid qualifier.\n", search_name);
		return -1;
	}

	printf("getting registered search\n");
	as = get_registered_search(search_qualifier);
	if (as == NULL) {
		fprintf(stderr, "ERROR managing search %s, unknown search: %s\n", search_name, search_qualifier);
		return -1;
	}
	free(search_qualifier);

	printf("init search\n");
	printf("doing init\n");
	as->read_search(search_name, &search_data);
	printf("success\n");
	ms = (MANAGED_SEARCH*)malloc(sizeof(MANAGED_SEARCH));
	ms->search_name = (char*)malloc(sizeof(char) * SEARCH_NAME_SIZE);
	strcpy(ms->search_name, search_name);
	ms->search_data = search_data;
	ms->search = as;

	number_searches++;
	searches = (MANAGED_SEARCH**)realloc(searches, sizeof(MANAGED_SEARCH*) * number_searches);
	for (i = number_searches; i > search_pos; i--) {
		searches[i] = searches[i-1];
	}
	searches[search_pos] = ms;
	printf("inserted search in position: %d\n", search_pos);
	return 1;
}

int generate_search_parameters(SEARCH_PARAMETERS **sp) {
	int generated, current, i, j;
	generated = 0;
	current = 0;
	for (i = 0; i < number_searches; i++) {
		generated = (generation_rate - current) / (number_searches - i);
		for (j = 0; j < generated; j++) {
			searches[i]->search->generate_parameters(searches[i]->search_name, searches[i]->search_data, &(sp[current]));
			current++;
		}
	}
	return current;
}

int insert_search_parameters(SEARCH_PARAMETERS *sp) {
	MANAGED_SEARCH *ms = get_search(sp->search_name);
	if (ms != NULL) {
		ms->search->insert_parameters(ms->search_name, ms->search_data, sp);
		return 1;
	}
	return -1;
}
