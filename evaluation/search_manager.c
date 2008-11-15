#include <string.h>
#include <iostream>
#include <sstream>

#include "search_manager.h"

#include "../util/settings.h"
#include "../searches/search_parameters.h"
#include "../searches/search.h"

/********
	*	Handle registered searches.
 ********/
int number_registered_searches = 0;
REGISTERED_SEARCH **registered_searches;

int get_registered_search_pos(char* search_qualifier) {
	int cmp, i;
	for (i = 0; i < number_registered_searches; i++) {
		cmp = strcmp(search_qualifier, registered_searches[i]->search_qualifier);
		if (cmp == 0)		return i;
		else if (cmp < 0)	return -1 - i;
	}
	return -1 - i;
}

REGISTERED_SEARCH* get_registered_search(char* search_qualifier) {
	int pos = get_registered_search_pos(search_qualifier);
	if (pos >= 0) return registered_searches[pos];
	else return NULL;
}

void register_search(char* search_qualifier, init_search_type is) {
	int i;
	int pos = get_registered_search_pos(search_qualifier);
	if (pos < 0) {
		/********
			*	Search is not known. Inorder position to put the search is -(pos + 1).
		 ********/
		pos = -(pos + 1);
		number_registered_searches++;
		registered_searches = (REGISTERED_SEARCH**)realloc(registered_searches, sizeof(REGISTERED_SEARCH*) * number_registered_searches);
		for (i = number_registered_searches-1; i > pos; i--) {
			registered_searches[i] = registered_searches[i-1];
		}
		registered_searches[pos] = (REGISTERED_SEARCH*)malloc(sizeof(REGISTERED_SEARCH));
		registered_searches[pos]->search_qualifier = (char*)malloc(sizeof(char) * SEARCH_QUALIFIER_SIZE);
		strcpy(registered_searches[pos]->search_qualifier, search_qualifier);
		registered_searches[pos]->init_search = is;
	} else {
		fprintf(stderr, "ERROR registering search %s, already known.\n", search_qualifier);
	}
}


/********
	*	Manage searches.
 ********/

int number_searches = 0;
SEARCH** searches;  

int get_search_pos(char* search_name) {
	int cmp, i;
	for (i = 0; i < number_searches; i++) {
		cmp = strcmp(search_name, searches[i]->search_name);
		if (cmp == 0)		return i;
		else if (cmp < 0)	return -1 - i;
	}
	return -1 - i;
}

SEARCH* get_search(char* search_name) {
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

void manage_search(char* search_name) {
	char *search_qualifier;
	REGISTERED_SEARCH *rs;
	SEARCH *s;
	int search_pos, i, success;
	/********
		*	Check to see if the search is already being managed.
	 ********/
	search_pos = get_search_pos(search_name);
	if (search_pos >= 0) {
		fprintf(stderr, "ERROR manging search %s, already being managed.\n", search_name);
		return;
	}

	success = get_qualifier_from_name(search_name, &search_qualifier);
	if (success < 0) {
		fprintf(stderr, "ERROR managing search %s, has invalid qualifier.\n", search_name);
		return;
	}

	rs = get_registered_search(search_qualifier);
	if (rs == NULL) {
		fprintf(stderr, "ERROR managing search %s, unknown search: %s\n", search_name, search_qualifier);
		return;
	}

	s = (SEARCH*)malloc(sizeof(SEARCH));
	rs->init_search(search_name, s);
	number_searches++;
	searches = (SEARCH**)realloc(searches, sizeof(SEARCH*) * number_searches);
	for (i = number_searches; i > search_pos; i--) {
		searches[i] = searches[i-1];
	}
	searches[search_pos] = s;
}

int generate_search_parameters(int number_workunits, SEARCH_PARAMETERS **sp) {
	int generated, current, i, j;
	generated = 0;
	current = 0;
	for (i = 0; i < number_searches; i++) {
		generated = (number_workunits - current) / (number_searches - i);
		for (j = 0; j < generated; j++) {
			searches[i]->generate_parameters(searches[i], sp[current]);
			strcpy(sp[current]->search_name, searches[i]->search_name);
			current++;
		}
	}
	return current;
}

int insert_search_parameters(SEARCH_PARAMETERS *sp) {
	SEARCH *s = get_search(sp->search_name);
	if (s != NULL) {
		s->insert_parameters(s, sp);
		return 1;
	}
	return -1;
}
