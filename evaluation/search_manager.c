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

#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "search_manager.h"
#include "boinc_search_manager.h"

#include "../util/settings.h"
#include "../searches/search_parameters.h"
#include "../searches/asynchronous_search.h"
#include "evaluator.h"

double (*evaluate)(double*);

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
		} else if (!strcmp(argv[i], "-gen")) {
			generation_rate = atoi(argv[++i]);
		}
	}
}


/********
	*	Handle registered searches.
 ********/
int number_registered_searches = 0;
ASYNCHRONOUS_SEARCH **registered_searches;

void print_registered_searches() {
	int i;
//	fprintf(stderr, "registered searches: \n");
//	for (i = 0; i < number_registered_searches; i++) fprintf(stderr, "\t%s\n", registered_searches[i]->search_qualifier);
}

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

void register_search(ASYNCHRONOUS_SEARCH *as) {
	int i;
	int pos = get_registered_search_pos(as->search_qualifier);
	if (pos < 0) {
		/********
			*	Search is not known. Inorder position to put the search is -(pos + 1).
		 ********/
		pos = -(pos + 1);
		number_registered_searches++;
		if (number_registered_searches == 0) {
			registered_searches = (ASYNCHRONOUS_SEARCH**)malloc(sizeof(ASYNCHRONOUS_SEARCH*));
		} else {
			registered_searches = (ASYNCHRONOUS_SEARCH**)realloc(registered_searches, sizeof(ASYNCHRONOUS_SEARCH*) * number_registered_searches);
		}
		for (i = number_registered_searches-1; i > pos; i--) {
			registered_searches[i] = registered_searches[i-1];
		}
		registered_searches[pos] = as;
	} else {
		fprintf(stderr, "ERROR registering search %s, already known.\n", as->search_qualifier);
	}
//	fprintf(stderr, "registered searches: \n");
//	for (i = 0; i < number_registered_searches; i++) fprintf(stderr, "\t%s\n", registered_searches[i]->search_qualifier);
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
	MANAGED_SEARCH **temp;
	ASYNCHRONOUS_SEARCH *as;
	void *search_data;
	int search_pos, i, success, current;

	/********
		*	Check to see if the search exists.
	 ********/
	if (!search_exists(search_name)) return -1;

	/********
		*	Check to see if the search is already being managed.
	 ********/
	search_pos = get_search_pos(search_name);
	if (search_pos >= 0) {
		fprintf(stderr, "ERROR manging search %s, already being managed.\n", search_name);
		return -1;
	}
	search_pos = -(search_pos + 1);

	success = get_qualifier_from_name(search_name, &search_qualifier);
	if (success < 0) {
		fprintf(stderr, "ERROR managing search %s, has invalid qualifier.\n", search_name);
		return -1;
	}

	as = get_registered_search(search_qualifier);
	if (as == NULL) {
		fprintf(stderr, "ERROR managing search %s, unknown search: %s\n", search_name, search_qualifier);
		fprintf(stderr, "registered searches: \n");
		for (i = 0; i < number_registered_searches; i++) fprintf(stderr, "\t%s\n", registered_searches[i]->search_qualifier);
		return -1;
	}
	free(search_qualifier);

	as->read_search(search_name, &search_data);

	ms = (MANAGED_SEARCH*)malloc(sizeof(MANAGED_SEARCH));
	ms->completed = 0;
	ms->search_name = (char*)malloc(sizeof(char) * SEARCH_NAME_SIZE);
	strcpy(ms->search_name, search_name);
	ms->search_data = search_data;
	ms->search = as;

	number_searches++;
	temp = (MANAGED_SEARCH**)malloc(sizeof(MANAGED_SEARCH*) * number_searches);
	current = 0;
	for (i = 0; i < number_searches; i++) {
		if (i == search_pos) {
			temp[i] = ms;
		} else {
			temp[i] = searches[current];
			current++;
		}
	}
	free(searches);
	searches = temp;

//	printf("inserted search [%s] to position: %d\n", ms->search_name, search_pos);
	return search_pos;
}

int search_exists(char* search_name) {
	int success;
	struct stat buf;
	char directory[FILENAME_SIZE];
	sprintf(directory, "%s/%s", get_working_directory(), search_name);

	success = stat(directory, &buf);
	return (success == 0);
}

int generate_search_parameters(SEARCH_PARAMETERS **sp) {
	int generated, current, i, j;
	generated = 0;
	current = 0;
	for (i = 0; i < number_searches; i++) {
		generated = (generation_rate - current) / (number_searches - i);
		for (j = 0; j < generated; j++) {
			sprintf(sp[current]->search_name, "%s", searches[i]->search_name);
			searches[i]->search->generate_parameters(searches[i]->search_name, searches[i]->search_data, sp[current]);
			current++;
		}
	}
	return current;
}

int insert_search_parameters(SEARCH_PARAMETERS *sp) {
	MANAGED_SEARCH *ms = get_search(sp->search_name);

	if (ms == NULL) {
		int pos = manage_search(sp->search_name);
		if (pos < 0) return -1;
		ms = searches[pos];
	}
	ms->search->insert_parameters(ms->search_name, ms->search_data, sp);
	return 0;
}
