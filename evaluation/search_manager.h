#ifndef FGDO_SEARCH_MANAGER_H
#define FGDO_SEARCH_MANAGER_H

#include "../searches/asynchronous_search.h"
#include "../searches/search_parameters.h"


/********
	*	Registry for known searches.
 ********/

typedef struct managed_search {
	int			completed;
	char*			search_name;
	void*			search_data;
	ASYNCHRONOUS_SEARCH*	search;
} MANAGED_SEARCH;


extern int number_registered_searches;
extern ASYNCHRONOUS_SEARCH **registered_searches;

extern int number_searches;
extern MANAGED_SEARCH** searches;


/********
	*	Search manager functions.
 ********/
int get_generation_rate();

void init_search_manager(int argc, char** argv);
void start_search_manager();

void register_search(ASYNCHRONOUS_SEARCH as);
int manage_search(char* search_name);
int search_exists(char* search_name);

MANAGED_SEARCH* get_search(char* search_name);

int generate_search_parameters(SEARCH_PARAMETERS **sp);
int insert_search_parameters(SEARCH_PARAMETERS *sp);

int get_qualifier_from_name(char* search_name, char **search_qualifier);

#endif
