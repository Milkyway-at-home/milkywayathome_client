#ifndef FGDO_SEARCH_H
#define FGDO_SEARCH_H

#include "../searches/search_parameters.h"

/********
	*	Generic search data.
 ********/

typedef int (*create_search_type)(char*, ...);
typedef int (*read_search_type)(char*, void**);
typedef int (*checkpoint_search_type)(char*, void*);
typedef int (*generate_parameters_type)(char*, void*, SEARCH_PARAMETERS**);
typedef int (*insert_parameters_type)(char*, void*, SEARCH_PARAMETERS*);

typedef struct asynchronous_search {
	char*				search_qualifier;

	create_search_type		create_search;
	read_search_type		read_search;
	checkpoint_search_type		checkpoint_search;
	generate_parameters_type	generate_parameters;
	insert_parameters_type		insert_parameters;
} ASYNCHRONOUS_SEARCH;

#endif
