#ifndef FGDO_SEARCH_H
#define FGDO_SEARCH_H

#include "../searches/search_parameters.h"

/********
	*	Generic search data.
 ********/

typedef struct search SEARCH;

typedef int (*generate_parameters_type)(SEARCH*, SEARCH_PARAMETERS**);
typedef int (*insert_parameters_type)(SEARCH*, SEARCH_PARAMETERS*);

struct search {
	char*				search_name;
	void*				search_data;
	generate_parameters_type	generate_parameters;
	insert_parameters_type		insert_parameters;
};

#endif
