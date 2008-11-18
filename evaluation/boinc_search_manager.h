#ifndef FGDO_BOINC_SEARCH_MANAGER_H
#define FGDO_BOINC_SEARCH_MANAGER_H

#include "../searches/search_parameters.h"

void init_boinc_search_manager(int argc, char** argv, void (*add_wu)(SEARCH_PARAMETERS*));
int generate_workunits();

#endif
