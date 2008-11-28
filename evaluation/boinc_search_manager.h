#ifndef FGDO_BOINC_SEARCH_MANAGER_H
#define FGDO_BOINC_SEARCH_MANAGER_H

#include "../searches/search_parameters.h"

void print_wu_info(int count);

void init_boinc_search_manager(int argc, char** argv);
int generate_workunits();

#endif
