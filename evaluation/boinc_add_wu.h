#ifndef FGDO_BOINC_ADD_WU_H
#define FGDO_BOINC_ADD_WU_H

/********
	*	BOINC includes
 ********/
#include "boinc_db.h"
#include "sched_config.h"

/********
	*	FGDO includes
 ********/
#include "search_manager.h"
#include "../searches/search_parameters.h"

void init_workunit_generator(DB_APP db_app, SCHED_CONFIG sched_config);

void add_workunit(SEARCH_PARAMETERS *parameters);

#endif
