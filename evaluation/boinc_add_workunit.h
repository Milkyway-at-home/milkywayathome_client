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


typedef struct workunit_info {
	int number_parameters;
	double rsc_fpops_est, rsc_fpops_bound, rsc_memory_bound, rsc_disk_bound;
	int delay_bound, min_quorum, target_nresults, max_error_results, max_total_results, max_success_results;
	char* credit_str, *result_xml_path;
	char** required_files;
} WORKUNIT_INFO;


void fwrite_workunit_info(FILE* out, WORKUNIT_INFO *wu_info);

void add_workunit(SEARCH_PARAMETERS *parameters, WORKUNIT_INFO *wu_info);

void init_workunit_info(char* search_name, WORKUNIT_INFO **wu_info, DB_APP db_app);

#endif
