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

#ifndef FGDO_BOINC_ADD_WU_H
#define FGDO_BOINC_ADD_WU_H

/********
	*	BOINC includes
 ********/
#include "boinc_db.h"
#include "sched_config.h"
#include "util.h"

/********
	*	FGDO includes
 ********/
#include "search_manager.h"
#include "../searches/search_parameters.h"


typedef struct workunit_info {
	int number_parameters;
	double rsc_fpops_est, rsc_fpops_bound, rsc_memory_bound, rsc_disk_bound;
	int delay_bound, min_quorum, target_nresults, max_error_results, max_total_results, max_success_results;
	int number_required_files;
	char *credit_str, *result_xml_path, *template_filename, *template_file;
	char** required_files;
} WORKUNIT_INFO;


int fwrite_workunit_info(FILE* out, WORKUNIT_INFO *wu_info);
int write_workunit_info(char* outfile, WORKUNIT_INFO *wu_info);

int fread_workunit_info(FILE* in, WORKUNIT_INFO **wu_info);
int read_workunit_info(char* infile, WORKUNIT_INFO **wu_info);

void add_workunit(SEARCH_PARAMETERS *parameters, WORKUNIT_INFO *wu_info, DB_APP app);

#endif
