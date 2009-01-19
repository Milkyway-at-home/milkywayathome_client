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

/********
	*	FGDO includes
 ********/
#include "../evaluation/search_manager.h"
#include "../evaluation/boinc_add_workunit.h"
#include "../searches/search_parameters.h"
#include "../util/settings.h"

#include <string.h>
#include <iostream>
#include <sstream>

/********
	*	BOINC includes
 ********/
#include "config.h"
#include "boinc_db.h"
#include "backend_lib.h"
#include "parse.h"
#include "util.h"
#include "str_util.h"

#include "sched_config.h"
#include "sched_util.h"
#include "sched_msgs.h"

/********
	*	ASTRONOMY includes
 ********/
#include "parameters.h"
#include "star_points.h"

using std::string;
using std::stringstream;

/*
 *	The following are needed to create workunits
 */
char* wu_template;
int generated_wus = 0;

#define config_dir "/export/www/boinc/milkyway"


/*
 *	This creates a workunit from the specified parameters
 */
void add_workunit(SEARCH_PARAMETERS* parameters, WORKUNIT_INFO *wu_info, DB_APP app) {
	DB_WORKUNIT wu;
	char wu_name[1024], wu_file[1024], wu_path[1024], ap_file[1024], sp_file[1024];
	time_t current_time;
	const char *required_files[3];
	int retval;

	check_stop_daemons();

	time(&current_time);
	if (parameters == NULL) printf("ERROR! parameters = null!\n");
	if (parameters->search_name == NULL) printf("ERROR! search name == null!!!\n");

	sprintf(wu_name, "%s_%d_%ld", parameters->search_name, generated_wus, (long)current_time);
	sprintf(wu_file, "%s_search_parameters_%d_%ld", parameters->search_name, generated_wus, current_time);
	generated_wus++;

	if (!config.download_path(wu_file, wu_path)) {
//		retval = fwrite_search_parameters(stdout, parameters);
		retval = write_search_parameters(wu_path, parameters);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open search parameters file for write: %s\n", wu_path);
			return;
		}
	} else {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator couldnt open fgdo workunit file for write: %s\n", wu_path);
		return;
	}

	wu.clear();
	wu.appid = app.id;
	strcpy(wu.name, wu_name);

	wu.rsc_fpops_est = wu_info->rsc_fpops_est;
	wu.rsc_fpops_bound = wu_info->rsc_fpops_bound;
	wu.rsc_memory_bound = wu_info->rsc_memory_bound;
	wu.rsc_disk_bound = wu_info->rsc_disk_bound;
	wu.delay_bound = wu_info->delay_bound;
	wu.min_quorum = wu_info->min_quorum;
	wu.target_nresults = wu_info->target_nresults;
	wu.max_error_results = wu_info->max_error_results;
	wu.max_total_results = wu_info->max_total_results;
	wu.max_success_results = wu_info->max_success_results;

	int number_required_files = 3;
	sprintf(ap_file, "%s", wu_info->required_files[0]);		//	wu_astronomy_file
	sprintf(sp_file, "%s", wu_info->required_files[1]);		//	wu_star_file

	required_files[0] = ap_file;
	required_files[1] = sp_file;
	required_files[2] = wu_file;

//	printf("required_files[0]: %s\n", required_files[0]);
//	printf("required_files[1]: %s\n", required_files[1]);
//	printf("required_files[2]: %s\n", required_files[2]);
//	printf("wu_info->template_file: %s\n", wu_info->template_file);
//	printf("wu_info->result_xml_path: %s\n", wu_info->result_xml_path);
//	printf("wu_info->credit_str: %s\n", wu_info->credit_str);


	create_work(wu, wu_info->template_file, "templates/a_result.xml", wu_info->result_xml_path, required_files, number_required_files, config, "", wu_info->credit_str);
}
