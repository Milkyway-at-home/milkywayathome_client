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
DB_APP app;
int generated_wus = 0;

#define config_dir "/export/www/boinc/milkyway"


/*
 *	This creates a workunit from the specified parameters
 */
void add_workunit(SEARCH_PARAMETERS* parameters, WORKUNIT_INFO *wu_info) {
	DB_WORKUNIT wu;
	char wu_name[1024], wu_file[1024], wu_path[1024];
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
	required_files[0] = wu_info->required_files[0];		//	wu_astronomy_file
	required_files[1] = wu_info->required_files[1];		//	wu_star_file
	required_files[2] = wu_file;

	create_work(wu, wu_template, "templates/a_result.xml", wu_info->result_xml_path, required_files, number_required_files, config, "", wu_info->credit_str);
}

void init_workunit_info(char* search_name, WORKUNIT_INFO **wu_info, DB_APP db_app) {
	char wu_template_file[1024];
	char star_file[1024], star_path[1024];
	char wu_astronomy_file[1024], wu_star_file[1024];
	char astronomy_file[1024], astronomy_path[1024];
	int retval, i;

	app = db_app;

	sprintf(wu_template_file, "%s/templates/milkyway_wu.xml", config_dir);
	if ( read_file_malloc(wu_template_file, wu_template) ) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "can't read milkyway WU template\n");
		exit(1);
	}

	sprintf(wu_astronomy_file, "%s_astronomy_parameters.txt", search_name);
	sprintf(astronomy_file, "%s/%s/astronomy_parameters.txt", get_working_directory(), search_name);
	ASTRONOMY_PARAMETERS *ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	retval = read_astronomy_parameters(astronomy_file, ap);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimlator could not open parameter file for read: %s\n", astronomy_file);
		exit(0);
	}
	if (!config.download_path(wu_astronomy_file, astronomy_path)) {
		/********
			*	Astronomy parameters not moved to the download diretory yet, do this.
		 ********/
		retval = write_astronomy_parameters(astronomy_path, ap);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open parameter file for write: %s\n", astronomy_path);
			exit(0);
		}
	}

	sprintf(wu_star_file, "%s_stars.txt", search_name);
	sprintf(star_file, "%s/%s/stars.txt", get_working_directory(), search_name);
	STAR_POINTS *sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
	retval = read_star_points(star_file, sp);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open star file for read: %s\n", star_file);
		exit(0);
	}
	if (!config.download_path(wu_star_file, star_path)) {
		/********
			*	Stars not moved to the download diretory yet, do this.
		 ********/
		retval = write_star_points(star_path, sp);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open star file for write: %s\n", star_path);
			exit(0);
		}
	}

	double calc_prob_count;
	double credit = ((double)ap->r_steps * (double)ap->mu_steps * (double)ap->nu_steps * (double)ap->convolve) + ((double)sp->number_stars * (double)ap->convolve);
	for (i = 0; i < ap->number_cuts; i++) {
		credit += (ap->r_cut[i][2] * ap->mu_cut[i][2] * ap->nu_cut[i][2]);
	}
	calc_prob_count = credit;
	credit /= 270000000.0;

	printf("awarded credit: %lf\n", credit);

	(*wu_info) = (WORKUNIT_INFO*)malloc(sizeof(WORKUNIT_INFO));
	(*wu_info)->credit_str = (char*)malloc(sizeof(char) * 1024);
	(*wu_info)->result_xml_path = (char*)malloc(sizeof(char) * 1024);

	sprintf((*wu_info)->credit_str, "<credit>%lf</credit>", credit);
	sprintf((*wu_info)->result_xml_path, "%s/templates/a_result.xml", config_dir);

	(*wu_info)->required_files = (char**)malloc(sizeof(char*) * 2);
	(*wu_info)->required_files[0] = (char*)malloc(sizeof(char) * 1024);
	(*wu_info)->required_files[1] = (char*)malloc(sizeof(char) * 1024);

	strcpy((*wu_info)->required_files[0], wu_astronomy_file);
	strcpy((*wu_info)->required_files[1], wu_star_file);

	(*wu_info)->rsc_fpops_est = calc_prob_count * (40 + (45 * ap->convolve) + ap->number_streams * (5 + ap->convolve * 30));
	printf("calculated fpops: %lf\n", (*wu_info)->rsc_fpops_est);

	(*wu_info)->rsc_fpops_bound = (*wu_info)->rsc_fpops_est * 100;
	(*wu_info)->rsc_memory_bound = 5e8;
	(*wu_info)->rsc_disk_bound = 15e6;			//15 MB
	(*wu_info)->delay_bound = 60 * 60 * 24 * 5;		//5 days
	(*wu_info)->min_quorum = 1;
	(*wu_info)->target_nresults = 1;
	(*wu_info)->max_error_results = 0;
	(*wu_info)->max_total_results = 4;
	(*wu_info)->max_success_results = 1;

	free_parameters(ap);
	free(ap);
	free_star_points(sp);
	free(sp);
}
