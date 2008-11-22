/********
	*	FGDO includes
 ********/
#include "../evaluation/search_manager.h"
#include "../searches/genetic_search.h"
#include "../searches/newton_method.h"
#include "../searches/differential_evolution.h"
#include "../searches/particle_swarm.h"
#include "../searches/search_parameters.h"

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
SCHED_CONFIG sched_config;
DB_APP app;
int generated_wus = 0;

#define config_dir "/export/www/boinc/milkyway"

/*
 *	This creates a workunit from the specified parameters
 */
void add_workunit(SEARCH_PARAMETERS* parameters) {
	DB_WORKUNIT wu;
	char wu_name[1024], wu_file[1024], wu_path[1024];
	char star_file[1024], star_path[1024];
	char astronomy_file[1024], astronomy_path[1024];
	char wu_astronomy_file[1024], wu_star_file[1024];
	char credit_str[512];
	time_t current_time;
	const char *required_files[3];
	int retval;

	printf("checking stop daemons\n");

	check_stop_daemons();

	printf("checked\n");

	time(&current_time);
	printf("got time: %ld\n", (long)current_time);
	printf("generated wus: %d\n", generated_wus);

	if (parameters == NULL) printf("ERROR! parameters = null!\n");
	if (parameters->search_name == NULL) printf("ERROR! search name == null!!!\n");
	printf("search name: %s\n", parameters->search_name);

	sprintf(wu_name, "%s_%d_%ld", parameters->search_name, generated_wus, (long)current_time);
	printf("made wu_name: %s\n", wu_name);

	sprintf(wu_file, "%s_search_parameters_%d_%ld", parameters->search_name, generated_wus, current_time);
	printf("made wu_file: %s\n", wu_file);

	sprintf(wu_astronomy_file, "%s_astronomy_parameters_%d_%ld", parameters->search_name, generated_wus, current_time);
	printf("made wu_astronomy_file: %s\n", wu_astronomy_file);

	sprintf(wu_star_file, "stars.txt");

	generated_wus++;


	if (!sched_config.download_path(wu_file, wu_path)) {
		printf("writing wu to path: %s\n", wu_path);
		retval = write_search_parameters(wu_path, parameters);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open search parameters file for write: %s\n", wu_path);
			return;
		}
	} else {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator couldnt open fgdo workunit file for write: %s\n", wu_path);
		return;
	}

	printf("reading astronomy file: %s/%s/astronomy_parameters.txt\n", get_working_directory(), parameters->search_name);

	sprintf(astronomy_file, "%s/%s/astronomy_parameters.txt", get_working_directory(), parameters->search_name);

	ASTRONOMY_PARAMETERS *ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	retval = read_astronomy_parameters(astronomy_file, ap);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimlator could not open parameter file for read: %s\n", astronomy_file);
		return;
	}
	if (!sched_config.download_path(wu_astronomy_file, astronomy_path)) {
		/********
			*	Astronomy parameters not moved to the download diretory yet, do this.
		 ********/
		printf("writing astronomy parameters to: %s\n", astronomy_path);
		retval = write_astronomy_parameters(astronomy_path, ap);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open parameter file for write: %s\n", astronomy_path);
			return;
		}
	}

	printf("reading star file: %s/%s/stars.txt\n", get_working_directory(), parameters->search_name);

	sprintf(star_file, "%s/%s/stars.txt", get_working_directory(), parameters->search_name);
	if (sched_config.download_path(wu_star_file, star_path)) {
		/********
			*	Stars not moved to the download diretory yet, do this.
		 ********/
		STAR_POINTS *sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
		int retval = read_star_points(star_file, sp);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open star file for read: %s\n", star_file);
			return;
		}
		printf("writing stars to: %s\n", star_path);
		retval = write_star_points(star_path, sp);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open star file for write: %s\n", star_path);
			return;
		}
		free_star_points(sp);
		free(sp);
	}

	wu.clear();
	wu.appid = app.id;
	strcpy(wu.name, wu_name);
	wu.rsc_fpops_est = 1e12;
	wu.rsc_fpops_bound = 1e14;
	wu.rsc_memory_bound = 5e8;
	wu.rsc_disk_bound = 15e6;		//15 MB
	wu.delay_bound = 60 * 60 * 24 * 5;	//5 days
	wu.min_quorum = 1;
	wu.target_nresults = 1;
	wu.max_error_results = 0;
	wu.max_total_results = 4;
	wu.max_success_results = 1;

	int number_required_files = 3;
	required_files[0] = wu_astronomy_file;
	required_files[1] = wu_star_file;
	required_files[2] = wu_file;

        double credit = ((double)ap->r_steps * (double)ap->mu_steps * (double)ap->nu_steps) / (350.0 * 800.0 * 80.0);
        if (ap->convolve > 0) {
                credit = 6.5 * (((double)ap->r_steps * (double)ap->mu_steps * (double)ap->nu_steps) / (175.0 * 400.0 * 20.0));
                credit += (credit / 6.0) * ((ap->convolve - 30)/15.0);
        }
        credit /= 1.6;

	sprintf(credit_str, "<credit>%lf</credit>", credit);

	char result_xml_path[1024];
	sprintf(result_xml_path, "%s/templates/a_result.xml", config_dir);

	create_work(wu, wu_template, "templates/a_result.xml", result_xml_path, required_files, number_required_files, sched_config, "", credit_str);

	free_parameters(ap);
	free(ap);
}

void init_add_workunit() {
	char buf[512];
	char wu_template_file[1024];
	int retval;

	retval = sched_config.parse_file(config_dir);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "Can't parse ../config.xml: %s\n", boincerror(retval));
		exit(1);
	}
	log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Starting\n");

	retval = boinc_db.open(sched_config.db_name, sched_config.db_host, sched_config.db_user, sched_config.db_passwd);
	if (retval) { 
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "Can't open DB\n");
		exit(1);
	}

	sprintf(buf, "where name='milkyway'");
	retval = app.lookup(buf);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "Can't find app: milkyway\n");
		exit(1);
	}

	sprintf(wu_template_file, "%s/templates/milkyway_wu.xml", config_dir);
	if ( read_file_malloc(wu_template_file, wu_template) ) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "can't read milkyway WU template\n");
		exit(1);
	}
}
