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

/*
 *	This creates a workunit from the specified parameters
 */
void add_workunit(SEARCH_PARAMETERS* parameters) {
	DB_WORKUNIT wu;
	char wu_name[1024], wu_file[1024], wu_path[1024];
	char star_file[1024], star_path[1024];
	char astronomy_file[1024], astronomy_path[1024];
	char credit_str[512];
	long current_time;
	const char *required_files[3];

	check_stop_daemons();

	current_time = time(NULL) + getpid(); 
	sprintf(wu_name, "%s_%d_%ld", parameters->search_name, generated_wus, current_time);
	sprintf(wu_file, "%s_parameters_%d_%ld", parameters->search_name, generated_wus, current_time);
	generated_wus++;

	if (sched_config.download_path(wu_file, wu_path)) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator couldnt open fgdo workunit file for write: %s\n", wu_path);
		return;
	}

	sprintf(astronomy_file, "%s/astronomy_parameters.txt", parameters->search_path);
	ASTRONOMY_PARAMETERS *ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	int retval = read_astronomy_parameters(astronomy_file, ap);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimlator could not open parameter file for read: %s\n", astronomy_file);
		return;
	}
	if (sched_config.download_path(astronomy_file, astronomy_path)) {
		/********
			*	Astronomy parameters not moved to the download diretory yet, do this.
		 ********/
		retval = write_astronomy_parameters(astronomy_path, ap);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open parameter file for write: %s\n", astronomy_path);
			return;
		}
	}
	sprintf(star_file, "%s/stars.txt", parameters->search_path);
	if (sched_config.download_path(star_file, star_path)) {
		/********
			*	Stars not moved to the download diretory yet, do this.
		 ********/
		STAR_POINTS *sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
		int retval = read_star_points(star_file, sp);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open star file for read: %s\n", star_file);
			return;
		}
		retval = write_star_points(star_path, sp);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "assimilator could not open star file for write: %s\n", star_path);
			return;
		}
		free_star_points(sp);
		free(sp);
	}


	FILE *out_p = fopen(wu_path, "w+");
        fprintf(out_p, "%s\n", parameters->search_name);

        fprintf(out_p, "%d:", parameters->number_parameters);
        for (int i = 0; i < parameters->number_parameters; i++) fprintf(out_p, " %lf", parameters->parameters[i]);

	fprintf(out_p, "\n%s\n", parameters->metadata);

	fflush(out_p);
	fclose(out_p);

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
	required_files[0] = wu_file;
	required_files[1] = astronomy_file;
	required_files[2] = star_file;

        double credit = ((double)ap->r_steps * (double)ap->mu_steps * (double)ap->nu_steps) / (350.0 * 800.0 * 80.0);
        if (ap->convolve > 0) {
                credit = 6.5 * (((double)ap->r_steps * (double)ap->mu_steps * (double)ap->nu_steps) / (175.0 * 400.0 * 20.0));
                credit += (credit / 6.0) * ((ap->convolve - 30)/15.0);
        }
        credit /= 1.6;

	sprintf(credit_str, "<credit>%lf</credit>", credit);
	create_work(wu, wu_template, "templates/a_result.xml", "../templates/a_result.xml", required_files, number_required_files, sched_config, "", credit_str);

	free_parameters(ap);
	free(ap);
}

int main(int argc, char** argv) {
	char working_directory[512], search_name[512], buf[512];
	int retval;
	/********
		*	Arguments
		*		-d <working_directory>
		*		-s <search_name>
		*		-wus <wus_to_generate>
	 ********/

	retval = config.parse_file("..");
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "Can't parse ../config.xml: %s\n", boincerror(retval));
		exit(1);
	}
	log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Starting\n");

	retval = boinc_db.open(config.db_name, config.db_host, config.db_user, config.db_passwd);
	if (retval) { 
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "Can't open DB\n");
		exit(1);
	}

	sprintf(buf, "where name='%s'", app.name);
	retval = app.lookup(buf);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "Can't find app\n");
		exit(1);
	}

	if ( read_file_malloc("../templates/a_wu2.xml", wu_template) ) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "can't read astronomy WU template\n");
		exit(1);
	}

	search_manager__register_search("nm", start_newton_method);
	search_manager__register_search("gs", start_genetic_search);
	search_manager__register_search("de", start_differential_evolution);
	search_manager__register_search("pso", start_particle_swarm);

	init_search_manager(working_directory, add_workunit);
	manage_search(search_name);
	generate_workunits(search_name, wus_to_generate);
}
