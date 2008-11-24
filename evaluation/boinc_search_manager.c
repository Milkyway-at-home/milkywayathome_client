#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <ctime>
#include <vector>
#include <string.h>
#include <iostream>
#include <sstream>

#include "time.h"

/********
	*	BOINC includes
 ********/
#include "backend_lib.h"
#include "boinc_db.h"
#include "config.h"
#include "error_numbers.h"
#include "parse.h"
#include "str_util.h"
#include "util.h"
#include "validate_util.h"

#include "sched_msgs.h"
#include "sched_config.h"
#include "sched_util.h"

/********
	*	FGDO includes
 ********/
#include "search_manager.h"
#include "boinc_add_workunit.h"
#include "../searches/asynchronous_search.h"
#include "../searches/search_parameters.h"
#include "../util/settings.h"

#define		LOCKFILE	"assimilator.out"
#define		PIDFILE		"assimilator.pid"
#define		SLEEP_INTERVAL	10

#define config_dir "/export/www/boinc/milkyway"

using std::vector;
using std::string;

DB_APP		bsm_app;
bool		update_db = true;
bool		noinsert = false;

int		wu_id_modulus = 0;
int		wu_id_remainder = 0;

int		sleep_interval = SLEEP_INTERVAL;

bool		one_pass = false;
int		one_pass_N_WU = 0;
long		processed_wus = 0;
int		unsent_wu_buffer = 400;

WORKUNIT_INFO** workunit_info;

void init_boinc_search_manager(int argc, char** argv) {
	int i, retval;
	char buf[256];

	init_search_manager(argc, argv);
	check_stop_daemons();
	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-one_pass_N_WU")) {
			one_pass_N_WU = atoi(argv[++i]);
			one_pass = true;
		} else if (!strcmp(argv[i], "-one_pass")) {
			one_pass = true;
		} else if (!strcmp(argv[i], "-sleep_interval")) {
			sleep_interval = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "-d")) {
			log_messages.set_debug_level(atoi(argv[++i]));
		} else if (!strcmp(argv[i], "-app")) {
			strcpy(bsm_app.name, argv[++i]);
		} else if (!strcmp(argv[i], "-dont_update_db")) {
			/********
				*	This option is for testing your assimilator.  When set, it ensures that the assimilator does not actually modify
				*	the assimilate_state of the workunits, so you can run your assimilator over and over again without affecting
				*	your project.
			 ********/
			update_db = false;
		} else if (!strcmp(argv[i], "-noinsert")) {
			/********
				*	This option is also for testing and is used to prevent the inserting of results into the *backend*
				*	(as opposed to the boinc) DB.
			 ********/
			noinsert = true;
		} else if (!strcmp(argv[i], "-mod")) {
			wu_id_modulus   = atoi(argv[++i]);
			wu_id_remainder = atoi(argv[++i]);
		} else if (!strcmp(argv[i], "-unsent_wu_buffer")) {
			/********
				*	Generate more workunits if less than this number are available on the server.
			 ********/
			unsent_wu_buffer = atoi(argv[++i]);
		}
	}

	if (wu_id_modulus) {
		log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "Using mod'ed WU enumeration.  modulus = %d  remainder = %d\n", wu_id_modulus, wu_id_remainder);
	}

	retval = config.parse_file(config_dir);
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

	printf("looking up app\n");

	sprintf(buf, "where name='%s'", bsm_app.name);
	retval = bsm_app.lookup(buf);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "Can't find app\n");
		exit(1);
	}

	printf("installing stop signal handler\n");

	install_stop_signal_handler();

	printf("initializing workunit info\n");
	workunit_info = (WORKUNIT_INFO**)malloc(sizeof(WORKUNIT_INFO*) * number_searches);
	for (i = 0; i < number_searches; i++) {
		printf("\tfor search: %s\n", searches[i]->search_name);
		init_workunit_info(searches[i]->search_name, &(workunit_info[i]));
		printf("\tsuccess.\n");
	}
	printf("finished.\n");
}

int generate_workunits() {
	int i, j, generated, current, generation_rate;
	SEARCH_PARAMETERS *sp;

	generation_rate = get_generation_rate();
	current = 0;
        for (i = 0; i < number_searches; i++) {
                generated = (generation_rate - current) / (number_searches - i);
                for (j = 0; j < generated; j++) {
                        searches[i]->search->generate_parameters(searches[i]->search_name, searches[i]->search_data, &sp);
			current++;
			add_workunit(sp, workunit_info[i]);
			free_search_parameters(sp);
			free(sp);
                }
        }
	return current;
}

int insert_workunit(WORKUNIT& wu, vector<RESULT>& /*results*/, RESULT& canonical_result) {
	int success;
	string output_file_name;
	SEARCH_PARAMETERS* sp = (SEARCH_PARAMETERS*)malloc(sizeof(SEARCH_PARAMETERS));
	SCOPE_MSG_LOG scope_messages(log_messages, SCHED_MSG_LOG::MSG_NORMAL);

	if (wu.canonical_resultid) {
		log_messages.printf_multiline(SCHED_MSG_LOG::MSG_DEBUG, canonical_result.xml_doc_out, "[%s] canonical result", wu.name);

		get_output_file_path(canonical_result, output_file_name);
		success = read_search_parameters(output_file_name.c_str(), sp);
		if (success <= 0) scope_messages.printf("[%s] Error parsing result file: %s [%d]\n", wu.name, output_file_name.c_str(), success);

		scope_messages.printf("[%s] Assimilating, search_name: [%s], fitness: [%lf], metadata: [%s]\n", wu.name, sp->search_name, sp->fitness, sp->metadata);

		success = insert_search_parameters(sp);
		if (success <= 0) scope_messages.printf("[%s] Error inserting results to search: %s [%d]\n", wu.name, sp->search_name, success);

		free_search_parameters(sp);
		free(sp);
	} else {
		scope_messages.printf("[%s] No canonical result\n", wu.name);
	}

	if (wu.error_mask&WU_ERROR_COULDNT_SEND_RESULT)		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: couldn't send a result\n", wu.name);
	if (wu.error_mask&WU_ERROR_TOO_MANY_ERROR_RESULTS)	log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many error results\n", wu.name);
	if (wu.error_mask&WU_ERROR_TOO_MANY_TOTAL_RESULTS)	log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many total results\n", wu.name);
	if (wu.error_mask&WU_ERROR_TOO_MANY_SUCCESS_RESULTS)	log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many success results\n", wu.name);
	return 0;
}

void start_search_manager() {
	DB_WORKUNIT wu;
	DB_RESULT canonical_result, result;
	bool did_something = false;
	int retval, num_generated, unsent_wus, processed_wus, num_assimilated;
	time_t start_time, current_time;
	double wus_per_second;

	time(&start_time);
	log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Starting at: %d\n", (int)start_time);

	processed_wus = 0;
	while (!one_pass) {
		char buf[256], mod_clause[256];

		check_stop_daemons();

		if (wu_id_modulus)	sprintf(mod_clause, " and workunit.id %% %d = %d ", wu_id_modulus, wu_id_remainder);
		else			strcpy(mod_clause, "");
		sprintf(buf, "where appid=%d and assimilate_state=%d %s limit %d", bsm_app.id, ASSIMILATE_READY, mod_clause, one_pass_N_WU ? one_pass_N_WU : 1000);

		num_assimilated = 0;
		while (!wu.enumerate(buf)) {
			vector<RESULT> results;     // must be inside while()!
			/********
				*	for testing purposes, pretend we did nothing
			 ********/
			if (update_db) did_something = true;
			log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "[%s] assimilating boinc WU %d; state=%d\n", wu.name, wu.id, wu.assimilate_state);

			sprintf(buf, "where workunitid=%d", wu.id);
			while (!result.enumerate(buf)) {
				results.push_back(result);
				if (result.id == wu.canonical_resultid) canonical_result = result;
			}
			insert_workunit(wu, results, canonical_result);

			if (update_db) {
				sprintf(buf, "assimilate_state=%d, transition_time=%d", ASSIMILATE_DONE, (int)time(0));

				retval = wu.update_field(buf);
				if (retval) {
					log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] update failed: %d\n", wu.name, retval);
					exit(1);
				}
			}
			num_assimilated++;
		}
		if (did_something) boinc_db.commit_transaction();

		if (num_assimilated)  {
			log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Assimilated %d workunits.\n", num_assimilated);
			count_unsent_results(unsent_wus, 0);

			processed_wus += num_assimilated;
			time(&current_time);
			wus_per_second = (double)processed_wus/((double)current_time-(double)start_time);
			log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "wus/sec: %lf, unsent wus: %d\n", wus_per_second, unsent_wus);

			if (unsent_wus < unsent_wu_buffer) {
				num_generated = generate_workunits();
				log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Generated %d new workunits.\n", num_generated);
			}
		}
		if (!one_pass) sleep(sleep_interval);
	}
}
