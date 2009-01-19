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
SEARCH_PARAMETERS **gen_sp, *insert_sp;

long checkpoint_time = 360;		//	1 hour

void update_workunit_info(int pos) {
	int i, current;
	WORKUNIT_INFO **temp_wu;
	SEARCH_PARAMETERS **temp_sp;

	temp_wu = (WORKUNIT_INFO**)malloc(sizeof(WORKUNIT_INFO*) * number_searches);
	temp_sp = (SEARCH_PARAMETERS**)malloc(sizeof(SEARCH_PARAMETERS*) * number_searches);
	current = 0;
	printf("inserting wu_info to pos %d of %d\n", pos, number_searches);
	for (i = 0; i < number_searches; i++) {
		if (i == pos) {
			char *workunit_info_file = (char*)malloc(sizeof(char) * FILENAME_SIZE);
			sprintf(workunit_info_file, "%s/%s/workunit_info", get_working_directory(), searches[pos]->search_name);
			read_workunit_info(workunit_info_file, &(temp_wu[pos]));
			init_search_parameters(&temp_sp[i], temp_wu[pos]->number_parameters);
			sprintf(temp_sp[i]->search_name, searches[i]->search_name);
		} else {
			temp_wu[i] = workunit_info[current];
			temp_sp[i] = gen_sp[current];
			current++;
		}
	}
	free(workunit_info);
	free(gen_sp);
	workunit_info = temp_wu;
	gen_sp = temp_sp;
}

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
		} else if (!strcmp(argv[i], "-cp_time")) {
			checkpoint_time = atoi(argv[++i]);
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

	sprintf(buf, "where name='%s'", bsm_app.name);
	retval = bsm_app.lookup(buf);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "Can't find app\n");
		exit(1);
	}

	install_stop_signal_handler();

	init_search_parameters(&insert_sp, 8);
}

int generate_workunits() {
	int i, j, generated, current, generation_rate, initial, result, number_completed, count_completed;
	SCOPE_MSG_LOG scope_messages(log_messages, SCHED_MSG_LOG::MSG_NORMAL);

	generation_rate = get_generation_rate();
	current = 0;
	number_completed = 0;
	count_completed = 0;
	for (i = 0; i < number_searches; i++) if (searches[i]->completed) number_completed++;

        for (i = 0; i < number_searches; i++) {
		if (searches[i]->completed) count_completed++;
                generated = (generation_rate - current) / ((number_searches - number_completed) - (i - count_completed));
		initial = current;
                for (j = 0; j < generated; j++) {
			result = searches[i]->search->generate_parameters(searches[i]->search_name, searches[i]->search_data, gen_sp[i]);
			if (result != AS_GEN_SUCCESS) {
				if (result == AS_GEN_OVER) searches[i]->completed = 1;
				scope_messages.printf("Not generating workunits: [%s]\n", AS_GEN_STR[result]);
				break;
			} 
			current++;
			add_workunit(gen_sp[i], workunit_info[i], bsm_app);
                }
		scope_messages.printf("[%s] Generated %d workunits.\n", searches[i]->search_name, (current-initial));
        }
	return current;
}

void generate_search_workunits(char *search_name) {
	manage_search(search_name);
	update_workunit_info(0);
	generate_workunits();
}

int insert_workunit(WORKUNIT& wu, vector<RESULT>& /*results*/, RESULT& canonical_result) {
	int result, hostid, userid, retval;
	string output_file_name;
	int sent_time, received_time, trip_time;
	double cpu_time;
	int v_num;
	char version[512];
	SCOPE_MSG_LOG scope_messages(log_messages, SCHED_MSG_LOG::MSG_NORMAL);
	DB_HOST host;

	if (wu.canonical_resultid) {
		MANAGED_SEARCH *ms;
		log_messages.printf_multiline(SCHED_MSG_LOG::MSG_DEBUG, canonical_result.xml_doc_out, "[%s] canonical result", wu.name);

		get_output_file_path(canonical_result, output_file_name);
		result = boinc_read_search_parameters2(output_file_name.c_str(), insert_sp, version);
		if (result != 0) {
			scope_messages.printf("[%-18s] could not read search parameters file: [%s]\n", wu.name, output_file_name.c_str());
			return result;
		}

		hostid = canonical_result.hostid;
		userid = canonical_result.userid;
		cpu_time = canonical_result.cpu_time;
		sent_time = canonical_result.sent_time;
		received_time = canonical_result.received_time;
		trip_time = received_time - sent_time;
		v_num = canonical_result.app_version_num;

		ms = get_search(insert_sp->search_name);

		retval = host.lookup_id(hostid);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] lookup of host %d failed %d\n", canonical_result.id, canonical_result.hostid, retval);
			return retval;   
		}

//		scope_messages.printf("[%-18s] OS[%s -- %s] CPU[%s -- %s]\n", insert_sp->search_name, host.os_name, host.os_version, host.p_vendor, host.p_model);
		if (ms == NULL) {
			int pos;
			pos = manage_search(insert_sp->search_name);
			if (pos < 0) {
		   		scope_messages.printf("[%-18s] [%-120s][%-25s] v[%-2d:%10s] trip/cpu[%*d/%*.2lf] u/h[%*d/%*d]\n", insert_sp->search_name, wu.name, "unknown search", v_num, version, 6, trip_time, 8, cpu_time, 6, userid, 6, hostid);
				return 1;
			}
			ms = searches[pos];
			update_workunit_info(pos);
		}

		sprintf(insert_sp->app_version, "%s", version);
		if (host.os_name[0] = 'M') sprintf(insert_sp->host_os, "Windows");
		else if (host.os_name[0] = 'D') sprintf(insert_sp->host_os, "Darwin");
		else if (host.os_name[0] = 'L') sprintf(insert_sp->host_os, "Linux");
		else sprintf(insert_sp->host_os, "?");
//		sprintf(insert_sp->host_os, "%s %s", host.os_name, host.os_version);

		printf("SET VERSION: %s, OS: %s\n", insert_sp->app_version, insert_sp->host_os);

		result = ms->search->insert_parameters(ms->search_name, ms->search_data, insert_sp);
		scope_messages.printf("[%-18s] [%-120s][%-25s] v[%-2d:%10s] [%-45s]\n", insert_sp->search_name, AS_MSG, AS_INSERT_STR[result], v_num, version, host.os_name);
//		scope_messages.printf("[%-18s] [%-120s][%-25s] v[%-2d:%10s] trip/cpu[%*d/%*.2lf] u/h[%*d/%*d]\n", insert_sp->search_name, AS_MSG, AS_INSERT_STR[result], v_num, version, 6, trip_time, 8, cpu_time, 6, userid, 6, hostid);
//		scope_messages.printf("[%-18s] [%-120s][%-25s] time[%*d/%*.2lf] u/h[%*d/%*d] v[%3d] OS[%s]\n", insert_sp->search_name, AS_MSG, AS_INSERT_STR[result], 6, trip_time, 8, cpu_time, 6, userid, 6, hostid, v_num, host.os_name);
		AS_MSG[0] = '\0';
	} else {
		scope_messages.printf("[%-18s] No canonical result\n", wu.name);
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
	int retval, num_generated, unsent_wus, processed_wus, num_assimilated, i;
	time_t start_time, current_time, last_checkpoint;
	double wus_per_second;

	time(&start_time);
	time(&last_checkpoint);
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
			count_unsent_results(unsent_wus, bsm_app.id);

			processed_wus += num_assimilated;
			time(&current_time);
			wus_per_second = (double)processed_wus/((double)current_time-(double)start_time);
			log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "[appid: %d] assimilated %d workunits, wus/sec: %lf, unsent wus: %d\n", bsm_app.id, num_assimilated, wus_per_second, unsent_wus);

			if (unsent_wus < unsent_wu_buffer) {
				log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Generating %d new workunits.\n", get_generation_rate());
				num_generated = generate_workunits();
				log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Generated %d new workunits.\n", num_generated);
			}

			if ((current_time - last_checkpoint) > checkpoint_time) {
				log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Checkpointing %d searches after %ld seconds.\n", number_searches, (current_time - last_checkpoint));
				{
					SCOPE_MSG_LOG scope_messages(log_messages, SCHED_MSG_LOG::MSG_NORMAL);
					for (i = 0; i < number_searches; i++) {
						retval = searches[i]->search->checkpoint_search(searches[i]->search_name, searches[i]->search_data);
						scope_messages.printf("[%-18s] checkpointed with result: [%s], msg: [%s]\n", searches[i]->search_name, AS_CP_STR[retval], AS_MSG);
						if (retval == AS_CP_OVER) searches[i]->completed = 1;
						AS_MSG[0] = '\0';
					}
				}
				log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Checkpointing completed.\n");
				last_checkpoint = current_time;
			}
		}
		if (!one_pass) sleep(sleep_interval);
	}
}
