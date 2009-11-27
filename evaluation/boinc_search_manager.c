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
#include "../searches/search_log.h"
#include "../searches/search_parameters.h"
#include "../searches/result.h"
#include "../util/settings.h"


#define		LOCKFILE	"assimilator.out"
#define		PIDFILE		"assimilator.pid"
#define		SLEEP_INTERVAL	10

#define config_dir "/export/www/boinc/milkyway"

using std::vector;
using std::string;

#define SELECT_LIMIT    1000
#define SLEEP_PERIOD    5

typedef enum { NEVER, DELAYED, IMMEDIATE, NO_CHANGE } TRANSITION_TIME;

DB_APP		bsm_app;
bool		update_db = true;
bool		noinsert = false;

int		wu_id_modulus = 0;
int		wu_id_remainder = 0;

int		sleep_interval = SLEEP_INTERVAL;

bool		one_pass = false;
int		one_pass_N_WU = 0;
long		validated_wus = 0;
long		assimilated_wus = 0;
int		unsent_wu_buffer = 600;
double		max_credit_per_cpu_second = 0;
double		max_granted_credit = 0;  
double		max_claimed_credit = 0;
bool		grant_claimed_credit = false;
bool		update_credited_job = false;
bool		credit_from_wu = false;


WORKUNIT_INFO** workunit_info;
SEARCH_PARAMETERS **gen_sp, *insert_sp;

long checkpoint_time = 360;		//	1 hour

void print_message(char *search_name, char *as_msg, const char *as_result, char *verify_msg, int v_num, char *version, char *host_os, double credit, RESULT& result) {
	SCOPE_MSG_LOG scope_messages(log_messages, SCHED_MSG_LOG::MSG_NORMAL);
	int trip_time = result.received_time - result.sent_time;
	if (trip_time < 0) trip_time = 0;

	//scope_messages.printf("[%-13s] [%-87s][%-17s][%-8s] v[%-37s][%-7s] c[%*.5lf], t[%*d/%*.2lf] h[%*d]\n", search_name, as_msg, as_result, verify_msg, version, host_os, 9, credit, 6, trip_time, 8, result.cpu_time, 6, result.hostid);
	scope_messages.printf("[%-20s] [%-87s][%-17s][%-8s] v[%-37s][%-7s] c[%*.3lf/%*.3lf], t[%*d/%*.2lf], u[%*d], h[%*d]\n", search_name, as_msg, as_result, verify_msg, version, host_os, 6, credit, 6, result.claimed_credit, 6, trip_time, 8, result.cpu_time, 6, result.userid, 6, result.hostid);
}

void update_workunit_info(int pos) {
	int i, current;
	WORKUNIT_INFO **temp_wu;
	SEARCH_PARAMETERS **temp_sp;

	temp_wu = (WORKUNIT_INFO**)malloc(sizeof(WORKUNIT_INFO*) * number_searches);
	temp_sp = (SEARCH_PARAMETERS**)malloc(sizeof(SEARCH_PARAMETERS*) * number_searches);
	current = 0;
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

void get_search_name(const char* wu_name, char *search_name) {
	int i, found_underscore;
	found_underscore = 0;
	memset(search_name, '\0', 64);
	for (i = strlen(wu_name)-1; i >= 0; i--) {
		if (found_underscore >= 2) {
			search_name[i] = wu_name[i];
		} else if (wu_name[i] == '_') found_underscore++;
	}
}

MANAGED_SEARCH* get_search_from_wu_name(const char* wu_name) {
	MANAGED_SEARCH *ms;
	char search_name[64];
	get_search_name(wu_name, search_name);
	ms = get_search(search_name);
	if (ms == NULL) {
		int pos;
		pos = manage_search(search_name);
		if (pos < 0) return NULL;
		ms = searches[pos];
		update_workunit_info(pos);
	}
        return ms;        
}

void update_error_rate(DB_HOST& host, bool valid) {
	if (host.error_rate > 1) host.error_rate = 1;
	if (host.error_rate <= 0) host.error_rate = 0.1;

	host.error_rate *= 0.95;
	if (!valid) host.error_rate += 0.05;
}


// Here when a result has been validated and its granted_credit has been set.
// Grant credit to host, user and team, and update host error rate.
int is_valid(RESULT& result, WORKUNIT& wu) {
	DB_USER user;
	DB_HOST host;
	DB_TEAM team;
	DB_CREDITED_JOB credited_job;
	int retval;
	char buf[256];

	retval = host.lookup_id(result.hostid);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] lookup of host %d failed %d\n", result.id, result.hostid, retval);
		return retval;
	}
	retval = user.lookup_id(host.userid);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] lookup of user %d failed %d\n", result.id, host.userid, retval);
		return retval;
	}

	update_average(result.sent_time, result.granted_credit, CREDIT_HALF_LIFE, user.expavg_credit, user.expavg_time);
	sprintf(buf, "total_credit=total_credit+%f, expavg_credit=%f, expavg_time=%f", result.granted_credit,  user.expavg_credit, user.expavg_time);

	retval = user.update_field(buf);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] update of user %d failed %d\n", result.id, host.userid, retval);
	}

	update_average(result.sent_time, result.granted_credit, CREDIT_HALF_LIFE, host.expavg_credit, host.expavg_time);

	double turnaround = result.received_time - result.sent_time;
	compute_avg_turnaround(host, turnaround);

	// compute new credit per CPU time
	retval = update_credit_per_cpu_sec(result.granted_credit, result.cpu_time, host.credit_per_cpu_sec);
	if (retval) {
		//log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d][HOST#%d] claimed too much credit (%f) in too little CPU time (%f)\n", result.id, result.hostid, result.granted_credit, result.cpu_time);
	}

	double old_error_rate = host.error_rate;
	update_error_rate(host, true);
	sprintf(buf, "total_credit=total_credit+%f, expavg_credit=%f, expavg_time=%f, avg_turnaround=%f, credit_per_cpu_sec=%f, error_rate=%f", result.granted_credit, host.expavg_credit, host.expavg_time, host.avg_turnaround, host.credit_per_cpu_sec, host.error_rate);

	retval = host.update_field(buf);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] update of host %d failed %d\n", result.id, result.hostid, retval);
	}
	log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "[HOST#%d] error rate %f->%f\n", host.id, old_error_rate, host.error_rate);

	if (user.teamid) {
		retval = team.lookup_id(user.teamid);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] lookup of team %d failed %d\n", result.id, user.teamid, retval);
			return retval;
		}
		update_average(result.sent_time, result.granted_credit, CREDIT_HALF_LIFE, team.expavg_credit, team.expavg_time);
		sprintf(buf, "total_credit=total_credit+%f, expavg_credit=%f, expavg_time=%f", result.granted_credit,  team.expavg_credit, team.expavg_time);
		retval = team.update_field(buf);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] update of team %d failed %d\n", result.id, team.id, retval);
		}
	}

	if (update_credited_job) {
		credited_job.userid = user.id;
		credited_job.workunitid = long(wu.opaque);
		retval = credited_job.insert();
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] Warning: credited_job insert failed (userid: %d workunit: %f err: %d)\n", result.id, user.id, wu.opaque, retval);
		} else {
			log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "[RESULT#%d %s] added credited_job record [WU#%d OPAQUE#%f USER#%d]\n", result.id, result.name, wu.id, wu.opaque, user.id);
		}
	}
	return 0;
}

int is_invalid(RESULT& result) {
	char buf[256];
	int retval;
	DB_HOST host;

	retval = host.lookup_id(result.hostid);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] lookup of host %d failed %d\n", result.id, result.hostid, retval);
		return retval;
	}
	double old_error_rate = host.error_rate;
	update_error_rate(host, false);
	sprintf(buf, "error_rate=%f", host.error_rate);
	retval = host.update_field(buf);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] update of host %d failed %d\n", result.id, result.hostid, retval);
		return retval;
	}
	log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "[HOST#%d] invalid result; error rate %f->%f\n", host.id, old_error_rate, host.error_rate);
	return 0;
}

double update_workunit(DB_VALIDATOR_ITEM_SET& validator, int valid_state, RESULT& result, WORKUNIT &wu) {
	SCOPE_MSG_LOG scope_messages(log_messages, SCHED_MSG_LOG::MSG_NORMAL);
	int retval, x;
	bool assimilate = false;
	vector<RESULT> results;
	double credit = 0;
	TRANSITION_TIME transition_time = NO_CHANGE;

	switch (valid_state) {
		case AS_VERIFY_VALID:
			results.push_back(result);
			// grant credit for valid results
			credit = get_credit_from_wu(wu, results);

			result.granted_credit = grant_claimed_credit ? result.claimed_credit : credit;
			if (max_granted_credit && result.granted_credit > max_granted_credit) result.granted_credit = max_granted_credit;
			//if (result.claimed_credit < credit) {
			//	printf("claimed [%lf] < granted [%lf]!\n", result.claimed_credit, credit);
			//	credit = result.claimed_credit;
			//}

//			if (cpu_time == trip_time) cpu_time = 0.01 * trip_time;
//			else if (cpu_time > 0.35 * trip_time) cpu_time = 0.15 * trip_time;

			//if (max_credit_per_cpu_second && (result.granted_credit / cpu_time) > max_credit_per_cpu_second) result.granted_credit = cpu_time * max_credit_per_cpu_second;
			//if (cpu_time > 5 && cpu_time < 30) result.granted_credit *= 3.0;

			credit = result.granted_credit;
			wu.canonical_credit = credit;
			wu.canonical_resultid = result.id;

			result.validate_state = VALIDATE_STATE_VALID;
			if (update_db) {
				retval = is_valid(result, wu);
				if (retval) log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "[RESULT#%d %s] is_valid() failed: %d\n", result.id, result.name, retval);
			}
//			log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "[RESULT#%d %s] Granted %f credit to valid result [HOST#%d]\n", result.id, result.name, result.granted_credit, result.hostid);
			assimilate = true;
			transition_time = IMMEDIATE;
			break;
		case AS_VERIFY_INVALID:
			result.validate_state = VALIDATE_STATE_INVALID;
			if (update_db) is_invalid(result);
			assimilate = true;
			transition_time = IMMEDIATE;
			break;
		case AS_VERIFY_IN_PROGRESS:
			result.validate_state = VALIDATE_STATE_INCONCLUSIVE;
			wu.need_validate = 1;
			transition_time = DELAYED;
			break;
	}

	if (update_db && assimilate) {
		DB_WORKUNIT db_wu;
//		log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "[RESULT#%d %s] granted_credit %f\n",  result.id, result.name, result.granted_credit);

		retval = validator.update_result(result);
		if (retval) log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d %s] Can't update result: %d\n", result.id, result.name, retval);

		wu.assimilate_state = ASSIMILATE_READY;
	}
	switch (transition_time) {
		case IMMEDIATE:
			wu.transition_time = time(0);
			break;
		case DELAYED:
			x = time(0) + 6*3600;
			if (x < wu.transition_time) wu.transition_time = x;
			break;
		case NEVER:
			wu.transition_time = INT_MAX;
			break;
		case NO_CHANGE:   
			break;
	}
	validator.update_workunit(wu);

	return credit;
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

//			for (k = 0; k < gen_sp[i]->number_parameters; k++) {
//				if (isnan(gen_sp[i]
//			}

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


int insert_workunit(DB_VALIDATOR_ITEM_SET& validator, std::vector<VALIDATOR_ITEM>& items) { 
	SCOPE_MSG_LOG scope_messages(log_messages, SCHED_MSG_LOG::MSG_NORMAL);
	int retval = 0, valid_state;
	double credit = 0, version = 0.0;
	unsigned int i;
	int j;
	char search_name[64];
	DB_HOST host;
	bool has_valid = false;

	WORKUNIT& wu = items[0].wu;

	for (i = 0; i < items.size(); i++) {
		RESULT& result = items[i].res;
		MANAGED_SEARCH *ms;
		string output_file_name;

		get_search_name(wu.name, search_name);

		if ((result.server_state != RESULT_SERVER_STATE_OVER) || (result.outcome != RESULT_OUTCOME_SUCCESS)) {
			print_message(search_name, "RESULT NOT COMPLETED", "", "in prog", result.app_version_num, "?", "?", credit, result);
			continue;
		}

		if (result.validate_state != VALIDATE_STATE_INIT) {
			print_message(search_name, "RESULT ALREADY VALIDATED", "", wu.name, result.app_version_num, "?", "?", credit, result);
			has_valid = true;
			continue;
		}

		/********
			*	Read the result file
		 ********/
		get_output_file_path(result, output_file_name);
		const char *output_string = output_file_name.c_str();
		retval = read_cpu_result__realloc(output_string, insert_sp->search_name, &(insert_sp->number_parameters), &(insert_sp->parameters), &(insert_sp->fitness), insert_sp->metadata, insert_sp->app_version);
		if (retval) {
			credit = update_workunit(validator, AS_VERIFY_INVALID, result, wu);
			print_message(search_name, "error reading result", "", "invalid", result.app_version_num, insert_sp->app_version, "?", credit, result);
			continue;
		}

		if (isnan(insert_sp->fitness) || insert_sp->fitness > -2) {
			credit = update_workunit(validator, AS_VERIFY_INVALID, result, wu);
			print_message(search_name, "invalid fitness (nan or > -2)", "", "invalid", result.app_version_num, insert_sp->app_version, "?", credit, result);
			continue;
		}

		insert_sp->hostid = result.hostid;
		retval = host.lookup_id(result.hostid);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d] lookup of host %d failed %d\n", result.id, result.hostid, retval);
			return retval;
		}       

		if (host.os_name[0] == 'M') sprintf(insert_sp->host_os, "Windows");
		else if (host.os_name[0] == 'D') sprintf(insert_sp->host_os, "Darwin");
		else if (host.os_name[0] == 'L') sprintf(insert_sp->host_os, "Linux");
		else sprintf(insert_sp->host_os, "?");

		if (insert_sp->app_version[0] == '?' || !strcmp(insert_sp->app_version, "#IND00000000000") || (insert_sp->app_version[0] == 't' && insert_sp->app_version[1] == 'i' && insert_sp->app_version[2] == 'm' && insert_sp->app_version[3] == 'e')) {
			credit = update_workunit(validator, AS_VERIFY_INVALID, result, wu);
			print_message(search_name, "invalid app_version", "", "invalid", result.app_version_num, insert_sp->app_version, "?", credit, result);
			continue;
		}
		for (j = strlen(insert_sp->app_version) - 1; j >= 0; j--) {
			if (insert_sp->app_version[j] == ':') break;
		}
		if ( !(j > 0 && sscanf(&(insert_sp->app_version[j]), ": %lf", &version)) ) {
			version = 0.0;
		}

		if (version < 0.16) {
			credit = update_workunit(validator, AS_VERIFY_INVALID, result, wu);
			print_message(search_name, "invalid app_version number", "", "invalid", result.app_version_num, insert_sp->app_version, insert_sp->host_os, credit, result);
			continue;
		}

		if (version == 0.16 && strlen(insert_sp->app_version) >= 11 && !strncmp(insert_sp->app_version, "mindc_linux", 11)) {
			credit = update_workunit(validator, AS_VERIFY_VALID, result, wu);
			print_message(search_name, "invalid app_version: mindc_linux: 0.16", "", "valid", result.app_version_num, insert_sp->app_version, insert_sp->host_os, credit, result);
			continue;
		}

		if ((NULL == strstr(insert_sp->app_version, "gpu") && NULL == strstr(insert_sp->app_version, "GPU")) && result.cpu_time < 600) {
			credit = update_workunit(validator, AS_VERIFY_INVALID, result, wu);
			print_message(search_name, "invalid CPU time", "", "invalid", result.app_version_num, insert_sp->app_version, "?", credit, result);
			continue;
//		} else if (result.cpu_time < 10) {
//			credit = update_workunit(validator, AS_VERIFY_INVALID, result, wu);
//			print_message(search_name, "invalid GPU time", "", "invalid", result.app_version_num, insert_sp->app_version, "?", credit, result);
//			continue;
		}

		ms = get_search_from_wu_name(wu.name);
		if (ms == NULL) {
			credit = update_workunit(validator, AS_VERIFY_VALID, result, wu);
			print_message(search_name, "completed search", "", "valid", result.app_version_num, insert_sp->app_version, insert_sp->host_os, credit, result);
			has_valid = true;
			continue;
		}

		if (ms->search == NULL) printf("search == NULL\n");
		if (ms->search_name == NULL) printf("search_name == NULL\n");
		if (ms->search_data == NULL) printf("search_data == NULL\n");
	

		FILE *error_file;

		retval = ms->search->insert_parameters(ms->search_name, ms->search_data, insert_sp);

		switch(retval) {
			case AS_INSERT_SUCCESS:
				valid_state = AS_VERIFY_IN_PROGRESS;
			break;
			case AS_INSERT_OVER:
			case AS_INSERT_OUT_OF_RANGE:
			case AS_INSERT_OUT_OF_ITERATION:
			case AS_INSERT_BAD_METADATA:
			case AS_INSERT_NOT_UNIQUE:
				valid_state = AS_VERIFY_VALID;
			break;
			case AS_INSERT_FITNESS_INVALID:		// NEED TO MOVE THIS TO INVALID
				//compare to previous iteration values
				error_file = error_log_open(ms->search_name);
				fprintf(error_file, "[%s] [%s] [invalid fitness]\n", insert_sp->app_version, AS_MSG);
				fwrite_search_parameters(error_file, insert_sp);
				fclose(error_file);
				valid_state = AS_VERIFY_VALID;
			break;
			case AS_INSERT_FITNESS_NAN:
			case AS_INSERT_PARAMETERS_NAN:
			case AS_INSERT_OUT_OF_BOUNDS:
			case AS_INSERT_ERROR:
			case AS_INSERT_OUTLIER:
				valid_state = AS_VERIFY_INVALID;
			break;
			default:
				valid_state = AS_VERIFY_VALID;
			break;
		}

		if (valid_state == AS_VERIFY_IN_PROGRESS) {
//			transition_time = DELAYED;
//			wu.need_validate = 1;
			valid_state = AS_VERIFY_VALID;
			sprintf(AS_VERIFY_MSG, "valid");
			has_valid = true;
		} else if (valid_state == AS_VERIFY_VALID) {
			sprintf(AS_VERIFY_MSG, "valid");
			has_valid = true;
		} else {
			sprintf(AS_VERIFY_MSG, "invalid");
			has_valid = true;
		}
		credit = update_workunit(validator, valid_state, result, wu);
		print_message(ms->search_name, AS_MSG, AS_INSERT_STR[retval], AS_VERIFY_MSG, result.app_version_num, insert_sp->app_version, insert_sp->host_os, credit, result);
		AS_VERIFY_MSG[0] = '\0';
	}

	if (has_valid) {
		/********
			*	Already have a result, don't sent any unsent results.
		 ********/
		for (i = 0; i < items.size(); i++) {
			RESULT& result = items[i].res;
			if (result.server_state != RESULT_SERVER_STATE_UNSENT && result.outcome != RESULT_OUTCOME_NO_REPLY && result.outcome != RESULT_OUTCOME_DIDNT_NEED) {
				wu.canonical_resultid = result.id;
				continue;
			}
			result.server_state = RESULT_SERVER_STATE_OVER;
			result.outcome = RESULT_OUTCOME_DIDNT_NEED;
			retval = validator.update_result(result);
			if (retval) {
				log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d %s] result.update() failed: %d\n", result.id, result.name, retval);
			}
		}
		wu.need_validate = 0;
		wu.assimilate_state = ASSIMILATE_DONE;
		//wu.transition_time = (int)time(0);
		retval = validator.update_workunit(wu);
		if (retval) {            
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[WU#%d %s] update_workunit() failed: %d; exiting\n", wu.id, wu.name, retval);
			return retval;
		}
	}

        if (wu.error_mask&WU_ERROR_COULDNT_SEND_RESULT)         log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: couldn't send a result\n", wu.name);
        if (wu.error_mask&WU_ERROR_TOO_MANY_ERROR_RESULTS)      log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many error results\n", wu.name);
        if (wu.error_mask&WU_ERROR_TOO_MANY_TOTAL_RESULTS)      log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many total results\n", wu.name);
        if (wu.error_mask&WU_ERROR_TOO_MANY_SUCCESS_RESULTS)    log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many success results\n", wu.name);

	return 0;
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
		} else if (!strcmp(argv[i], "-max_credit_per_cpu_second")) {
			max_credit_per_cpu_second = atof(argv[++i]);
		} else if (!strcmp(argv[i], "-max_granted_credit")) {
			max_granted_credit = atof(argv[++i]);
		} else if (!strcmp(argv[i], "-max_claimed_credit")) {
			max_claimed_credit = atof(argv[++i]);
		} else if (!strcmp(argv[i], "-grant_claimed_credit")) {
			grant_claimed_credit = true;
		} else if (!strcmp(argv[i], "-update_credited_job")) {
			update_credited_job = true;
		} else if (!strcmp(argv[i], "-credit_from_wu")) {
			credit_from_wu = true;
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

void start_search_manager() {
	DB_WORKUNIT wu;
	DB_RESULT result;
	bool did_something = false;
	int retval, num_generated, unsent_wus, num_validated, num_assimilated;
	time_t start_time, current_time, last_checkpoint;
	double wus_per_second;

	time(&start_time);
	time(&last_checkpoint);
	log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Starting at: %d\n", (int)start_time);

	DB_VALIDATOR_ITEM_SET validator;
	while (!one_pass) {
		char buf[256], mod_clause[256];
		std::vector<VALIDATOR_ITEM> items;
        
		check_stop_daemons();

		num_validated = 0;
		/********
			*       loop over entries that need to be checked
		 ********/
	        while (!validator.enumerate(bsm_app.id, one_pass_N_WU?one_pass_N_WU:SELECT_LIMIT, wu_id_modulus, wu_id_remainder, items)) {
			/********
				*	for testing purposes, pretend we did nothing
			 ********/
			if (update_db) did_something = true;

			insert_workunit(validator, items);

			num_validated++;
		}
		validated_wus += num_validated;
		time(&current_time);
		wus_per_second = (double)validated_wus/((double)current_time-(double)start_time);
		count_unsent_results(unsent_wus, bsm_app.id);
		log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "[appid: %d] validated %d results, wus/sec: %lf, unsent wus: %d\n", bsm_app.id, num_validated, wus_per_second, unsent_wus);

                if (wu_id_modulus)      sprintf(mod_clause, " and workunit.id %% %d = %d ", wu_id_modulus, wu_id_remainder);
                else                    strcpy(mod_clause, "");
                sprintf(buf, "where appid=%d and assimilate_state=%d %s limit %d", bsm_app.id, ASSIMILATE_READY, mod_clause, one_pass_N_WU ? one_pass_N_WU : 1000);
                         
		num_assimilated = 0;
		while (!wu.enumerate(buf)) {
//			log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "[%s] assimilating boinc WU %d; state=%d\n", wu.name, wu.id, wu.assimilate_state);

			if (update_db) {
				char buf[256];
				did_something = true;
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

		assimilated_wus += num_assimilated;
		time(&current_time);
		wus_per_second = (double)assimilated_wus/((double)current_time-(double)start_time);
		log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "[appid: %d] assimilated %d workunits, wus/sec: %lf, unsent wus: %d\n", bsm_app.id, num_assimilated, wus_per_second, unsent_wus);

		if (unsent_wus < unsent_wu_buffer) {
			log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Generating %d new workunits.\n", get_generation_rate());
			num_generated = generate_workunits();
			log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "Generated %d new workunits.\n", num_generated);
		}

/*		if (num_validated && (current_time - last_checkpoint) > checkpoint_time) {
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
*/		if (!one_pass) sleep(sleep_interval);
	}
}

void generate_search_workunits(char *search_name) {
        manage_search(search_name);
        update_workunit_info(0);
        generate_workunits();
}
