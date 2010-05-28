// Berkeley Open Infrastructure for Network Computing
// http://boinc.berkeley.edu
// Copyright (C) 2005 University of California
//
// This is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation;
// either version 2.1 of the License, or (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Lesser General Public License for more details.
//
// To view the GNU Lesser General Public License visit
// http://www.gnu.org/copyleft/lesser.html
// or write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// validator - check and validate results, and grant credit
//  -app appname
//  [-d debug_level]
//  [-one_pass_N_WU N]      // Validate only N WU in one pass, then exit
//  [-one_pass]             // make one pass through WU table, then exit
//  [-mod n i]              // process only WUs with (id mod n) == i
//  [-max_credit_per_cpu_second X] // limit maximum credit per compute time
//  [-max_granted_credit X] // limit maximum granted credit to X
//  [-max_claimed_credit Y] // invalid if claims more than Y
//  [-grant_claimed_credit] // just grant whatever is claimed
//  [-update_credited_job]  // add userid/wuid pair to credited_job table
//  [-credit_from_wu]       // get credit from WU XML
//
// This program must be linked with two project-specific functions:
// check_set() and check_pair().
// See doc/validate.php for a description.

using namespace std;

#include "config.h"
#include <unistd.h>
#include <cmath>
#include <vector>

#include "boinc_db.h"
#include "util.h"
#include "str_util.h"
#include "error_numbers.h"

#include "sched_config.h"
#include "sched_util.h"
#include "sched_msgs.h"
#include "validate_util.h"

#define LOCKFILE "validate.out"
#define PIDFILE  "validate.pid"

#define SELECT_LIMIT    1000
#define SLEEP_PERIOD    5

int sleep_interval = SLEEP_PERIOD;

typedef enum { NEVER, DELAYED, IMMEDIATE, NO_CHANGE } TRANSITION_TIME;

//SCHED_CONFIG config;
char app_name[256];
int wu_id_modulus=0;
int wu_id_remainder=0;
int one_pass_N_WU=0;
bool one_pass = false;
double max_credit_per_cpu_second = 0;
double max_granted_credit = 0;
double max_claimed_credit = 0;
bool grant_claimed_credit = false;
bool update_credited_job = false;
bool credit_from_wu = false;

void update_error_rate(DB_HOST& host, bool valid) {
	if (host.error_rate > 1) host.error_rate = 1;
	if (host.error_rate <= 0) host.error_rate = 0.1;

	host.error_rate *= 0.95;
	if (!valid) {
		host.error_rate += 0.05;
	}
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
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d][HOST#%d] claimed too much credit (%f) in too little CPU time (%f)\n", result.id, result.hostid, result.granted_credit, result.cpu_time);
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

double update_workunit(DB_VALIDATOR_ITEM_SET& validator, int valid_state, RESULT& result) {
	bool assimilate = false;
	double cpu_time = result.cpu_time;

	switch (valid_state) {
		case AS_VERIFY_VALID:
			// grant credit for valid results
			result.granted_credit = grant_claimed_credit ? result.claimed_credit : credit;
			if (max_granted_credit && result.granted_credit > max_granted_credit) result.granted_credit = max_granted_credit

			if (result.cpu_time > 0.3 * round_trip_time) cpu_time = 0.2 * round_trip_time;
			if (max_credit_per_cpu_second && (result.granted_credit / cpu_time) > max_credit_per_cpu_second) result.granted_credit = cpu_time * max_credit_per_cpu_second;

			result.validate_state = VALIDATE_STATE_VALID;
			if (update_db) {
				retval = is_valid(result, wu);
				if (retval) log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "[RESULT#%d %s] is_valid() failed: %d\n", result.id, result.name, retval);
			}
			log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "[RESULT#%d %s] Granted %f credit to valid result [HOST#%d]\n", result.id, result.name, result.granted_credit, result.hostid);
			assimilate = true;
			break;
		case AS_VERIFY_INVALID:
			result.validate_state = VALIDATE_STATE_INVALID;
			if (update_db) is_invalid(result);
			assimilate = true;
			break;
		case AS_VERIFY_INCONCLUSIVE:
			ms->search->insert_parameters(ms->search_name, ms->search_data, insert_sp);
			result.validate_state = VALIDATE_STATE_INCONCLUSIVE;
			break;
	}

	if (update_db && assimilate) {
		log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL, "[RESULT#%d %s] granted_credit %f\n",  result.id, result.name, result.granted_credit);
		retval = validator.update_result(result);
		if (retval) log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[RESULT#%d %s] Can't update result: %d\n", result.id, result.name, retval);

		sprintf(buf, "assimilate_state=%d, transition_time=%d", ASSIMILATE_DONE, (int)time(0));
		retval = wu.update_field(buf);
		if (retval) {
			log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] update failed: %d\n", wu.name, retval);
			exit(1);
		}
	}
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


int insert_workunit(DB_VALIDATOR_ITEM_SET& validator, std::vector<VALIDATOR_ITEM>& items) {
	int canonical_result_index = -1;
	bool update_result, retry;
	TRANSITION_TIME transition_time = NO_CHANGE;
	int retval = 0, canonicalid = 0, x;
	double credit = 0;
	unsigned int i;

	WORKUNIT& wu = items[0].wu;

	if (items.size() > 1) printf("MORE THAN 1 ITEM IN WORKUNIT!\n");

	for (i = 0; i < items.size(); i++) {
		WORKUNIT& wu = items[i].wu;
		RESULT& result = items[i].res;
		MANAGED_SEARCH *ms;
		string output_file_name;
		double credit;

		ms = get_search_from_wu_name(wu.name);
		/********
			*	Check if this is a valid search.
		 ********/
		if (ms == NULL) {
			credit = update_workunit(validator, AS_VERIFY_VALID, result);
			print_message(scope_messages, ms->search_name, "completed search", VALIDATE_MSG, credit, result);
			transition_time = IMMEDIATE;
			continue;
		}

		/********
			*	Read the result file
		 ********/
		get_output_file_path(result, output_file_name);
		retval = boinc_read_search_parameters2(output_file_name.c_str(), insert_sp);
		if (retval) {
			credit = update_workunit(validator, AS_VERIFY_INVALID, result);
			print_message(scope_messages, ms->search_name, "error reading result", VALIDATE_MSG, credit, result);
			transition_time = IMMEDIATE;
			continue;
		}

		valid_state = ms->search->verify(ms->search_name, ms->search_data, insert_sp), result);
		/********
			*
		 ********/
		credit = update_workunit(validator, valid_state, result);
		if (valid_state == AS_VERIFY_IN_PROGRESS) {
			retval = ms->search->insert_parameters(ms->search_name, ms->search_data, insert_sp);
			transition_time = DELAYED;
			wu.need_validate = 1;
		} else {
			transition_time = IMMEDIATE;
		}
		print_message(scope_messages, ms->search_name, AS_MSG, VALIDATE_MSG, credit, result);
	}

        if (wu.error_mask&WU_ERROR_COULDNT_SEND_RESULT)         log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: couldn't send a result\n", wu.name);
        if (wu.error_mask&WU_ERROR_TOO_MANY_ERROR_RESULTS)      log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many error results\n", wu.name);
        if (wu.error_mask&WU_ERROR_TOO_MANY_TOTAL_RESULTS)      log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many total results\n", wu.name);
        if (wu.error_mask&WU_ERROR_TOO_MANY_SUCCESS_RESULTS)    log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many success results\n", wu.name);

	return 0;
}

// make one pass through the workunits with need_validate set.
// return true if there were any
//
bool do_validate_scan(APP& app) {
	DB_VALIDATOR_ITEM_SET validator;
	std::vector<VALIDATOR_ITEM> items;
	bool found=false;
	int retval;

	/********
		*	loop over entries that need to be checked
	 ********/
	while (1) {
		retval = validator.enumerate(app.id, one_pass_N_WU?one_pass_N_WU:SELECT_LIMIT, wu_id_modulus, wu_id_remainder, items);
		if (retval) break;
		retval = handle_wu(validator, items);
		if (!retval) found = true;
	}
	return found;
}

int main_loop() {
	int retval;
	DB_APP app;
	bool did_something;
	char buf[256];

	retval = boinc_db.open(config.db_name, config.db_host, config.db_user, config.db_passwd);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "boinc_db.open failed: %d\n", retval);
		exit(1);
	}

	sprintf(buf, "where name='%s'", app_name);
	retval = app.lookup(buf);
	if (retval) {
		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "can't find app %s\n", app_name);
		exit(1);
	}

	while (1) {
		check_stop_daemons();
		did_something = do_validate_scan(app);
		if (!did_something) {
			if (one_pass) break;
			sleep(sleep_interval);
		}
	}
	return 0;
}

// For use by user routines check_set() and check_match() that link to
// this code.
int boinc_validator_debuglevel=0;

int main(int argc, char** argv) {
    int i, retval;

#if 0
    int mypid=getpid();
    char debugcmd[512];
    sprintf(debugcmd, "ddd %s %d &", argv[0], mypid);
    system(debugcmd);
    sleep(30);
#endif

    const char *usage =
      "\nUsage: %s -app <app-name> [OPTIONS]\n"
      "Start validator for application <app-name>\n\n"
      "Optional arguments:\n"
      "  -one_pass_N_WU N                Validate at most N WUs, then exit\n"
      "  -one_pass                       Make one pass through WU table, then exit\n"
      "  -mod n i                        Process only WUs with (id mod n) == i\n"
      "  -max_credit_per_cpu_second X    Limit credit per compute second\n"
      "  -max_claimed_credit X           If a result claims more credit than this, mark it as invalid\n"
      "  -max_granted_credit X           Grant no more than this amount of credit to a result\n"
      "  -grant_claimed_credit           Grant the claimed credit, regardless of what other results for this workunit claimed\n"
      "  -update_credited_job            Add record to credited_job table after granting credit\n"
      "  -credit_from_wu                 Credit is specified in WU XML\n"
      "  -sleep_interval n               Set sleep-interval to n\n"
      "  -d level                        Set debug-level\n\n";

    if ( (argc > 1) && ( !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help") ) ) {
      printf (usage, argv[0] );
      exit(1);
    }


    check_stop_daemons();

    for (i=1; i<argc; i++) {
        if (!strcmp(argv[i], "-one_pass_N_WU")) {
            one_pass_N_WU = atoi(argv[++i]);
            one_pass = true;
        } else if (!strcmp(argv[i], "-sleep_interval")) {
            sleep_interval = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-one_pass")) {
            one_pass = true;
        } else if (!strcmp(argv[i], "-app")) {
            strcpy(app_name, argv[++i]);
        } else if (!strcmp(argv[i], "-d")) {
            boinc_validator_debuglevel=atoi(argv[++i]);
            log_messages.set_debug_level(boinc_validator_debuglevel);
        } else if (!strcmp(argv[i], "-mod")) {
            wu_id_modulus = atoi(argv[++i]);
            wu_id_remainder = atoi(argv[++i]);
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
        } else {
            fprintf(stderr, "Invalid option '%s'\nTry `%s --help` for more information\n", argv[i], argv[0]);
            log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "unrecognized arg: %s\n", argv[i]);
            exit(1);
        }
    }

    // -app is required
    if ( app_name[0] == 0 ) {
      fprintf (stderr, "\nERROR: use '-app' to specify the application to run the validator for.\n");
      printf (usage, argv[0] );
      exit(1);
    }

    retval = config.parse_file("..");
    if (retval) {
        log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL,
            "Can't parse ../config.xml: %s\n", boincerror(retval)
        );
        exit(1);
    }

    log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL,
        "Starting validator, debug level %d\n", log_messages.debug_level
    );
    if (wu_id_modulus) {
        log_messages.printf(SCHED_MSG_LOG::MSG_NORMAL,
            "Modulus %d, remainder %d\n", wu_id_modulus, wu_id_remainder
        );
    }

    install_stop_signal_handler();

    main_loop();
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

	DB_VALIDATOR_ITEM_SET validator;
	while (!one_pass) {
		char buf[256], mod_clause[256];

		check_stop_daemons();

		if (wu_id_modulus)	sprintf(mod_clause, " and workunit.id %% %d = %d ", wu_id_modulus, wu_id_remainder);
		else			strcpy(mod_clause, "");
		sprintf(buf, "where appid=%d and assimilate_state=%d %s limit %d", bsm_app.id, ASSIMILATE_READY, mod_clause, one_pass_N_WU ? one_pass_N_WU : 1000);

		num_assimilated = 0;


		std::vector<VALIDATOR_ITEM> items;
		bool found=false;
		int retval;

		/********
			*       loop over entries that need to be checked
		 ********/
	        while (!validator.enumerate(app.id, one_pass_N_WU?one_pass_N_WU:SELECT_LIMIT, wu_id_modulus, wu_id_remainder, items) {
			/********
				*	for testing purposes, pretend we did nothing
			 ********/
			if (update_db) did_something = true;
			log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "[%s] validating/assimilating boinc WU %d; state=%d\n", wu.name, wu.id, wu.assimilate_state);

			sprintf(buf, "where workunitid=%d", wu.id);
			while (!result.enumerate(buf)) {
				results.push_back(result);
				if (result.id == wu.canonical_resultid) canonical_result = result;
			}

			insert_workunit(validator, items);

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
