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

void init_boinc_search_manager(int argc, char** argv) {
	int i, retval;
	char buf[256];

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
}

void start_search_manager() {
	DB_WORKUNIT wu;
	DB_RESULT canonical_result, result;
	bool did_something = false;
	int retval;

	while (!one_pass) {
		char buf[256], mod_clause[256];

		check_stop_daemons();

		if (wu_id_modulus)	sprintf(mod_clause, " and workunit.id %% %d = %d ", wu_id_modulus, wu_id_remainder);
		else			strcpy(mod_clause, "");
		sprintf(buf, "where appid=%d and assimilate_state=%d %s limit %d", bsm_app.id, ASSIMILATE_READY, mod_clause, one_pass_N_WU ? one_pass_N_WU : 1000);

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

			if (update_db) {
				sprintf(buf, "assimilate_state=%d, transition_time=%d", ASSIMILATE_DONE, (int)time(0));

				retval = wu.update_field(buf);
				if (retval) {
					log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] update failed: %d\n", wu.name, retval);
					exit(1);
				}
			}
		}
		if (did_something) boinc_db.commit_transaction();

		if (!one_pass) sleep(sleep_interval);
	}
}
