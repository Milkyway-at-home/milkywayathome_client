#include "config.h"
#include <cstdio>
#include "time.h"

#include "boinc_db.h"
#include "backend_lib.h"
#include "util.h"
#include "sched_config.h"
#include "sched_msgs.h"
#include "sched_util.h"
#include "assimilate_handler.h"
#include "validate_util.h"
#include "string.h"

/********
	*	FGDO includes
 ********/
#include "../settings.h"
#include "search_manager.h"


using std::vector;
using std::string;

int parse_result_file(const char* output_file_name, char** search_name, double** parameters, char** metadata, double* fitness) {
	int number_parameters, i, retval;
	FILE* file = fopen(output_file_name, "r");

	if (file == NULL) return -1;

	(*search_name) = (char*)malloc(sizeof(char) * METADATA_SIZE);
	retval = fscanf(file, "%s\n", search_name);
	if (retval != 1) return -1;

	retval = fscanf(file, "%d:", &number_parameters);
	if (retval != 1) return -1;

	(*parameters) = (double*)malloc(sizeof(double) * number_parameters);
	for (i = 0; i < number_parameters; i++) {
		retval = fscanf(file, " %lf", &(*parameters)[i]);
		if (retval != 1) return -1;
	}
	fscanf(file, "\n");

	(*metadata) = (char*)malloc(sizeof(char) * METADATA_SIZE);
	for (i = 0; i < METADATA_SIZE && (*metadata)[i-1] != '\n'; i++) (*metadata)[i] = fgetc(file);

	retval = fscanf(file, "fitness: %lf", fitness);
	if (retval != 1) return -1;
	fclose(file);
	return 1;
}

/*
 *	A general handler for multiple types of searches
 */
int assimilate_handler(WORKUNIT& wu, vector<RESULT>& /*results*/, RESULT& canonical_result) {
	string output_file_name;
	char* search_name;
	char* metadata;
	double fitness;
	double* parameters;
	int success;
	SCOPE_MSG_LOG scope_messages(log_messages, SCHED_MSG_LOG::MSG_NORMAL);

	if (wu.canonical_resultid) {
		log_messages.printf_multiline(SCHED_MSG_LOG::MSG_DEBUG, canonical_result.xml_doc_out, "[%s] canonical result", wu.name);

		get_output_file_path(canonical_result, output_file_name);
		success = parse_result_file(output_file_name.c_str(), &search_name, &parameters, &metadata, &fitness);
		if (success <= 0) scope_messages.printf("Error parsing result file: %s [%s]\n", wu.name, search_name);

		scope_messages.printf("[%s] assimilating, search_name: [%s], fitness: [%lf], metadata: [%s]\n", wu.name, search_name, fitness, metadata);

		success = insert_workunit(search_name, fitness, parameters, metadata);
		if (success <= 0) scope_messages.printf("Unknown search: [%s]\n", search_name);

		free(search_name);
		free(metadata);
		free(parameters);
	} else {
		scope_messages.printf("[%s] No canonical result\n", wu.name);
	}

	if (wu.error_mask&WU_ERROR_COULDNT_SEND_RESULT)		log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: couldn't send a result\n", wu.name);
	if (wu.error_mask&WU_ERROR_TOO_MANY_ERROR_RESULTS)	log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many error results\n", wu.name);
	if (wu.error_mask&WU_ERROR_TOO_MANY_TOTAL_RESULTS)	log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many total results\n", wu.name);
	if (wu.error_mask&WU_ERROR_TOO_MANY_SUCCESS_RESULTS)	log_messages.printf(SCHED_MSG_LOG::MSG_CRITICAL, "[%s] Error: too many success results\n", wu.name);
	return 0;
}
const char *BOINC_RCSID_8f6a5a2d27 = "$Id: boinc_assimilate_handler.c,v 1.1 2008/10/30 21:07:02 deselt Exp $";
