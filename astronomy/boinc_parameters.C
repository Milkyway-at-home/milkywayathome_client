/*
 *  parameters.c
 *  Astronomy
 *
 *  Created by Travis Desell on 2/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

/****
         *      BOINC includes
****/
#ifdef _WIN32
	#include "boinc_win.h"
#else
	#include "config.h"
#endif

#ifndef _WIN32
	#include <cstdio>
	#include <cctype>
	#include <ctime>
	#include <cstring>
	#include <cstdlib>
	#include <csignal>
	#include <unistd.h>
#endif

#include "diagnostics.h"
#include "util.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"

#include <stdlib.h>
#include <math.h>

/****
	*	Astronomy includes
*****/
#include "parameters.h"

int boinc_read_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS *ap) {
	char input_path[512];
	int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));

	if (retval) {
		fprintf(stderr, "APP: error resolving parameters file [%s], %d\n", filename, retval);
		return retval;
	}

	FILE* data_file = boinc_fopen(input_path, "r");
	if (!data_file) {
		fprintf(stderr, "Couldn't find input file [%s] to read astronomy parameters.\n", filename);
		return 1;
	}

	fread_astronomy_parameters(data_file, ap);
	fclose(data_file);
	return 0;
}

int boinc_write_astronomy_parameters(const char* filename, ASTRONOMY_PARAMETERS *ap) {
	char input_path[512];
	int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));

	if (retval) {
		fprintf(stderr, "APP: error writing astronomy parameters [%s], %d\n", filename, retval);
		return retval;
	}

	FILE* data_file = boinc_fopen(input_path, "w");
	if (!data_file) {
		fprintf(stderr, "Couldn't find output file [%s] to write astronomy parameters.\n", filename);
		return 1;
	}

	fwrite_astronomy_parameters(data_file, ap);
	fclose(data_file);
	return 0;
}
