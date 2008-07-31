/*
 *  star_points.c
 *  Astronomy
 *
 *  Created by Travis Desell on 2/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

/****
         *     BOINC includes
*****/
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

using std::string;

/****
	*	Astronomy includes
*****/
#include <stdio.h>
#include <stdlib.h>
#include "star_points.h"

int boinc_read_star_points(const char* filename, STAR_POINTS* sp) {
	char input_path[512];
	int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
	if (retval) {
		fprintf(stderr, "APP: error resolving star points file %d\n", retval);
		return retval;
	}

	FILE* data_file = boinc_fopen(input_path, "r");
	return fread_star_points(data_file, filename, sp);
}
