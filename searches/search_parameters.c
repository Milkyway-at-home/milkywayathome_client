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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "search_parameters.h"
#include "../util/settings.h"

/****
         *     BOINC includes
*****/
#ifdef BOINC_APPLICATION
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
#endif


void free_search_parameters(SEARCH_PARAMETERS *parameters) {
	free(parameters->search_name);
	free(parameters->parameters);
	free(parameters->metadata);
}

void init_search_parameters(SEARCH_PARAMETERS **p, int number_parameters) {
	(*p) = (SEARCH_PARAMETERS*)malloc(sizeof(SEARCH_PARAMETERS));
	(*p)->search_name = (char*)malloc(sizeof(char) * 1024);
	(*p)->metadata = (char*)malloc(sizeof(char) * METADATA_SIZE);
	(*p)->number_parameters = number_parameters;
	(*p)->parameters = (double*)malloc(sizeof(double) * number_parameters);

	memset((*p)->host_os, '\0', 256);
	memset((*p)->app_version, '\0', 128);
}

void set_search_parameters(SEARCH_PARAMETERS *p, char *search_name, int number_parameters, double* parameters, char* metadata) {
	strcpy(p->search_name, search_name);
	strcpy(p->metadata, metadata);

	if (p->number_parameters != number_parameters) {
		p->number_parameters = number_parameters;
		p->parameters = (double*)realloc(p->parameters, sizeof(double) * number_parameters);
	}
	memcpy(p->parameters, parameters, sizeof(double) * number_parameters);
}

int fread_metadata(FILE* file, char* metadata) {
	int c, i;
	fscanf(file, "metadata: ");
	c = fgetc(file);
	i = 0;
	while (i < METADATA_SIZE && c != '\n' && c != '\0') {
		if (c == 13) {
			c = fgetc(file);
			continue;
		}
		metadata[i] = c;
		c = fgetc(file);
		i++;
	}
	metadata[i] = '\0';
	if (c == '\0' || i == METADATA_SIZE) return 1;
	return 0;
}

int fread_metadata__alloc(FILE* file, char **metadata) {
	(*metadata) = (char*)malloc(sizeof(char) * METADATA_SIZE);
	return fread_metadata(file, (*metadata));
}


int fread_search_parameters(FILE* file, SEARCH_PARAMETERS *p) {
	int i, number_parameters;
	if (fscanf(file, "%s\n", p->search_name) != 1) return 1;
	if (fscanf(file, "parameters [%d]:", &number_parameters) != 1) return 1;

	if (p->number_parameters != number_parameters) {
		p->number_parameters = number_parameters;
		p->parameters = (double*)realloc(p->parameters, sizeof(double) * number_parameters);
	}
	for (i = 0; i < number_parameters; i++) {
		if (fscanf(file, " %lf", &(p->parameters[i])) != 1) return 1;
	}
	if (fscanf(file, "\n") != 0) return 1;
	return fread_metadata(file, p->metadata);
}

int fwrite_search_parameters(FILE* file, SEARCH_PARAMETERS *parameters) {
	int i;
	if (fprintf(file, "%s\n", parameters->search_name) < 0) return 1;
	if (fprintf(file, "parameters [%d]:", parameters->number_parameters) < 0) return 1;

	for (i = 0; i < parameters->number_parameters; i++) {
		if (fprintf(file, " %.15lf", parameters->parameters[i]) < 0) return 1;
	}
	if (fprintf(file, "\n") < 0) return 1;
	if (fprintf(file, "metadata: %s\n", parameters->metadata) < 0) return 1;

	return 0;
}

int read_search_parameters(const char* filename, SEARCH_PARAMETERS *parameters) {
	int retval;

#ifdef BOINC_APPLICATION
	char input_path[512];
	retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
	if (retval) {
		fprintf(stderr, "APP: error resolving search parameters file (for read): %d\n", retval);
		fprintf(stderr, "\tfilename: %s\n", filename);
		fprintf(stderr, "\tresolved input path: %s\n", input_path);
		return retval;
	}

	FILE* data_file = boinc_fopen(input_path, "r");
#else
	FILE* data_file = fopen(filename, "r");
#endif

	if (!data_file) {
		fprintf(stderr, "Error reading search parameters: Couldn't find input file %s.\n", filename);
		return 1;
	}
        
	retval = fread_search_parameters(data_file, parameters);
	fclose(data_file);
	return retval;
}

int write_search_parameters(const char* filename, SEARCH_PARAMETERS *parameters) {
	int retval;
	FILE* data_file = fopen(filename, "w+");

	if (!data_file) {
		fprintf(stderr, "Error writing search parameters: Couldn't find input file %s.\n", filename);
		return 1;
	}
        
	retval = fwrite_search_parameters(data_file, parameters);
	fclose(data_file);
	return retval;
}
