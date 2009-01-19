#ifndef FGDO_SEARCH_PARAMETERS_H
#define FGDO_SEARCH_PARAMETERS_H

#include <stdio.h>

typedef struct search_parameters {
	char *search_name;

	int number_parameters;
	double* parameters;
	char* metadata;
	char app_version[128];
	char host_os[256];

	int has_result;
	double fitness;
} SEARCH_PARAMETERS;

void free_search_parameters(SEARCH_PARAMETERS *parameters);
void init_search_parameters(SEARCH_PARAMETERS **p, int number_parameters);
void set_search_parameters(SEARCH_PARAMETERS *p, char *search_name, int number_parameters, double* parameters, char* metadata);

int fread_search_parameters(FILE* file, SEARCH_PARAMETERS *parameters);
int fwrite_search_parameters(FILE* file, SEARCH_PARAMETERS *parameters);

int read_search_parameters(const char* filename, SEARCH_PARAMETERS *parameters);
int write_search_parameters(const char* filename, SEARCH_PARAMETERS *parameters);

#ifdef GMLE_BOINC
	int boinc_read_search_parameters(const char* filename, SEARCH_PARAMETERS* parameters);
	int boinc_read_search_parameters2(const char* filename, SEARCH_PARAMETERS* parameters, char* version);
	int boinc_write_search_parameters(const char* filename, SEARCH_PARAMETERS* parameters, double fitness);
#endif

#endif
