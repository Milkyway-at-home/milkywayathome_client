#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/********
	*	FGDO includes
 ********/
#include "../evaluation/search_manager.h"
#include "../evaluation/boinc_search_manager.h"
#include "../searches/asynchronous_newton_method.h"
#include "../searches/search_parameters.h"

/********
	*	Astronomy includes
 ********/
#include "add_workunit.h"
#include "parameters.h"
#include "star_points.h"


void print_arguments() {
	printf("Usage:\n");
	printf("\t-d <working_directory>, default: ./\n");
	printf("\t-s <search_name>, required.\n");
	printf("\t-wus <wus_to_generate>, required.\n");
	exit(1);
}

int main(int argc, char** argv) {
	int i, retval;
	char *search_name, *parameters_name, *stars_name, *search_qualifier, *outfile;
	double *point, *step;
	search_name = NULL;
	parameters_name = NULL;
	stars_name = NULL;

	printf("registering search\n");
	register_search(asynchronous_newton_method);

	printf("sending arguments to boinc search manager\n");
	init_boinc_search_manager(argc, argv, add_workunit);
	printf("parsing arguments\n");
        for (i = 0; i < argc; i++) {
                if (!strcmp(argv[i], "-s")) {
			search_name = (char*)malloc(sizeof(char) * SEARCH_NAME_SIZE);
			strcpy(search_name, argv[++i]);
                } else if (!strcmp(argv[i], "-mw_parameters")) {
			parameters_name = (char*)malloc(sizeof(char) * FILENAME_SIZE);
			strcpy(parameters_name, argv[++i]);
		} else if (!strcmp(argv[i], "-mw_stars")) {
			stars_name = (char*)malloc(sizeof(char) * FILENAME_SIZE);
			strcpy(stars_name, argv[++i]);
		}
        }

	ASTRONOMY_PARAMETERS *ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
	retval = read_astronomy_parameters(parameters_name, ap);
	if (retval) {
		fprintf(stderr, "ERROR: reading astronomy parameters\n");
		exit(0);
	}
        get_search_parameters(ap, &point);
        get_step(ap, &step);

	get_qualifier_from_name(search_name, &search_qualifier);
	printf("qualifier: %s\n", search_qualifier);
	if (!strcmp(search_qualifier, "nm")) {
		printf("creating newton method...\n");
		create_newton_method(search_name, 10, 300, ap->number_parameters, point, step);
		printf("created.\n");
	} else if (!strcmp(search_qualifier, "gs")) {
	} else if (!strcmp(search_qualifier, "de")) {
	} else if (!strcmp(search_qualifier, "pso")) {
	}
	free(search_qualifier);

	if (search_name == NULL) {
		fprintf(stderr, "ERROR: search name not specified for copying parameters\n");
		exit(0);
	}
	
	printf("copying parameters\n");
	outfile = (char*)malloc(sizeof(char) * FILENAME_SIZE);
	sprintf(outfile, "%s/%s/astronomy_parameters.txt", get_working_directory(), search_name);
	retval = write_astronomy_parameters(outfile, ap);
	if (retval) {
		fprintf(stderr, "ERROR: writing astronomy parameters to search directory\n");
		exit(0);
	}
	free(outfile);
	outfile = (char*)malloc(sizeof(char) * FILENAME_SIZE);


	printf("copying stars\n");
	STAR_POINTS *sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
	retval = read_star_points(stars_name, sp);
	if (retval) {
		fprintf(stderr, "ERROR: reading star points\n");
		exit(0);
	}
	sprintf(outfile, "%s/%s/stars.txt", get_working_directory(), search_name);
	retval = write_star_points(outfile, sp);
	if (retval) {
		fprintf(stderr, "ERROR: writing star points to search directory\n");
		exit(0);
	}
	free(outfile);

	printf("adding workunits\n");

	init_add_workunit();
	generate_workunits();
}
