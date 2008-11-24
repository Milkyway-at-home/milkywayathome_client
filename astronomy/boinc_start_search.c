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
#include "../util/settings.h"

/********
	*	Astronomy includes
 ********/
#include "../evaluation/boinc_add_workunit.h"
#include "parameters.h"
#include "star_points.h"


void print_usage() {
	printf("Usage:\n");
	printf("\t-app <application\n");
	printf("\t\tspecifies the application: (should be milkyway) -- required\n");
	printf("\t-cwd <working_directory>\n");
	printf("\t\tspecifies the directory to put the search in -- default: ./\n");
	printf("\t-s <search_name>\n");
	printf("\t\tspecifies the name of the search -- required.\n");
	printf("\t-gen #\n");
	printf("\t\tspecifies how many workunits to generate -- required\n");
	printf("\t-mw_stars <star_file>\n");
	printf("\t\tspecifies the stars file -- required\n");
	printf("\t-mw_parameters <parameters_file>\n");
	printf("\t\tspecifies the parameter file -- required\n");
	printf("\t-h\n");
	printf("\t\tprints this help message.\n");
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

	printf("parsing arguments\n");
	if (argc == 1) {
		print_usage();
		exit(0);
	}
        for (i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-h")) {
			print_usage();
			exit(0);
		} else if (!strcmp(argv[i], "-s")) {
			if (search_name != NULL) {
				printf("ERROR: multiple searches specified, can only create workunits for one search at a time.\n");
				return 0;
			}
			search_name = (char*)malloc(sizeof(char) * SEARCH_NAME_SIZE);
			strcpy(search_name, argv[++i]);
                } else if (!strcmp(argv[i], "-mw_parameters")) {
			parameters_name = (char*)malloc(sizeof(char) * FILENAME_SIZE);
			strcpy(parameters_name, argv[++i]);
		} else if (!strcmp(argv[i], "-mw_stars")) {
			stars_name = (char*)malloc(sizeof(char) * FILENAME_SIZE);
			strcpy(stars_name, argv[++i]);
		} else if (!strcmp(argv[i], "-cwd")) {
			set_working_directory(argv[++i]);
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

	printf("sending arguments to boinc search manager\n");
	init_boinc_search_manager(argc, argv);

	printf("generating workunits\n");
	generate_workunits();
}
