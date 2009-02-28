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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/********
	*	FGDO includes
 ********/
#include "../evaluation/search_manager.h"
#include "../evaluation/boinc_search_manager.h"
#include "../searches/asynchronous_newton_method.h"
#include "../searches/asynchronous_particle_swarm.h"
#include "../searches/bounds.h"
#include "../searches/search_parameters.h"
#include "../util/settings.h"
#include "../util/io_util.h"

/********
	*	Astronomy includes
 ********/
#include "../evaluation/boinc_add_workunit.h"
#include "parameters.h"
#include "star_points.h"

/********
	*	BOINC includes
 ********/
#include "util.h"


void print_usage() {
	printf("Usage:\n");
	printf("\t-app <application\n");
	printf("\t\tspecifies the application: (should be milkyway) -- required\n");
	printf("\t-cwd <working_directory>\n");
	printf("\t\tspecifies the directory to put the search in -- default: ./\n");
	printf("\t-s <search_name>\n");
	printf("\t\tspecifies the name of the search -- required\n");
	printf("\t-gen #\n");
	printf("\t\tspecifies how many workunits to generate -- required\n");
	printf("\t-mw_stars <star_file>\n");
	printf("\t\tspecifies the stars file -- required for a new search\n");
	printf("\t-mw_parameters <parameters_file>\n");
	printf("\t\tspecifies the parameter file -- required for a new search\n");
	printf("\t-h\n");
	printf("\t\tprints this help message.\n");
	exit(1);
}


void get_filename(char *filepath, char *filename) {
	int i, length;
	for (i = 0; i < (int)strlen(filepath); i++) {
		if (filepath[i] == '/') break;
	}
	if (i == (int)strlen(filepath)) i = 0;
	else i++;
	length = strlen(filepath) - i;
	strncpy(filename, &(filepath[i]), length);
	filename[length] = '\0';
}

int main(int argc, char** argv) {
	int i, retval, app_specified;
	char *astronomy_name, *astronomy_path, *wu_astronomy_path;
	char *star_name, *star_path, *wu_star_path;
	char *search_name, *wu_info_file;
	double *point, *step, *min_bound, *max_bound;

	search_name = NULL;
	app_specified = 0;

	printf("registering search\n");
	register_search(get_asynchronous_newton_method());
	register_search(get_asynchronous_particle_swarm());

	astronomy_path = NULL;
	astronomy_name = (char*)malloc(sizeof(char) * FILENAME_SIZE);
	wu_astronomy_path = (char*)malloc(sizeof(char) * FILENAME_SIZE);
	star_path = NULL;
	star_name = (char*)malloc(sizeof(char) * FILENAME_SIZE);
	wu_star_path = (char*)malloc(sizeof(char) * FILENAME_SIZE);
	wu_info_file = (char*)malloc(sizeof(char) * FILENAME_SIZE);

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
			astronomy_path = (char*)malloc(sizeof(char) * FILENAME_SIZE);
			strcpy(astronomy_path, argv[++i]);
			get_filename(astronomy_path, astronomy_name);
		} else if (!strcmp(argv[i], "-mw_stars")) {
			star_path = (char*)malloc(sizeof(char) * FILENAME_SIZE);
			strcpy(star_path, argv[++i]);
			get_filename(star_path, star_name);
		} else if (!strcmp(argv[i], "-cwd")) {
			set_working_directory(argv[++i]);
		} else if (!strcmp(argv[i], "-app")) {
			app_specified = 1;
			i++;
		}
        }

	if (!app_specified) {
		fprintf(stderr, "ERROR: application not specified.\n");
		print_usage();
		exit(0);
	}

	if (search_name == NULL) {
		fprintf(stderr, "ERROR: search name not specified.\n");
		print_usage();
		exit(0);
	}
	
	printf("sending arguments to boinc search manager\n");
	init_boinc_search_manager(argc, argv);

	if (!search_exists(search_name)) {
		printf("search doesnt exist\n");

		if (astronomy_path == NULL) {
			fprintf(stderr, "[start search] parameters file not specified, required for starting a new search.\n");
			print_usage();
			exit(0);
		}
		if (star_path == NULL) {
			fprintf(stderr, "[start search] stars file not specified, required for starting a new search.\n");
			print_usage();
			exit(0);
		}

		ASTRONOMY_PARAMETERS *ap = (ASTRONOMY_PARAMETERS*)malloc(sizeof(ASTRONOMY_PARAMETERS));
		retval = read_astronomy_parameters(astronomy_path, ap);
		if (retval) {
			fprintf(stderr, "ERROR: reading astronomy parameters\n");
			exit(0);
		}
	        get_search_parameters(ap, &point);
	        get_step(ap, &step);
		get_min_parameters(ap, &min_bound);
		get_max_parameters(ap, &max_bound);

		int number_parameters = get_optimized_parameter_count(ap);
		BOUNDS *bounds;
		int *in_radians = (int*)malloc(sizeof(int) * number_parameters);
		for (i = 0; i < number_parameters; i++) {
			in_radians[i] = 0;
		}
		for (i = 0; i < ap->number_streams; i++) {
			in_radians[2 + (i * 6) + 3] = 1;
			in_radians[2 + (i * 6) + 4] = 1;
		}
		new_bounds(&bounds, number_parameters, min_bound, max_bound, in_radians);
		asynchronous_search__init(argc, argv, number_parameters, point, step, bounds);

		/********
			*	Move input files
		 ********/
		config.download_path(astronomy_name, wu_astronomy_path);
		if (fopen(wu_astronomy_path, "r") == NULL) {
			/********
				*	Astronomy parameters not moved to the download diretory yet, do this.
			********/
			retval = write_astronomy_parameters(wu_astronomy_path, ap);
			printf("[start search] copied [%s] to [%s]\n", astronomy_path, wu_astronomy_path);
			if (retval) {
				fprintf(stderr, "[start search] could not open parameter file for write: %s\n", astronomy_path);
				exit(0);
			}
		} else {
			printf("[start search] [%s] already existed at [%s]\n", astronomy_path, wu_astronomy_path);
		}

		STAR_POINTS *sp = (STAR_POINTS*)malloc(sizeof(STAR_POINTS));
		retval = read_star_points(star_path, sp);
		if (retval) {
			fprintf(stderr, "[start search] could not open star file for read: %s\n", star_path);
			exit(0);
		}
		config.download_path(star_name, wu_star_path);
		if (fopen(wu_star_path, "r") == NULL) {
			/********   
				*	Stars not moved to the download diretory yet, do this.
			 ********/
			retval = write_star_points(wu_star_path, sp);
			printf("[start search] copied [%s] to [%s]\n", star_path, wu_star_path);
			if (retval) {
				fprintf(stderr, "[start search] could not open star file for write: %s\n", star_path);
				exit(0);
			}
		} else {
			printf("[start search] [%s] already existed at [%s]\n", star_path, wu_star_path);
		}

		/********
			*	Create workunit info
		 ********/
		WORKUNIT_INFO *wu_info = (WORKUNIT_INFO*)malloc(sizeof(WORKUNIT_INFO));

		double calc_prob_count = 0;
		for (i = 0; i < ap->number_integrals; i++) calc_prob_count += ap->integral[i]->r_steps * ap->integral[i]->mu_steps * ap->integral[i]->nu_steps;

		double integral_flops = calc_prob_count * (4.0 + 2.0 * ap->number_streams + ap->convolve * (56 + 58 * ap->number_streams));
		double likelihood_flops = sp->number_stars * (ap->convolve * (100.0 + ap->number_streams * 58.0) + 251.0 + ap->number_streams * 12.0 + 54.0);

		double flops = integral_flops + likelihood_flops;
		double credit = flops / 1000000000000.0;
		double multiplier = 7.0;
		credit *= multiplier;

		printf("awarded credit: %lf\n", credit);

		wu_info->number_parameters = number_parameters;
		printf("number parameters: %d\n", wu_info->number_parameters);

		wu_info->credit_str = (char*)malloc(sizeof(char) * 1024);
		wu_info->result_xml_path = (char*)malloc(sizeof(char) * 1024);

		sprintf(wu_info->credit_str, "<credit>%lf</credit>", credit);
		sprintf(wu_info->result_xml_path, "/export/www/boinc/milkyway/templates/a_result.xml");

		wu_info->template_filename = (char*)malloc(sizeof(char) * 1024);
		sprintf(wu_info->template_filename, "/export/www/boinc/milkyway/templates/milkyway_wu.xml");
		if (read_file_malloc(wu_info->template_filename, wu_info->template_file)) {
			fprintf(stderr, "[start search] could not read workunit result template: %s\n", wu_info->template_filename);
			exit(0);
		}

		wu_info->number_required_files = 2;
		wu_info->required_files = (char**)malloc(sizeof(char*) * 2);
		wu_info->required_files[0] = (char*)malloc(sizeof(char) * 1024);
		wu_info->required_files[1] = (char*)malloc(sizeof(char) * 1024);

		strcpy(wu_info->required_files[0], astronomy_name);
		strcpy(wu_info->required_files[1], star_name);

		wu_info->rsc_fpops_est = flops;
		printf("calculated fpops: %lf\n", wu_info->rsc_fpops_est);

		wu_info->rsc_fpops_bound = wu_info->rsc_fpops_est * 100;
		wu_info->rsc_memory_bound = 5e8;
		wu_info->rsc_disk_bound = 15e6;                      //15 MB
		wu_info->delay_bound = 60 * 60 * 24 * 5;             //5 days
		wu_info->min_quorum = 1;
		wu_info->target_nresults = 1;
		wu_info->max_error_results = 0;
		wu_info->max_total_results = 4;
		wu_info->max_success_results = 1;

		sprintf(wu_info_file, "%s/%s/workunit_info", get_working_directory(), search_name);
		retval = write_workunit_info(wu_info_file, wu_info);
		if (retval) {
			fprintf(stderr, "[start search] could not write workunit info: %s\n", wu_info_file);
			exit(0);
		}

		free_parameters(ap);
		free(ap);
		free_star_points(sp);
		free(sp);
	}
	generate_search_workunits(search_name);
}
