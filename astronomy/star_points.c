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

/****
	*	Astronomy includes
*****/
#include <stdio.h>
#include <stdlib.h>
#include "star_points.h"

/****
         *     BOINC includes
*****/
#ifdef GMLE_BOINC
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


int read_star_points(const char* filename, STAR_POINTS* sp) {
	int retval;
	FILE* data_file = fopen(filename, "r");

	if (!data_file) {
		fprintf(stderr, "Couldn't find input file %s.\n", filename);
		return 1;
	}

	retval = fread_star_points(data_file, sp);
	fclose(data_file);
	return retval;
}

int write_star_points(const char* filename, STAR_POINTS* sp) {
	int retval;
	FILE* data_file = fopen(filename, "w");
	if (!data_file) {
		fprintf(stderr, "Couldn't find input file %s.\n", filename);
		return 1;
	}

	retval = fwrite_star_points(data_file, sp);
	fclose(data_file);
	return retval;
}

int fread_star_points(FILE* data_file, STAR_POINTS* sp) {
	int i;
	fscanf(data_file, "%d\n", &sp->number_stars);

	sp->stars = (double**)malloc(sizeof(double*) * sp->number_stars);
	for (i = 0; i < sp->number_stars; i++) {
		sp->stars[i] = (double*)malloc(sizeof(double)*3);
		fscanf(data_file, "%lf %lf %lf\n", &sp->stars[i][0], &sp->stars[i][1], &sp->stars[i][2]);
	}
	return 0;
}

int fwrite_star_points(FILE* data_file, STAR_POINTS* sp) {
	int i;
	fprintf(data_file, "%d\n", sp->number_stars);

	for (i = 0; i < sp->number_stars; i++) {
		fprintf(data_file, "%lf %lf %lf\n", sp->stars[i][0], sp->stars[i][1], sp->stars[i][2]);
	}
	return 0;
}

void free_star_points(STAR_POINTS* sp) {
	int i;
	for (i = 0; i < sp->number_stars; i++) {
		free(sp->stars[i]);
	}
	free(sp->stars);
}

void split_star_points(STAR_POINTS* sp, int rank, int max_rank) {
	int first_star, last_star, num_stars;
	int i;
	double** new_stars;

	if (rank == 0 && max_rank == 0) return;

	first_star = (int) (((double)sp->number_stars) * (((double)rank)/((double)max_rank)));
	last_star = (int) (((double)sp->number_stars) * (((double)rank+1.0)/((double)max_rank)));
	num_stars = last_star-first_star;
	new_stars = (double**)malloc(sizeof(double*) * num_stars);

	for (i = 0; i < num_stars; i++) {
		new_stars[i] = (double*)malloc(sizeof(double) * 3);
		new_stars[i][0] = sp->stars[i+first_star][0];
		new_stars[i][1] = sp->stars[i+first_star][1];
		new_stars[i][2] = sp->stars[i+first_star][2];
	}
	free_star_points(sp);
	sp->stars = new_stars;
	sp->number_stars = num_stars;

	printf("[worker: %d] using [%d] stars\n", rank, sp->number_stars);
}

#ifdef GMLE_BOINC
	int boinc_read_star_points(const char* filename, STAR_POINTS* sp) {
		char input_path[512];
		int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
			if (retval) {
			fprintf(stderr, "APP: error resolving star points file %d\n", retval);
			return retval;
		}

		FILE* data_file = boinc_fopen(input_path, "r");
		retval = fread_star_points(data_file, sp);
		fclose(data_file);
		return retval;
	}
#endif
