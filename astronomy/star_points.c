/*
 *  star_points.c
 *  Astronomy
 *
 *  Created by Travis Desell on 2/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

/****
	*	Astronomy includes
*****/
#include <stdio.h>
#include <stdlib.h>
#include "star_points.h"

int read_star_points(const char* filename, STAR_POINTS* sp) {
	FILE* data_file = fopen(filename, "r");
	return fread_star_points(data_file, filename, sp);
}

int fread_star_points(FILE* data_file, const char* filename, STAR_POINTS* sp) {
	int i;
	if (!data_file) {
		fprintf(stderr, "Couldn't find input file %s.\n", filename);
		return 1;
	}

	fscanf(data_file, "%d\n", &sp->number_stars);

	sp->stars = (double**)malloc(sizeof(double*) * sp->number_stars);
	for (i = 0; i < sp->number_stars; i++) {
		sp->stars[i] = (double*)malloc(sizeof(double)*3);
		fscanf(data_file, "%lf %lf %lf\n", &sp->stars[i][0], &sp->stars[i][1], &sp->stars[i][2]);
	}
	fclose(data_file);
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

	first_star = ((double)sp->number_stars) * (((double)rank)/((double)max_rank));
	last_star = ((double)sp->number_stars) * (((double)rank+1.0)/((double)max_rank));
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

	printf("rank [%d / %d] stars [%d]\n", rank, max_rank, sp->number_stars);
}
