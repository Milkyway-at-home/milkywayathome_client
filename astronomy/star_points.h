/*
 *  star_points.h
 *  Astronomy
 *
 *  Created by Travis Desell on 2/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef ASTRONOMY_STAR_POINTS_H
#define ASTRONOMY_STAR_POINTS_H

#include <stdio.h>

typedef struct star_points {
	int number_stars;
	double** stars;
} STAR_POINTS;

int	read_star_points(const char* file, STAR_POINTS* sp);
int	fread_star_points(FILE* data_file, const char* file, STAR_POINTS* sp);
void	free_star_points(STAR_POINTS* sp);

void	split_star_points(STAR_POINTS* sp, int rank, int max_rank);

#endif
