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

#ifndef ASTRONOMY_STAR_POINTS_H
#define ASTRONOMY_STAR_POINTS_H

#include <stdio.h>

typedef struct star_points {
	int number_stars;
	double** stars;
} STAR_POINTS;

int	read_star_points(const char* file, STAR_POINTS* sp);
int	fread_star_points(FILE* data_file, STAR_POINTS* sp);
int	write_star_points(const char* file, STAR_POINTS* sp);
int	fwrite_star_points(FILE* data_file, STAR_POINTS* sp);
void	free_star_points(STAR_POINTS* sp);

void	split_star_points(STAR_POINTS* sp, int rank, int max_rank);

#ifdef GMLE_BOINC
	int	boinc_read_star_points(const char* file, STAR_POINTS* sp);
#endif

#endif
