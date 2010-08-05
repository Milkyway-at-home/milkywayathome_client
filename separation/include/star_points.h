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

#ifndef _STAR_POINTS_H_
#define _STAR_POINTS_H_

#include <stdio.h>

#define VECTOR_SIZE 3

/* Get the xth component of the nth item in STAR_POINTS */
#define VN(sp, n) ((sp->stars)[VECTOR_SIZE * (n)])
#define XN(sp, n) VN(sp,n)
#define YN(sp, n) ((sp->stars)[VECTOR_SIZE * (n) + 1])
#define ZN(sp, n) ((sp->stars)[VECTOR_SIZE * (n) + 2])

typedef struct star_points
{
    unsigned int number_stars;
    double* stars;
} STAR_POINTS;

#define EMPTY_STAR_POINTS { 0, NULL }

int read_star_points(const char* file, STAR_POINTS* sp);
int fread_star_points(FILE* data_file, STAR_POINTS* sp);
int write_star_points(const char* file, STAR_POINTS* sp);
int fwrite_star_points(FILE* data_file, STAR_POINTS* sp);
void free_star_points(STAR_POINTS* sp);

void split_star_points(STAR_POINTS* sp, int rank, int max_rank);

#endif /* _STAR_POINTS_H_ */

