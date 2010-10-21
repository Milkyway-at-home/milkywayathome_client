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

#include <stdio.h>
#include <stdlib.h>
#include "star_points.h"
#include "milkyway_util.h"

#if DOUBLEPREC
  #define STAR_POINTS_READ_STR "%lf %lf %lf\n"
#else
  #define STAR_POINTS_READ_STR "%f %f %f\n"
#endif /* DOUBLEPREC */

static int fread_star_points(FILE* data_file, STAR_POINTS* sp)
{
    real x, y, z;
    unsigned int i;

    if (fscanf(data_file, "%u\n", &sp->number_stars) != 1)
    {
        perror("reading sp->number_stars");
        warn("Failed to read number of star points from file\n");
        return 1;
    }

    sp->stars = (mwvector*) mwMallocAligned(sizeof(mwvector) * sp->number_stars, sizeof(mwvector));
    for (i = 0; i < sp->number_stars; ++i)
    {
        if (fscanf(data_file, STAR_POINTS_READ_STR, &x, &y, &z) != 3)
        {
            perror("star points");
            warn("Failed to read star points item\n");
            return 1;
        }

        SET_VECTOR(sp->stars[i], x, y, z);
    }

    return 0;
}

int read_star_points(STAR_POINTS* sp, const char* filename)
{
    int rc;
    FILE* f;

    f = mwOpenResolved(filename, "r");
    if (!f)
    {
        perror("Opening star points file");
        return 1;
    }

    rc = fread_star_points(f, sp);
    fclose(f);

    return rc;
}

void free_star_points(STAR_POINTS* sp)
{
    mwAlignedFree(sp->stars);
}

