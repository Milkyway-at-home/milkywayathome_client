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

#include "star_points.h"
#include "milkyway_util.h"

#if DOUBLEPREC
  #define STARPOINTS_READ_STR "%lf %lf %lf\n"
#else
  #define STARPOINTS_READ_STR "%f %f %f\n"
#endif /* DOUBLEPREC */

static int freadStarPoints(FILE* data_file, StarPoints* sp)
{
    real x, y, z;
    unsigned int i;

    if (fscanf(data_file, "%u\n", &sp->number_stars) != 1)
    {
        perror("reading sp->number_stars");
        warn("Failed to read number of star points from file\n");
        return 1;
    }

    sp->stars = (mwvector*) mwMallocA(sizeof(mwvector) * sp->number_stars);
    for (i = 0; i < sp->number_stars; ++i)
    {
        if (fscanf(data_file, STARPOINTS_READ_STR, &x, &y, &z) != 3)
        {
            perror("star points");
            warn("Failed to read star points item\n");
            return 1;
        }

        SET_VECTOR(sp->stars[i], x, y, z);
    }

    return 0;
}

int readStarPoints(StarPoints* sp, const char* filename)
{
    int rc;
    FILE* f;

    f = mwOpenResolved(filename, "r");
    if (!f)
    {
        perror("Opening star points file");
        return 1;
    }

    rc = freadStarPoints(f, sp);
    fclose(f);

    return rc;
}

void freeStarPoints(StarPoints* sp)
{
    mwFreeA(sp->stars);
}

