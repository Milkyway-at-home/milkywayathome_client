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

int write_star_points(const char* filename, STAR_POINTS* sp)
{
    int retval;
    FILE* f;

    f = mw_fopen(filename, "w");
    if (!f)
    {
        perror("Writing star file");
        return 1;
    }

    retval = fwrite_star_points(f, sp);
    fclose(f);
    return retval;
}

#if DOUBLEPREC
  #define STAR_POINTS_READ_STR "%lf %lf %lf\n"
#else
  #define STAR_POINTS_READ_STR "%f %f %f\n"
#endif /* DOUBLEPREC */

int fread_star_points(FILE* data_file, STAR_POINTS* sp)
{
    unsigned int i;
    if (!fscanf(data_file, "%u\n", &sp->number_stars))
    {
        fprintf(stderr, "Failed to read star points file\n");
        mw_finish(EXIT_FAILURE);
    }

    sp->stars = (real*) mallocSafe(sizeof(real) * VECTOR_SIZE * sp->number_stars);
    for (i = 0; i < sp->number_stars; ++i)
    {
        fscanf(data_file, STAR_POINTS_READ_STR, &XN(sp, i), &YN(sp, i), &ZN(sp, i));
    }

    return 0;
}

int fwrite_star_points(FILE* data_file, STAR_POINTS* sp)
{
    unsigned int i;
    fprintf(data_file, "%u\n", sp->number_stars);

    for (i = 0; i < sp->number_stars; i++)
        fprintf(data_file, STAR_POINTS_READ_STR, XN(sp, i), YN(sp, i), ZN(sp, i));

    return 0;
}

void free_star_points(STAR_POINTS* sp)
{
    free(sp->stars);
}

