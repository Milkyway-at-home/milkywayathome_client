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

#include "separation_priv.h"
#include "separation.h"
#include "star_points.h"

int read_star_points(const char* filename, STAR_POINTS* sp)
{
    int retval;
    FILE* data_file;
#if BOINC_APPLICATION
    char input_path[512];

    retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
    if (retval)
    {
        fprintf(stderr, "APP: error resolving star points file %d\n", retval);
        return retval;
    }

    data_file = boinc_fopen(input_path, "r");
#else
    data_file = fopen(filename, "r");
#endif

    if (!data_file)
    {
        fprintf(stderr, "Couldn't find input file %s.\n", filename);
        return 1;
    }

    retval = fread_star_points(data_file, sp);
    fclose(data_file);
    return retval;
}

int write_star_points(const char* filename, STAR_POINTS* sp)
{
    int retval;
    FILE* data_file = fopen(filename, "w");
    if (!data_file)
    {
        fprintf(stderr, "Couldn't find input file %s.\n", filename);
        return 1;
    }

    retval = fwrite_star_points(data_file, sp);
    fclose(data_file);
    return retval;
}

int fread_star_points(FILE* data_file, STAR_POINTS* sp)
{
    unsigned int i;
    if (!fscanf(data_file, "%u\n", &sp->number_stars))
    {
        fprintf(stderr, "Failed to read star points file\n");
        mw_finish(EXIT_FAILURE);
    }

    sp->stars = (double*) malloc(sizeof(double) * VECTOR_SIZE * sp->number_stars);
    for (i = 0; i < sp->number_stars; i++)
        fscanf(data_file, "%lf %lf %lf\n", &XN(sp, i), &YN(sp, i), &ZN(sp, i));

    return 0;
}

int fwrite_star_points(FILE* data_file, STAR_POINTS* sp)
{
    unsigned int i;
    fprintf(data_file, "%u\n", sp->number_stars);

    for (i = 0; i < sp->number_stars; i++)
        fprintf(data_file, "%lf %lf %lf\n", XN(sp, i), YN(sp, i), ZN(sp, i));

    return 0;
}

void free_star_points(STAR_POINTS* sp)
{
    free(sp->stars);
}

