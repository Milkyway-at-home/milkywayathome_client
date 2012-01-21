/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole, Dave Przybylo
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2010 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2010 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "star_points.h"
#include "milkyway_util.h"

static int freadStarPoints(FILE* data_file, StarPoints* sp)
{
    double x, y, z;
    unsigned int i;

    if (fscanf(data_file, "%u\n", &sp->number_stars) != 1)
    {
        mwPerror("Failed to read number of star points from file\n");
        return 1;
    }

    sp->stars = (mwvector*) mwMallocA(sizeof(mwvector) * sp->number_stars);
    for (i = 0; i < sp->number_stars; ++i)
    {
        if (fscanf(data_file, "%lf %lf %lf\n", &x, &y, &z) != 3)
        {
            mwPerror("Failed to read star points item\n");
            return 1;
        }

        SET_VECTOR(sp->stars[i], (real) x, (real) y, (real) z);
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
        mwPerror("Opening star points file '%s'", filename);
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

