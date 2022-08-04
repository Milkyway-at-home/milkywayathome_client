/*
 *  Copyright (c) 2011 Matthew Arsenault
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

#include "milkyway_show.h"
#include "milkyway_util.h"

char* showVector(mwvector v)
{
    char* buf;

    if (asprintf(&buf, "{ %g, %g, %g }", X(v), Y(v), Z(v)) < 0)
        mw_fail("asprintf() failed\n");

    return buf;
}

void printVector(mwvector v)
{
    char* buf = showVector(v);
    puts(buf);
    free(buf);
}

char* showMatrix(mwmatrix m)
{
    char* buf;

    if (0 > asprintf(&buf,
                     "{\n"
                     "  { %g, %g, %g }\n"
                     "  { %g, %g, %g }\n"
                     "  { %g, %g, %g }\n"
                     "}\n",
                     X(m[0]), Y(m[0]), Z(m[0]),
                     X(m[1]), Y(m[1]), Z(m[1]),
                     X(m[2]), Y(m[2]), Z(m[2])))
    {
        mw_fail("asprintf() failed\n");
    }

    return buf;
}

void printMatrix(mwmatrix m)
{
    char* buf = showMatrix(m);
    puts(buf);
    free(buf);
}

