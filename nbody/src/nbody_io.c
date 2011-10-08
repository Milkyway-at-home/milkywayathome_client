/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
   Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#include "nbody_io.h"
#include "milkyway_util.h"
#include "nbody_coordinates.h"
#include "nbody_curses.h"

static void out_2vectors(FILE* f, mwvector v1, mwvector v2)
{
    fprintf(f,
            " %21.15f %21.15f %21.15f %21.15f %21.15f %21.15f\n",
            X(v1), Y(v1), Z(v1),
            X(v2), Y(v2), Z(v2));
}

static void printHeader(FILE* f, int cartesian)
{
    fprintf(f, "# ignore %21s %21s %21s %21s %21s %21s\n",
            cartesian ? "x" : "l",
            cartesian ? "y" : "b",
            cartesian ? "z" : "r",
            "v_x",
            "v_y",
            "v_z"
        );
}

/* output: Print bodies */
static int outputBodies(FILE* f, const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf)
{
    Body* p;
    mwvector lbR;
    const Body* endp = st->bodytab + st->nbody;

    printHeader(f, nbf->outputCartesian);

    for (p = st->bodytab; p < endp; p++)
    {
        fprintf(f, "%8d", ignoreBody(p));  /* Print if model it belongs to is ignored */
        if (nbf->outputCartesian)
        {
            out_2vectors(f, Pos(p), Vel(p));
        }
        else
        {
            lbR = cartesianToLbr(Pos(p), ctx->sunGCDist);
            out_2vectors(f, lbR, Vel(p));
        }
    }

    if (fflush(f))
    {
        perror("Body output flush");
        return TRUE;
    }

    return FALSE;
}

int finalOutput(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf, real chisq)
{
    FILE* f = DEFAULT_OUTPUT_FILE;
    int rc = 0;

    if (nbf->outFileName)
    {
        f = mwOpenResolved(nbf->outFileName, "w");
        if (!f)
        {
            mw_printf("Failed to open output file '%s'. Using default output\n", nbf->outFileName);
            f = DEFAULT_OUTPUT_FILE;
        }
    }

    /* Printing out the bodies will food the server. */
    if (nbf->printBodies)
    {
        mw_boinc_print(f, "<bodies>\n");
        rc = outputBodies(f, ctx, st, nbf);
        mw_boinc_print(f, "</bodies>\n");
    }

    fprintf(f, "<search_likelihood>%.15f</search_likelihood>\n", chisq);

    if (fclose(f) < 0)
    {
        perror("Error closing output file");
    }

    return rc;
}

