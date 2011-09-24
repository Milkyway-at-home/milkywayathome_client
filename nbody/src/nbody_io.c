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

/* Low-level input and output operations. */

static void out_2vectors(FILE* str, mwvector vec1, mwvector vec2)
{
    fprintf(str,
            " %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E\n",
            X(vec1), Y(vec1), Z(vec1),
            X(vec2), Y(vec2), Z(vec2));
}

/* output: Print bodies */
static int outputBodies(FILE* f, const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf)
{
    Body* p;
    mwvector lbR;
    const Body* endp = st->bodytab + st->nbody;

    for (p = st->bodytab; p < endp; p++)
    {
        fprintf(f, "%d ", ignoreBody(p));  /* Print if model it belongs to is ignored */
        if (nbf->outputCartesian)     /* Probably useful for making movies and such */
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
    FILE* f;
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

