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

static void printHeader(FILE* f, int cartesian)
{
    fprintf(f, "# ignore %22s %22s %22s %22s %22s %22s %22s\n",
            cartesian ? "x" : "l",
            cartesian ? "y" : "b",
            cartesian ? "z" : "r",
            "v_x",
            "v_y",
            "v_z",
            "lambda"
        );
}

/* output: Print bodies */
static int outputBodies(FILE* f, const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf)
{
    Body* p;
    mwvector lbr;
    const Body* endp = st->bodytab + st->nbody;
    NBHistTrig ht;
    real lambda;

    nbGetHistTrig(&ht, &ctx->histogramParams);
    printHeader(f, nbf->outputCartesian);

    for (p = st->bodytab; p < endp; p++)
    {
        lambda = nbXYZToLambda(&ht, Pos(p), ctx->sunGCDist);

        fprintf(f, "%8d,", ignoreBody(p));  /* Print if model it belongs to is ignored */
        if (nbf->outputCartesian)
        {
            fprintf(f,
                    " %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f\n",
                    X(Pos(p)), Y(Pos(p)), Z(Pos(p)),
                    X(Vel(p)), Y(Vel(p)), Z(Vel(p)),
                    lambda);
        }
        else
        {
            lbr = cartesianToLbr(Pos(p), ctx->sunGCDist);
            fprintf(f,
                    " %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f\n",
                    L(lbr), B(lbr), R(lbr),
                    X(Vel(p)), Y(Vel(p)), Z(Vel(p)),
                    lambda);
        }
    }

    if (fflush(f))
    {
        perror("Body output flush");
        return TRUE;
    }

    return FALSE;
}

int nbWriteBodies(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
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

    if (f != DEFAULT_OUTPUT_FILE)
    {
        if (fclose(f) < 0)
        {
            perror("Error closing output file");
        }
    }

    return rc;
}

