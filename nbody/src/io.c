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

#include <string.h>
#include "nbody_priv.h"
#include "io.h"
#include "milkyway_util.h"

int initOutput(NBodyCtx* ctx)
{
    ctx->outfile = ctx->outfilename ? mwOpenResolved(ctx->outfilename, "w") : DEFAULT_OUTPUT_FILE;
    if (ctx->outfile == NULL)
    {
        warn("initOutput: cannot open output file %s\n", ctx->outfilename);
        return TRUE;
    }

    return FALSE;
}

/* Low-level input and output operations. */

static void out_2vectors(FILE* str, vector vec1, vector vec2)
{
    fprintf(str, " %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E\n", vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2]);
}


/* Output with the silly xml stuff that BOINC uses */
void boincOutput(const NBodyCtx* ctx, const NBodyState* st, const real chisq)
{
    //not printing out the bodies because it will food the server.

    if (ctx->outputBodies)
    {
        fprintf(ctx->outfile, "<bodies>\n");
        output(ctx, st);
        fprintf(ctx->outfile, "</bodies>\n");
    }

    fprintf(ctx->outfile, "<search_likelihood>%.20g</search_likelihood>\n", chisq);
}

/* output: Print bodies */
int output(const NBodyCtx* ctx, const NBodyState* st)
{
    bodyptr p;
    vector lbR;
    const bodyptr endp = st->bodytab + ctx->model.nbody;

    for (p = st->bodytab; p < endp; p++)
    {
        if (ctx->outputCartesian)     /* Probably useful for making movies and such */
            out_2vectors(ctx->outfile, Pos(p), Vel(p));
        else
        {
            cartesianToLbr(ctx, lbR, Pos(p));
            out_2vectors(ctx->outfile, lbR, Vel(p));
        }
    }

    if (fflush(ctx->outfile))
    {
        perror("Body output flush");
        return TRUE;
    }

    return FALSE;
}

int nbodyCtxDestroy(NBodyCtx* ctx)
{
    free(ctx->headline);
    if (ctx->outfile && ctx->outfile != DEFAULT_OUTPUT_FILE)
    {
        if (fclose(ctx->outfile))
        {
            perror("closing output\n");
            return TRUE;
        }
    }

    return FALSE;
}

static void freeTree(Tree* t)
{
    nodeptr p, tmp;

    p = (nodeptr) t->root;
    while (p != NULL)
    {
        if (Type(p) == CELL)
        {
            tmp = More(p);
            free(p);
            p = tmp;
        }
        else                        /* skip over bodies */
            p = Next(p);
    }

    t->root = NULL;
    t->cellused = 0;
}

void nbodyStateDestroy(NBodyState* st)
{
    freeTree(&st->tree);
    free(st->bodytab);
    free(st->acctab);

  #if NBODY_OPENCL
    cleanupNBodyCL(st);
  #endif /* NBODY_OPENCL */
}


