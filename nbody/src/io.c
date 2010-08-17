/* ************************************************************************** */
/* IO.C: I/O routines for export version of hierarchical N-body code. */
/* Public routines: inputdata(), initoutput(), stopoutput(), output(). */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include <string.h>
#include "nbody_priv.h"
#include "nbody_util.h"
#include "io.h"
#include "nbody_boinc.h"

void initOutput(NBodyCtx* ctx)
{
    ctx->outfile = ctx->outfilename ? nbodyOpenResolved(ctx->outfilename, "w") : DEFAULT_OUTPUT_FILE;
    if (ctx->outfile == NULL)
        fail("initOutput: cannot open output file %s\n", ctx->outfilename);
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
    fprintf(ctx->outfile, "<search_application>" BOINC_NBODY_APP_VERSION "</search_application>\n");
}

/* output: compute diagnostics and output data. */
void output(const NBodyCtx* ctx, const NBodyState* st)
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

    fflush(ctx->outfile);             /* drain output buffer */
}

void nbodyCtxDestroy(NBodyCtx* ctx)
{
    free(ctx->headline);
    if (ctx->outfile && ctx->outfile != stdout)
        fclose(ctx->outfile);
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


