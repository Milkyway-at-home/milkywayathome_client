/* ************************************************************************** */
/* CODE.C: hierarchical N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "nbody.h"
#include "types.h"
#include "stdinc.h"
#include "json_params.h"

Tree t = EMPTY_TREE;

/* main: toplevel routine for hierarchical N-body code. */
int main(int argc, char* argv[])
{
    NBodyCtx ctx         = EMPTY_CTX;
    InitialConditions ic = EMPTY_INITIAL_CONDITIONS;
    NBodyState st        = EMPTY_STATE;

    float chisqans = 0.0;

    initNBody(&ctx, &ic, argc, (const char**) argv);
    initoutput(&ctx);

    printContext(&ctx);
    printInitialConditions(&ic);

    /* FIXME */
    t.rsize = 4.0;

    // Calculate the reverse orbit
    printf("Calculating reverse orbit...");
    integrate(&ctx, &ic);
    printf("done\n");

    printInitialConditions(&ic);

    printf("Running nbody system\n");
    runSystem(&ctx, &ic, &st);

    printf("Running system done\n");
    // Get the likelihood
    //chisqans = chisq();
    //printf("Run finished. chisq = %f\n", chisqans);

    nbody_ctx_destroy(&ctx);               /* finish up output */

    return 0;
}

