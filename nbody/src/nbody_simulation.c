/* ************************************************************************** */
/* nbody_simulation.c: hierarchical N-body code. */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "json_params.h"
#include "nbody_priv.h"
#include "nbody.h"

inline static void gravmap(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    const bodyptr endp = st->bodytab + ctx->model.nbody;
    unsigned int nfcalc  = 0;         /* count force calculations */
    unsigned int n2bcalc = 0;         /* count body-body interactions */
    unsigned int nbccalc = 0;         /* count body-cell interactions */

    maketree(ctx, st);                /* build tree structure */

    for (p = st->bodytab; p < endp; p++)
    {
        /* loop over all bodies */
        hackgrav(ctx, st, p, Mass(p) > 0.0);     /* get force on each */
        ++nfcalc;                                /* count force calcs */
        n2bcalc += st->n2bterm;                  /* and 2-body terms */
        nbccalc += st->nbcterm;                  /* and body-cell terms */
    }

}

/* startrun: startup hierarchical N-body code. */
inline static void startrun(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st)
{
    srand48(ctx->seed);              /* set random generator */
    generatePlummer(ctx, ic, st);    /* make test model */

    st->tout       = st->tnow;            /* schedule first output */
    st->tree.rsize = ctx->tree_rsize;

    /* CHECKME: Why does maketree get used twice for the first step? */
    gravmap(ctx, st);               /* Take 1st step */
    output(ctx, st);                /* do initial output */
    st->nstep = 1;                  /* start counting steps */
}

/* stepsystem: advance N-body system one time-step. */
inline static void stepsystem(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    vector dvel, dpos;

    const bodyptr bodytab = st->bodytab;
    const bodyptr endp    = bodytab + ctx->model.nbody;
    const real dt         = inv(ctx->freq);         /* set basic time-step */

    for (p = bodytab; p < endp; p++)    /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);  /* get velocity increment */
        INCADDV(Vel(p), dvel);          /* advance v by 1/2 step */
        MULVS(dpos, Vel(p), dt);        /* get positon increment */
        INCADDV(Pos(p), dpos);          /* advance r by 1 step */
    }

    gravmap(ctx, st);

    for (p = bodytab; p < endp; p++)          /* loop over all bodies */
    {
        MULVS(dvel, Acc(p), 0.5 * dt);        /* get velocity increment */
        INCADDV(Vel(p), dvel);                /* advance v by 1/2 step */
    }

    st->nstep++;           /* count another time step */
    st->tnow += dt;        /* finally, advance time */
    output(ctx, st);       /* do major or minor output */
}


static void runSystem(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st)
{
    const real tstop = ctx->model.time_dwarf - 1.0 / (1024.0 * ctx->freq);

    startrun(ctx, ic, st);

    while (st->tnow < tstop)
        stepsystem(ctx, st);               /* advance N-body system */
}

/* Takes parsed json and run the simulation, using outFileName for
 * output */
#ifdef DYNAMIC_PRECISION
  #ifdef DOUBLEPREC
    void runNBodySimulation_double(json_object* obj, const char* outFileName)
  #else
    void runNBodySimulation_float(json_object* obj, const char* outFileName)
  #endif /* DOUBLEPREC */
#else
  void runNBodySimulation(json_object* obj, const char* outFileName)
#endif /* DYNAMIC_PRECISION */
{
    NBodyCtx ctx         = EMPTY_CTX;
    InitialConditions ic = EMPTY_INITIAL_CONDITIONS;
    NBodyState st        = EMPTY_STATE;

    ctx.outfilename = outFileName;

    get_params_from_json(&ctx, &ic, obj);
    initoutput(&ctx);

    printContext(&ctx);
    printInitialConditions(&ic);

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

}

