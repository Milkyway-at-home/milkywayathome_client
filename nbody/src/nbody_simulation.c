/* ************************************************************************** */
/* nbody_simulation.c: hierarchical N-body code. */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include <stdlib.h>

#include "nbody.h"
#include "json_params.h"
#include "nbody_priv.h"

inline static void initState(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st)
{
    printf("Starting nbody system\n");

    st->tout       = st->tnow;       /* schedule first output */
    st->tree.rsize = ctx->tree_rsize;

    st->tnow = 0.0;                 /* reset elapsed model time */
    st->bodytab = (bodyptr) mallocSafe(ctx->model.nbody * sizeof(body));
    st->acctab  = (vector*) mallocSafe(ctx->model.nbody * sizeof(vector));

    generatePlummer(ctx, ic, st);    /* make test model */

  #if NBODY_OPENCL
    setupNBodyCL(ctx, st);
  #endif /* NBODY_OPENCL */

    /* CHECKME: Why does makeTree get used twice for the first step? */
    /* Take 1st step */

  #if !NBODY_OPENCL
    gravMap(ctx, st);
  #else
    gravMapCL(ctx, st);
  #endif /* !NBODY_OPENCL */

}

static void startRun(const NBodyCtx* ctx, InitialConditions* ic, NBodyState* st)
{
    InitialConditions fc;

    printf("Starting fresh nbody run\n");
    reverseOrbit(&fc, ctx, ic);
    initState(ctx, &fc, st);
}

/* stepSystem: advance N-body system one time-step. */
inline static void stepSystem(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    vector* a;
    vector dvel, dpos;
    const bodyptr endp = st->bodytab + ctx->model.nbody;
    const real dt = ctx->model.timestep;

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)    /* loop over all bodies */
    {
        MULVS(dvel, (vectorptr) a, 0.5 * dt);      /* get velocity increment */
        INCADDV(Vel(p), dvel);                     /* advance v by 1/2 step */
        MULVS(dpos, Vel(p), dt);                   /* get positon increment */
        INCADDV(Pos(p), dpos);                     /* advance r by 1 step */
    }

  #if !NBODY_OPENCL
    gravMap(ctx, st);
  #else
    gravMapCL(ctx, st);
  #endif /* !NBODY_OPENCL */

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)      /* loop over all bodies */
    {
        MULVS(dvel, (vectorptr) a, 0.5 * dt);   /* get velocity increment */
        INCADDV(Vel(p), dvel);                  /* advance v by 1/2 step */
    }

    st->tnow += dt;                           /* finally, advance time */
}

static void runSystem(const NBodyCtx* ctx, NBodyState* st)
{
    const real tstop = ctx->model.time_dwarf - ctx->model.timestep / 1024.0;

    while (st->tnow < tstop)
    {
        stepSystem(ctx, st);   /* advance N-body system */
        #if BOINC_APPLICATION
          nbodyCheckpoint(ctx, st);
        #endif

        #if 0 /* TODO: Some day this will allow printing at intervals, for making movies etc. */
          /* TODO: organize use of this output better since it only
           * half makes sense now with boinc */

          if (ctx->model.time_dwarf - st->tnow < 0.01 * ctx->model.timestep)
              output(ctx, st);

          st->tout += 1.0 / ctx->freqout;     /* schedule next data out */
        #endif
    }
}

static void endRun(NBodyCtx* ctx, NBodyState* st)
{
  /* Make final output */
  #if BOINC_APPLICATION && !BOINC_DEBUG
    boincOutput(ctx, st);
  #else
    output(ctx, st);
  #endif /* BOINC_APPLICATION && !BOINC_DEBUG */

  #if BOINC_APPLICATION
    closeCheckpoint(ctx);       /* We finished so kill the checkpoint */
  #endif

    nbodyCtxDestroy(ctx);     /* finish up output */
    nbodyStateDestroy(st);
}

/* Takes parsed json and run the simulation, using outFileName for
 * output. The mess with the different names is for the hacky way we
 * can switch precision easily */
#if DYNAMIC_PRECISION
  #ifdef DOUBLEPREC
    #define RUN_NBODY_SIMULATION runNBodySimulation_double
  #else
    #define RUN_NBODY_SIMULATION runNBodySimulation_float
  #endif /* DOUBLEPREC */
#else
  #define RUN_NBODY_SIMULATION runNBodySimulation
#endif /* DYNAMIC_PRECISION */


void RUN_NBODY_SIMULATION(json_object* obj,
                          const FitParams* fitParams,
                          const char* outFileName,
                          const char* checkpointFileName,
                          const char* histogramFileName,
                          const char* histoutFileName,
                          const int outputCartesian,
                          const int printTiming,
                          const int verifyOnly)
{
    NBodyCtx ctx         = EMPTY_CTX;
    InitialConditions ic = EMPTY_INITIAL_CONDITIONS;
    NBodyState st        = EMPTY_STATE;

    double ts = 0.0, te = 0.0;
    int rc;

    rc = getParamsFromJSON(&ctx, &ic, obj, fitParams);
    if (verifyOnly)
    {
        if (rc)
            printf("File failed\n");
        else
            printf("File is OK\n");
        nbody_finish(rc);
    }

    if (rc)
        fail("Failed to read input parameters file\n");

    ctx.outputCartesian = outputCartesian;
    ctx.outfilename     = outFileName;
    ctx.histogram       = histogramFileName;
    ctx.histout         = histoutFileName;
    ctx.cp.filename     = checkpointFileName;

    initOutput(&ctx);

    printContext(&ctx);

  #if BOINC_APPLICATION
    /* If the checkpoint exists, try to use it */
    if (boinc_file_exists(ctx.cp.filename))
    {
        printf("Checkpoint exists. Attempting to resume from it.\n");
        openCheckpoint(&ctx);

        /* When the resume fails, start a fresh run */
        if (thawState(&ctx, &st))
        {
            warn("Failed to resume checkpoint\n");
            closeCheckpoint(&ctx);     /* Something is wrong with this file */
            openCheckpoint(&ctx);      /* Make a new one */
            nbodyStateDestroy(&st);
            startRun(&ctx, &ic, &st);
        }
    }
    else   /* Otherwise, just start a fresh run */
    {
        openCheckpoint(&ctx);
        startRun(&ctx, &ic, &st);
    }
  #else
    startRun(&ctx, &ic, &st);
  #endif /* BOINC_APPLICATION */

    if (printTiming)     /* Time the body of the calculation */
        ts = get_time();

    runSystem(&ctx, &st);

    if (printTiming)
    {
        te = get_time();
        printf("Elapsed time for run = %g\n", te - ts);
    }

    // Get the likelihood
    if (printTiming)
        ts = get_time();

    real chisqans = chisq(&ctx, &st);

    if (!isnan(chisqans))
        printf("Run finished. chisq = %f\n", chisqans);
    else
        warn("Failed to calculate chisq\n");

    if (printTiming)
    {
        te = get_time();
        printf("Elapsed time for chisq = %g\n", te - ts);
    }

    endRun(&ctx, &st);

}

