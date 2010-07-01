/* ************************************************************************** */
/* nbody_simulation.c: hierarchical N-body code. */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "json_params.h"
#include "nbody_priv.h"
#include "nbody.h"

#define DEFAULT_CHECKPOINT_FILE "nbody_checkpoint"

inline static void gravmap(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    const bodyptr endp = st->bodytab + ctx->model.nbody;

    maketree(ctx, st);                /* build tree structure */

    for (p = st->bodytab; p < endp; p++)        /* loop over all bodies */
        hackgrav(ctx, st, p, Mass(p) > 0.0);     /* get force on each */
}

/* startRun: startup hierarchical N-body code. */
inline static void startRun(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st)
{
    srand48(ctx->seed);              /* set random generator */
    st->tout       = st->tnow;       /* schedule first output */
    st->tree.rsize = ctx->tree_rsize;

    st->tnow = 0.0;                 /* reset elapsed model time */
    st->bodytab = (bodyptr) allocate(ctx->model.nbody * sizeof(body));

    generatePlummer(ctx, ic, st);    /* make test model */

    /* CHECKME: Why does maketree get used twice for the first step? */
    gravmap(ctx, st);               /* Take 1st step */
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

    st->tnow += dt;        /* finally, advance time */
}

static void runSystem(const NBodyCtx* ctx, NBodyState* st)
{
    const real tstop = ctx->model.time_dwarf - 1.0 / (1024.0 * ctx->freq);

    while (st->tnow < tstop)
    {
        stepsystem(ctx, st);   /* advance N-body system */
        #if BOINC_APPLICATION
          nbody_boinc_output(ctx, st);
        #else
          /* TODO: organize use of this output better since it only
           * half makes sense now with boinc */

          if (ctx->model.time_dwarf - st->tnow < 0.01 / ctx->freq)
          {
              output(ctx, st);
          }

          st->tout += 1.0 / ctx->freqout;     /* schedule next data out */
        #endif
    }
}

void endRun(NBodyCtx* ctx, NBodyState* st)
{
    /* Make final output */

  #if BOINC_APPLICATION && !BOINC_DEBUG
    boincOutput(ctx, st);
  #else
    output(ctx, st);
  #endif /* BOINC_APPLICATION && !BOINC_DEBUG */

    // Get the likelihood
    //chisqans = chisq();
    //printf("Run finished. chisq = %f\n", chisqans);

    nbody_ctx_destroy(ctx);               /* finish up output */
    nbody_state_destroy(st);

    /* We finished so kill the checkpoint */
    nbody_remove("nbody_checkpoint");

}

/* Takes parsed json and run the simulation, using outFileName for
 * output */
#ifdef DYNAMIC_PRECISION
  #ifdef DOUBLEPREC
    void runNBodySimulation_double(json_object* obj, const char* outFileName, const char* checkpointFileName)
  #else
    void runNBodySimulation_float(json_object* obj, const char* outFileName, const char* checkpointFileName)
  #endif /* DOUBLEPREC */
#else
    void runNBodySimulation(json_object* obj, const char* outFileName, const char* checkpointFileName)
#endif /* DYNAMIC_PRECISION */
{
    NBodyCtx ctx         = EMPTY_CTX;
    InitialConditions ic = EMPTY_INITIAL_CONDITIONS;
    NBodyState st        = EMPTY_STATE;

    ctx.outfilename = outFileName;
    ctx.cpFile = checkpointFileName ? checkpointFileName : DEFAULT_CHECKPOINT_FILE;

    get_params_from_json(&ctx, &ic, obj);
    initoutput(&ctx);

    printContext(&ctx);
    printInitialConditions(&ic);

    // Calculate the reverse orbit
    printf("Calculating reverse orbit...\n");
    integrate(&ctx, &ic);
    printf("done\n");

    printInitialConditions(&ic);

    printf("Running nbody system\n");

    startRun(&ctx, &ic, &st);
    runSystem(&ctx, &st);
    endRun(&ctx, &st);
}

#if BOINC_APPLICATION
/* Similar to runNBodySimulation, but resume from a checkpointed state
 * and don't integrate the orbit, etc. */
void resumeCheckpoint(json_object* obj, const char* outFileName, const char* checkpointFile)
{
    NBodyCtx ctx         = EMPTY_CTX;
    InitialConditions ic = EMPTY_INITIAL_CONDITIONS;  /* Don't actually use these now since not at start */
    NBodyState st        = EMPTY_STATE;

    ctx.outfilename = outFileName;
    ctx.cpFile      = checkpointFile;

    get_params_from_json(&ctx, &ic, obj);
    initoutput(&ctx);

    printf("Resuming nbody system\n");

    thawState(&ctx, &st);
    st.tree.rsize = ctx.tree_rsize;

    printf("System thawed. tnow = %g\n", st.tnow);
    runSystem(&ctx, &st);

    printf("Ran thawed system\n");
    endRun(&ctx, &st);
}
#endif /* BOINC_APPLICATION */


