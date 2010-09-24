/*
  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
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

#include <stdlib.h>
#include "nbody.h"
#include "json_params.h"
#include "calc_params.h"
#include "nbody_priv.h"
#include "milkyway_util.h"

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */


static inline void initState(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st)
{
    warn("Starting nbody system\n");

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

    warn("Starting fresh nbody run\n");
    reverseOrbit(&fc, ctx, ic);
    initState(ctx, &fc, st);
}

#ifndef _OPENMP

static inline void advancePosAcc(NBodyState* st, const unsigned int nbody, const real dt)
{
    bodyptr p;
    vector* a;
    vector dvel;
    vector dpos;
    const bodyptr endp = st->bodytab + nbody;

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)    /* loop over all bodies */
    {
        MULVS(dvel, (vectorptr) a, 0.5 * dt);      /* get velocity increment */
        INCADDV(Vel(p), dvel);                     /* advance v by 1/2 step */
        MULVS(dpos, Vel(p), dt);                   /* get positon increment */
        INCADDV(Pos(p), dpos);                     /* advance r by 1 step */
    }
}


static inline void advanceVelocities(NBodyState* st, const unsigned int nbody, const real dt)
{
    bodyptr p;
    vector* a;
    vector dvel;

    const bodyptr endp = st->bodytab + nbody;

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)      /* loop over all bodies */
    {
        MULVS(dvel, (vectorptr) a, 0.5 * dt);   /* get velocity increment */
        INCADDV(Vel(p), dvel);                  /* advance v by 1/2 step */
    }
}

#else

static inline void advancePosAcc(NBodyState* st, const unsigned int nbody, const real dt)
{
    unsigned int i;
    vector dvel;
    vector dpos;
    bodyptr p;

  #pragma omp parallel for private(i, p, dvel, dpos) schedule(static)
    for (i = 0; i < nbody; ++i)     /* loop over all bodies */
    {
        p = &st->bodytab[i];
        MULVS(dvel, (vectorptr) &st->acctab[i], 0.5 * dt);      /* get velocity increment */
        INCADDV(Vel(p), dvel);                     /* advance v by 1/2 step */
        MULVS(dpos, Vel(p), dt);                   /* get positon increment */
        INCADDV(Pos(p), dpos);                     /* advance r by 1 step */
    }
}

static inline void advanceVelocities(NBodyState* st, const unsigned int nbody, const real dt)
{
    unsigned int i;
    vector dvel;

    #pragma omp parallel for private(i, dvel) schedule(static)
    for (i = 0; i < nbody; ++i)      /* loop over all bodies */
    {
        MULVS(dvel, (vectorptr) &st->acctab[i], 0.5 * dt);   /* get velocity increment */
        INCADDV(Vel(&st->bodytab[i]), dvel);                  /* advance v by 1/2 step */
    }
}

#endif /* _OPENMP */

/* stepSystem: advance N-body system one time-step. */
static inline void stepSystem(const NBodyCtx* ctx, NBodyState* st)
{
    const real dt = ctx->model.timestep;

    advancePosAcc(st, ctx->model.nbody, dt);

  #if !NBODY_OPENCL
    gravMap(ctx, st);
  #else
    gravMapCL(ctx, st);
  #endif /* !NBODY_OPENCL */

    advanceVelocities(st, ctx->model.nbody, dt);

    st->tnow += dt;                           /* finally, advance time */
}

#if BOINC_APPLICATION
static inline void nbodyCheckpoint(const NBodyCtx* ctx, const NBodyState* st)
{
    if (boinc_time_to_checkpoint())
    {
        if (writeCheckpoint(ctx, st))
            fail("Failed to write checkpoint\n");

        boinc_checkpoint_completed();
    }

    boinc_fraction_done(st->tnow / ctx->model.time_dwarf);
}
#endif /* BOINC_APPLICATION */

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

  #if BOINC_APPLICATION
    mw_report("Making final checkpoint\n");
    if (writeCheckpoint(ctx, st))
        fail("Failed to write final checkpoint\n");
  #endif

}

static void endRun(NBodyCtx* ctx, NBodyState* st, const real chisq)
{
  /* Make final output */
  #if BOINC_APPLICATION && !BOINC_DEBUG
    boincOutput(ctx, st, chisq);
  #else
    output(ctx, st);
    printf("chisq = %.15g\n", chisq);
  #endif /* BOINC_APPLICATION && !BOINC_DEBUG */

    nbodyCtxDestroy(ctx);     /* finish up output */
    nbodyStateDestroy(st);
}


#if BOINC_APPLICATION

/* Setup the run, taking care of checkpointing things when using BOINC */
static void setupRun(NBodyCtx* ctx, InitialConditions* ic, NBodyState* st)
{
    /* If the checkpoint exists, try to use it */
    if (boinc_file_exists(ctx->cp_resolved))
    {
        mw_report("Checkpoint exists. Attempting to resume from it.\n");
        /* When the resume fails, start a fresh run */
        if (readCheckpoint(ctx, st))
        {
            mw_report("Failed to read checkpoint\n");
            nbodyStateDestroy(st);
            startRun(ctx, ic, st);
        }
        else
        {
            mw_report("Successfully read checkpoint\n");
            /* We restored the useful state. Now still need to create
             * the workspace where new accelerations are
             * calculated. */
            st->acctab  = (vector*) mallocSafe(ctx->model.nbody * sizeof(vector));
          #if !NBODY_OPENCL
            gravMap(ctx, st);
          #else
            gravMapCL(ctx, st);
          #endif /* !NBODY_OPENCL */
        }
    }
    else   /* Otherwise, just start a fresh run */
    {
        startRun(ctx, ic, st);
    }
}

#else

/* When not using BOINC, we don't need to deal with the checkpointing */
static void setupRun(NBodyCtx* ctx, InitialConditions* ic, NBodyState* st)
{
    startRun(ctx, ic, st);
}

#endif /* BOINC_APPLICATION */

/* Set context fields read from command line flags */
static inline void nbodySetCtxFromFlags(NBodyCtx* ctx, const NBodyFlags* nbf)
{
    ctx->outputCartesian = nbf->outputCartesian;
    ctx->outputBodies    = nbf->printBodies;
    ctx->outputHistogram = nbf->printHistogram;
    ctx->outfilename     = nbf->outFileName;
    ctx->histogram       = nbf->histogramFileName;
    ctx->histout         = nbf->histoutFileName;
    ctx->cp_filename     = nbf->checkpointFileName;
}

/* Takes parsed json and run the simulation, using outFileName for
 * output. */
void runNBodySimulation(json_object* obj,                 /* The main configuration */
                        const FitParams* fitParams,       /* For server's arguments */
                        const NBodyFlags* nbf)            /* Misc. parameters to control output */
{
    NBodyCtx ctx         = EMPTY_CTX;
    InitialConditions ic = EMPTY_INITIAL_CONDITIONS;
    NBodyState st        = EMPTY_STATE;

    real chisq;
    double ts = 0.0, te = 0.0;
    int rc = 0;

    rc |= getParamsFromJSON(&ctx, &ic, obj);
    rc |= setCtxConsts(&ctx, fitParams, &ic, nbf->setSeed);

    if (nbf->verifyOnly)
    {
        if (rc)
            warn("File failed\n");
        else
            warn("File is OK\n");
        mw_finish(rc);
    }

    if (rc)
        fail("Failed to read input parameters file\n");

    nbodySetCtxFromFlags(&ctx, nbf);

  #if BOINC_APPLICATION
    if (resolveCheckpoint(&ctx))
        fail("Failed to resolve checkpoint\n");
  #endif /* BOINC_APPLICATION */

    if (initOutput(&ctx))
        fail("Failed to open output files\n");

    setupRun(&ctx, &ic, &st);

    if (nbf->printTiming)     /* Time the body of the calculation */
        ts = mwGetTime();

    runSystem(&ctx, &st);
    mw_report("Simulation complete\n");

    if (nbf->printTiming)
    {
        te = mwGetTime();
        printf("Elapsed time for run = %g\n", te - ts);
    }

    /* Get the likelihood */
    chisq = nbodyChisq(&ctx, &st);
    if (isnan(chisq))
        warn("Failed to calculate chisq\n");

    endRun(&ctx, &st, chisq);
}

