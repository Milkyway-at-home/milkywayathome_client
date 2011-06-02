/*
  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
  Copyright (c) 2010 Matthew Arsenault, Travis Desell, Boleslaw
    Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
    Rensselaer Polytechnic Institute.
  Copyright (c) 2002-2006 John M. Fregeau, Richard Campbell, Jeff Molofee

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

#include "nbody.h"
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "nbody_grav.h"
#include "nbody_show.h"
#include "nbody_lua.h"
#include "nbody_shmem.h"


static inline int nbodyTimeToCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
  #if BOINC_APPLICATION
    return boinc_time_to_checkpoint();
  #else
    time_t now;

    if (ctx->checkpointT < 0 || ((now = time(NULL)) - st->lastCheckpoint) < ctx->checkpointT)
        return FALSE;
    else
    {
        st->lastCheckpoint = now;
        return TRUE;
    }
  #endif /* BOINC_APPLICATION */
}

static inline void nbodyCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
    if (nbodyTimeToCheckpoint(ctx, st))
    {
        if (writeCheckpoint(ctx, st))
            fail("Failed to write checkpoint\n");

      #if BOINC_APPLICATION
        boinc_checkpoint_completed();
      #endif /* BOINC_APPLICATION */
    }

  #if BOINC_APPLICATION
    boinc_fraction_done(st->tnow / ctx->timeEvolve);
  #endif /* BOINC_APPLICATION */
}


static int runSystem(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    const real tstop = ctx->timeEvolve - ctx->timestep / 1024.0;

    if (nbf->visualizer)
    {
        launchVisualizer(st, nbf->visArgs);
    }

    while (st->tnow < tstop)
    {
        updateDisplayedBodies(st);

        if (stepSystem(ctx, st))   /* advance N-body system */
            return 1;

        nbodyCheckpoint(ctx, st);
    }

    if (BOINC_APPLICATION || ctx->checkpointT >= 0)
    {
        mw_report("Making final checkpoint\n");
        if (writeCheckpoint(ctx, st))
            return warn1("Failed to write final checkpoint\n");
    }

    return 0;
}

static NBodyStatus setupRun(NBodyCtx* ctx, NBodyState* st, HistogramParams* hp, const NBodyFlags* nbf)
{
    if (resolveCheckpoint(st, nbf->checkpointFileName))
    {
        warn("Failed to resolve checkpoint\n");
        return NBODY_ERROR;
    }

    /* If the checkpoint exists, try to use it */
    if (nbf->ignoreCheckpoint || !resolvedCheckpointExists(st))
    {
        if (setupNBody(ctx, st, hp, nbf))
        {
            warn("Failed to read input parameters file\n");
            return NBODY_ERROR;
        }
    }
    else
    {
        mw_report("Checkpoint exists. Attempting to resume from it.\n");

        if (nbf->inputFile && !BOINC_APPLICATION)
            warn("Warning: input file '%s' unused\n", nbf->inputFile);

        if (readCheckpoint(ctx, st))
        {
            mw_report("Failed to read checkpoint\n");
            destroyNBodyState(st);
            return NBODY_CHECKPOINT_ERROR;
        }
        else
        {
            mw_report("Successfully read checkpoint\n");
        }
    }

    return gravMap(ctx, st); /* Start 1st step */
}

/* Set context fields read from command line flags */
static inline void nbodySetCtxFromFlags(NBodyCtx* ctx, const NBodyFlags* nbf)
{
    ctx->checkpointT = nbf->checkpointPeriod;
}

int verifyFile(const NBodyFlags* nbf)
{
    int rc;
    NBodyCtx ctx  = EMPTY_NBODYCTX;
    NBodyState st = EMPTY_NBODYSTATE;

    rc = setupNBody(&ctx, &st, &ctx.histogramParams, nbf);
    if (rc)
        warn("File failed\n");
    else
    {
        warn("File is OK\n");
        printNBodyCtx(&ctx);
        printHistogramParams(&ctx.histogramParams);
    }

    destroyNBodyState(&st);

    return rc;
}

/* FIXME */
static NBodyCtx _ctx = EMPTY_NBODYCTX;
static NBodyState _st = EMPTY_NBODYSTATE;


int runNBodySimulation(const NBodyFlags* nbf)
{
    NBodyCtx* ctx = &_ctx;
    NBodyState* st = &_st;

    int rc = 0;
    real chisq;
    double ts = 0.0, te = 0.0;

    nbodySetCtxFromFlags(ctx, nbf);
    if (setupRun(ctx, st, &ctx->histogramParams, nbf))
    {
        return warn1("Failed to setup run\n");
    }

    if (initOutput(st, nbf))
    {
        return warn1("Failed to open output files\n");
    }

    if (createSharedScene(st, ctx, nbf->inputFile))
    {
        return warn1("Failed to create shared scene\n");
    }

    ts = mwGetTime();
    rc = runSystem(ctx, st, nbf);
    if (rc)
        return warn1("Error running system\n");
    te = mwGetTime();

    if (nbf->printTiming)
    {
        printf("<run_time> %f </run_time>\n", te - ts);
    }

    /* Get the likelihood */
    chisq = nbodyChisq(ctx, st, nbf, &ctx->histogramParams);
    if (isnan(chisq))
    {
        warn("Failed to calculate chisq\n");
        rc = 1;
    }

    finalOutput(ctx, st, nbf, chisq);
    destroyNBodyState(st);

    return rc;
}

