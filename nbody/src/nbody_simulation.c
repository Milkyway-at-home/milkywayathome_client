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

#if NBODY_OPENCL
  #include "nbody_cl.h"
#endif


static inline int nbodyTimeToCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
  #if BOINC_APPLICATION
    return boinc_time_to_checkpoint();
  #else
    time_t now;

    if (ctx->checkpointT < 0)
        return FALSE;

    now = time(NULL);
    if ((now - st->lastCheckpoint) > ctx->checkpointT)
    {
        st->lastCheckpoint = now;
        return TRUE;
    }

    return FALSE;
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


static NBodyStatus runSystem(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    const real tstop = ctx->timeEvolve - ctx->timestep / 1024.0;
    NBodyStatus rc = NBODY_SUCCESS;

    if (nbf->visualizer)
    {
        launchVisualizer(st, nbf->visArgs);
    }

    rc |= gravMap(ctx, st); /* Calculate accelerations for 1st step this episode */
    if (nbodyStatusIsFatal(rc))
        return rc;

    while (st->tnow < tstop)
    {
        updateDisplayedBodies(st);
        rc |= stepSystem(ctx, st);
        if (nbodyStatusIsFatal(rc))   /* advance N-body system */
            return rc;

        nbodyCheckpoint(ctx, st);
    }

    if (BOINC_APPLICATION || ctx->checkpointT >= 0)
    {
        mw_report("Making final checkpoint\n");
        if (writeCheckpoint(ctx, st))
        {
            warn("Failed to write final checkpoint\n");
            return NBODY_CHECKPOINT_ERROR;
        }
    }

    return rc;
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

    return NBODY_SUCCESS;
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

    NBodyStatus rc = NBODY_SUCCESS;
    real chisq;
    double ts = 0.0, te = 0.0;

    if (setupRun(ctx, st, &ctx->histogramParams, nbf))
    {
        warn("Failed to setup run\n");
        return NBODY_ERROR;
    }

    nbodySetCtxFromFlags(ctx, nbf); /* Do this after setup to avoid the setup clobbering the flags */
    if (initOutput(st, nbf))
    {
        warn("Failed to open output files\n");
        return NBODY_ERROR;
    }

    warn("allow incest = %d\n", ctx->allowIncest);

    if (createSharedScene(st, ctx, nbf->inputFile))
    {
        warn("Failed to create shared scene\n");
        return NBODY_ERROR;
    }

    ts = mwGetTime();

  #if NBODY_OPENCL
    if (nbf->noCL)
    {
        rc = runSystem(ctx, st, nbf);
    }
    else
    {
        rc = runSystemCL(ctx, st, nbf);
    }

  #else
    rc = runSystem(ctx, st, nbf);
  #endif /* NBODY_OPENCL */

    if (nbodyStatusIsFatal(rc))
    {
        warn("Error running system: %s (%d)\n", showNBodyStatus(rc), rc);
        return rc;
    }
    else if (nbodyStatusIsWarning(rc))
    {
        warn("System complete with warnings: %s (%d)\n", showNBodyStatus(rc), rc);
    }
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
        rc = NBODY_ERROR;
    }

    finalOutput(ctx, st, nbf, chisq);
    destroyNBodyState(st);

    return rc;
}

