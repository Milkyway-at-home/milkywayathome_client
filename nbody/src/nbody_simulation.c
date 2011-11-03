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
#include "nbody_curses.h"
#include "nbody_defaults.h"

#if NBODY_OPENCL
  #include "nbody_cl.h"
#endif


static inline int nbTimeToCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
    time_t now;

    if (BOINC_APPLICATION)
    {
        return mw_time_to_checkpoint();
    }

    if (ctx->checkpointT < 0)
        return FALSE;

    now = time(NULL);
    if ((now - st->lastCheckpoint) > ctx->checkpointT)
    {
        st->lastCheckpoint = now;
        return TRUE;
    }

    return FALSE;
}

static inline NBodyStatus nbCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
    if (nbTimeToCheckpoint(ctx, st))
    {
        if (nbWriteCheckpoint(ctx, st))
            return NBODY_CHECKPOINT_ERROR;

        mw_checkpoint_completed();
    }

    return NBODY_SUCCESS;
}

static inline void nbReportProgress(const NBodyCtx* ctx, NBodyState* st, int reportProgress)
{
    double frac = (double) st->step / (double) ctx->nStep;

    mw_fraction_done(frac);

    if (reportProgress)
    {
        mw_mvprintw(0, 0,
                    "Running: %f / %f (%f%%)\n",
                    frac * ctx->timeEvolve,
                    ctx->timeEvolve,
                    100.0 * frac
            );

        mw_refresh();
    }
}

/* If enough time has passed, record the next center of mass position */
static void nbAddTracePoint(const NBodyCtx* ctx, NBodyState* st)
{
    int i = st->step * N_ORBIT_TRACE_POINTS / ctx->nStep;

    if (st->usesExact) /* FIXME?. We don't get the CM without the tree */
        return;

    if (i >= N_ORBIT_TRACE_POINTS) /* Just in case */
        return;

    if (X(st->orbitTrace[i]) < REAL_MAX)
        return;

    st->orbitTrace[i] = Pos(st->tree.root);
}

static NBodyStatus nbRunSystem(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    NBodyStatus rc = NBODY_SUCCESS;

    if (nbf->visualizer)
    {
        nbLaunchVisualizer(st, nbf->visArgs);
    }

    rc |= nbGravMap(ctx, st); /* Calculate accelerations for 1st step this episode */
    if (nbStatusIsFatal(rc))
        return rc;

    while (st->step < ctx->nStep)
    {
        nbAddTracePoint(ctx, st);
        nbUpdateDisplayedBodies(ctx, st);
        rc |= nbStepSystem(ctx, st);
        if (nbStatusIsFatal(rc))   /* advance N-body system */
            return rc;

        rc |= nbCheckpoint(ctx, st);
        if (nbStatusIsFatal(rc))
            return rc;

        nbReportProgress(ctx, st, nbf->reportProgress);
    }

    if (BOINC_APPLICATION || ctx->checkpointT >= 0)
    {
        mw_report("Making final checkpoint\n");
        if (nbWriteCheckpoint(ctx, st))
        {
            mw_printf("Failed to write final checkpoint\n");
            return NBODY_CHECKPOINT_ERROR;
        }
    }

    return rc;
}

/* If possible, resume from a checkpoint. Otherwise do the necessary
 * initialization for a new run */
static NBodyStatus nbResumeOrNewRun(NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    if (nbResolveCheckpoint(st, nbf->checkpointFileName))
    {
        mw_printf("Failed to resolve checkpoint\n");
        return NBODY_ERROR;
    }

    /* If the checkpoint exists (and we want to use it), try to use it */
    if (nbf->ignoreCheckpoint || !nbResolvedCheckpointExists(st))
    {
        if (nbSetup(ctx, st, nbf))
        {
            mw_printf("Failed to read input parameters file\n");
            return NBODY_PARAM_FILE_ERROR;
        }
    }
    else /* Resume from checkpoint */
    {
        if (nbf->inputFile && !BOINC_APPLICATION)
        {
            mw_printf("Warning: input file '%s' unused\n", nbf->inputFile);
        }

        if (nbReadCheckpoint(ctx, st))
        {
            mw_report("Failed to read checkpoint\n");
            destroyNBodyState(st);
            return NBODY_CHECKPOINT_ERROR;
        }
        else
        {
            mw_report("Resumed from checkpoint '%s'\n", nbf->checkpointFileName);
        }
    }

    if (ctx->potentialType == EXTERNAL_POTENTIAL_CUSTOM_LUA)
    {
        /* We're using a custom potential, so we'll reevaluate the
         * script. We must do this once per thread.
         */
        if (nbOpenPotentialEvalStatePerThread(st, nbf))
        {
            return NBODY_PARAM_FILE_ERROR;
        }
    }

    return NBODY_SUCCESS;
}

/* Set context fields read from command line flags */
static void nbSetCtxFromFlags(NBodyCtx* ctx, const NBodyFlags* nbf)
{
    ctx->checkpointT = nbf->checkpointPeriod;
}

static void nbSetStateFromFlags(NBodyState* st, const NBodyFlags* nbf)
{
    st->ignoreResponsive = nbf->ignoreResponsive;
}

static void nbSetCLRequestFromFlags(CLRequest* clr, const NBodyFlags* nbf)
{
    clr->platform = nbf->platform;
    clr->devNum = nbf->devNum;
    clr->verbose = nbf->verbose;
    clr->reportProgress = nbf->reportProgress;
    clr->enableCheckpointing = FALSE;
    clr->enableProfiling = TRUE;
}

/* Try evaluating everything in the file to make sure it's OK */
int nbVerifyFile(const NBodyFlags* nbf)
{
    NBodyCtx ctx  = EMPTY_NBODYCTX;
    NBodyState st = EMPTY_NBODYSTATE;
    HistogramParams hp;

    if (nbSetup(&ctx, &st, nbf) || nbHistogramParamsCheck(nbf, &hp))
    {
        mw_printf("File failed\n");
        destroyNBodyState(&st);
        return FALSE;
    }
    else
    {
        mw_printf("File is OK\n");
        printNBodyCtx(&ctx);
        printHistogramParams(&hp);
        destroyNBodyState(&st);
        return TRUE;
    }
}


static int nbOutputIsUseful(const NBodyFlags* nbf)
{
    if (!nbf->outFileName && !nbf->histogramFileName && !nbf->histoutFileName && !nbf->verifyOnly)
    {
        mw_printf("Don't you want some kind of result?\n");
        return FALSE;
    }

    return TRUE;
}


/* FIXME */
static NBodyCtx _ctx = EMPTY_NBODYCTX;
static NBodyState _st = EMPTY_NBODYSTATE;

int nbMain(const NBodyFlags* nbf)
{
    NBodyCtx* ctx = &_ctx;
    NBodyState* st = &_st;
    CLRequest clr;

    NBodyStatus rc = NBODY_SUCCESS;
    double chisq;
    double ts = 0.0, te = 0.0;

    if (!nbOutputIsUseful(nbf))
    {
        return NBODY_USER_ERROR;
    }

    rc = nbResumeOrNewRun(ctx, st, nbf);
    if (nbStatusIsFatal(rc))
    {
        return rc;
    }

    nbSetCtxFromFlags(ctx, nbf); /* Do this after setup to avoid the setup clobbering the flags */
    nbSetStateFromFlags(st, nbf);
    nbSetCLRequestFromFlags(&clr, nbf);

    if (nbCreateSharedScene(st, ctx))
    {
        mw_printf("Failed to create shared scene\n");
    }

    if (nbf->reportProgress)
    {
        nbSetupCursesOutput();
    }

    ts = mwGetTime();
  #if NBODY_OPENCL
    if (nbf->noCL)
    {
        rc = nbRunSystem(ctx, st, nbf);
    }
    else
    {
        rc = initCLNBodyState(st, ctx, &clr);
        if (!nbStatusIsFatal(rc))
        {
            rc = nbRunSystemCL(ctx, st, nbf);
        }
    }
  #else
    rc = nbRunSystem(ctx, st, nbf);
  #endif /* NBODY_OPENCL */
    te = mwGetTime();

    if (nbf->reportProgress)
    {
        nbCleanupCursesOutput();
    }

    if (nbf->printTiming)
    {
        printf("<run_time> %f </run_time>\n", te - ts);
        fflush(stdout); /* Odd things happen with the OpenCL one where stdout starts disappearing */
    }


    if (nbStatusIsFatal(rc))
    {
        mw_printf("Error running system: %s (%d)\n", showNBodyStatus(rc), rc);
        destroyNBodyState(st);
        return rc;
    }

    if (nbStatusIsWarning(rc))
    {
        mw_printf("System complete with warnings: %s (%d)\n", showNBodyStatus(rc), rc);
    }


    if (nbf->histogramFileName || nbf->histoutFileName)  /* We want to match or produce a histogram */
    {
        /* Get the likelihood */
        chisq = nbChisq(ctx, st, nbf);
    }

    if (nbf->histogramFileName) /* The likelihood only means something when matching a histogram */
    {
        if (isnan(chisq))
        {
            mw_printf("Failed to calculate chisq\n");
            rc = NBODY_ERROR;
        }

        mw_printf("<search_likelihood>%.15f</search_likelihood>\n", chisq);
    }

    if (nbf->printBodies)
    {
        nbWriteBodies(ctx, st, nbf);
    }

    destroyNBodyState(st);

    return rc;
}

