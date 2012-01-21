/*
 * Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
 * Copyright (c) 2010, 2011 Matthew Arsenault
 * Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute.
 * Copyright (c) 2002-2006 John M. Fregeau, Richard Campbell, Jeff Molofee
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody.h"
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "nbody_show.h"
#include "nbody_lua.h"
#include "nbody_curses.h"
#include "nbody_shmem.h"
#include "nbody_defaults.h"
#include "nbody_plain.h"
#include "nbody_chisq.h"

#if NBODY_OPENCL
  #include "nbody_cl.h"
#endif

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
        if (!nbf->inputFile)
        {
            mw_printf("No input file and no checkpoint\n");
            return NBODY_USER_ERROR;
        }

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
    st->reportProgress = nbf->reportProgress;
    st->ignoreResponsive = nbf->ignoreResponsive;
}

static void nbSetCLRequestFromFlags(CLRequest* clr, const NBodyFlags* nbf)
{
    clr->platform = nbf->platform;
    clr->devNum = nbf->devNum;
    clr->verbose = nbf->verbose;
    clr->enableCheckpointing = FALSE;
    clr->enableProfiling = TRUE;
}

/* Try to run a potential function and see if it fails. Return TRUE on failure. */
static int nbVerifyPotentialFunction(const NBodyFlags* nbf, const NBodyCtx* ctx, NBodyState* st)
{
    mwvector acc;
    mwvector pos = mw_vec(1.0, 1.0, 0.0);

    if (ctx->potentialType != EXTERNAL_POTENTIAL_CUSTOM_LUA)
    {
        return FALSE;
    }

    /* Try to use it once to make sure it is OK */
    if (nbOpenPotentialEvalStatePerThread(st, nbf))
    {
        return TRUE;
    }

    nbEvalPotentialClosure(st, pos, &acc);
    return st->potentialEvalError;
}

/* Try evaluating everything in the file to make sure it's OK */
int nbVerifyFile(const NBodyFlags* nbf)
{
    NBodyCtx ctx  = EMPTY_NBODYCTX;
    NBodyState st = EMPTY_NBODYSTATE;
    HistogramParams hp;

    if (nbSetup(&ctx, &st, nbf) || nbHistogramParamsCheck(nbf, &hp) || nbVerifyPotentialFunction(nbf, &ctx, &st))
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
    if (   !nbf->outFileName
        && !nbf->histogramFileName
        && !nbf->histoutFileName
        && !nbf->printHistogram
        && !nbf->verifyOnly
        && !nbf->printTiming)
    {
        mw_printf("Don't you want some kind of result?\n");
        return FALSE;
    }

    return TRUE;
}

NBodyStatus nbStepSystem(const NBodyCtx* ctx, NBodyState* st)
{
  #if NBODY_OPENCL
    if (st->usesCL)
    {
        return nbStepSystemCL(ctx, st);
    }
  #endif

    return nbStepSystemPlain(ctx, st);
}

NBodyStatus nbRunSystem(const NBodyCtx* ctx, NBodyState* st)
{
  #if NBODY_OPENCL
    if (st->usesCL)
    {
        return nbRunSystemCL(ctx, st);
    }
  #endif

    return nbRunSystemPlain(ctx, st);
}

/* Output appropriate things depending on whether raw output, a
 * histogram, or just a likelihood is wanted.
 */
static NBodyStatus nbReportResults(const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf)
{
    NBodyHistogram* data = NULL;
    NBodyHistogram* histogram = NULL;
    double likelihood = NAN;
    NBodyLikelihoodMethod method;

    /* The likelihood only means something when matching a histogram */
    mwbool calculateLikelihood = (nbf->histogramFileName != NULL);

    if (nbf->outFileName)
    {
        nbWriteBodies(ctx, st, nbf);
    }

    if (calculateLikelihood || nbf->histoutFileName || nbf->printHistogram)
    {
        HistogramParams hp;

        if (nbGetLikelihoodInfo(nbf, &hp, &method) || method == NBODY_INVALID_METHOD)
        {
            mw_printf("Failed to get likelihood information\n");
            return NBODY_LIKELIHOOD_ERROR;
        }

        histogram = nbCreateHistogram(ctx, st, &hp);
        if (!histogram)
        {
            mw_printf("Failed to create histogram\n");
            return NBODY_LIKELIHOOD_ERROR;
        }
    }

    /* We want to write something whether or not the likelihood can be
     * calculated (i.e. given a histogram) so write this first */
    if (nbf->histoutFileName)
    {
        nbWriteHistogram(nbf->histoutFileName, ctx, st, histogram);
    }

    if (nbf->printHistogram)
    {
        nbPrintHistogram(DEFAULT_OUTPUT_FILE, histogram);
    }

    if (calculateLikelihood)   /* We want to match or produce a histogram */
    {
        data = nbReadHistogram(nbf->histogramFileName);
        if (!data)
        {
            free(histogram);
            return NBODY_LIKELIHOOD_ERROR;
        }

        likelihood = nbSystemChisq(st, data, histogram, method);
    }

    free(histogram);
    free(data);

    if (calculateLikelihood)
    {
        /* Reported negated distance since the search maximizes this */
        mw_printf("<search_likelihood>%.15f</search_likelihood>\n", -likelihood);
        if (isnan(likelihood))
        {
            mw_printf("Failed to calculate likelihood\n");
            return NBODY_LIKELIHOOD_ERROR;
        }
    }

    return NBODY_SUCCESS;
}

static NBodyCtx _ctx = EMPTY_NBODYCTX;
static NBodyState _st = EMPTY_NBODYSTATE;

int nbMain(const NBodyFlags* nbf)
{
    NBodyCtx* ctx = &_ctx;
    NBodyState* st = &_st;
    CLRequest clr;

    NBodyStatus rc = NBODY_SUCCESS;
    double ts = 0.0, te = 0.0;

    if (!nbOutputIsUseful(nbf))
    {
        return NBODY_USER_ERROR;
    }

    rc = nbResumeOrNewRun(ctx, st, nbf);
    if (nbStatusIsFatal(rc))
    {
        destroyNBodyState(st);
        return rc;
    }

    nbSetCtxFromFlags(ctx, nbf); /* Do this after setup to avoid the setup clobbering the flags */
    nbSetStateFromFlags(st, nbf);
    nbSetCLRequestFromFlags(&clr, nbf);

    if (nbCreateSharedScene(st, ctx))
    {
        mw_printf("Failed to create shared scene\n");
    }

    if (nbf->visualizer)
    {
        nbLaunchVisualizer(st, nbf->visArgs);
    }

  #if NBODY_OPENCL
    if (!nbf->noCL)
    {
        rc = nbInitNBodyStateCL(st, ctx, &clr);
        if (nbStatusIsFatal(rc))
        {
            destroyNBodyState(st);
            return rc;
        }
    }
  #endif /* NBODY_OPENCL */

    if (nbf->reportProgress)
    {
        nbSetupCursesOutput();
    }

    ts = mwGetTime();
    rc = nbRunSystem(ctx, st);
    te = mwGetTime();

    if (nbf->reportProgress)
    {
        nbCleanupCursesOutput();
    }

    if (nbStatusIsFatal(rc))
    {
        mw_printf("Error running system: %s (%d)\n", showNBodyStatus(rc), rc);
        destroyNBodyState(st);
        return rc;
    }
    else
    {
        if (nbStatusIsWarning(rc))
        {
            mw_printf("System complete with warnings: %s (%d)\n", showNBodyStatus(rc), rc);
        }

        if (nbf->printTiming)
        {
            printf("<run_time> %f </run_time>\n", te - ts);
        }
    }

    rc = nbReportResults(ctx, st, nbf);

    destroyNBodyState(st);

    return rc;
}

