/*
 * Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
 * Copyright (c) 2010, 2011 Matthew Arsenault
 * Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute.
 * Copyright (c) 2002-2006 John M. Fregeau, Richard Campbell, Jeff Molofee
 * Copyright (c) 2016-2018 Siddhartha Shelton
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
#include "nbody_likelihood.h"
#include "nbody_histogram.h"
#include "nbody_types.h"
#include "nbody_autodiff.h"

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
    memset(clr, 0, sizeof(*clr));

    clr->platform = nbf->platform;
    clr->devNum = nbf->devNum;
    clr->verbose = nbf->verbose;
    clr->enableCheckpointing = !nbf->disableGPUCheckpointing;
    clr->enableProfiling = TRUE;
    clr->pollingMode = MW_POLL_CL_WAIT_FOR_EVENTS;
}

/* Try to run a potential function and see if it fails. Return TRUE on failure. */
static int nbVerifyPotentialFunction(const NBodyFlags* nbf, const NBodyCtx* ctx, NBodyState* st)
{
    mwvector acc;
    mwvector pos = mw_vec(mw_real_const(1.0), mw_real_const(1.0), ZERO_REAL);

    if (ctx->potentialType != EXTERNAL_POTENTIAL_CUSTOM_LUA)
    {
        return FALSE;
    }

    /* Try to use it once to make sure it is OK */
    if (nbOpenPotentialEvalStatePerThread(st, nbf))
    {
        return TRUE;
    }

    nbEvalPotentialClosure(st, &pos, &acc);
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
    if(!ctx->LMC)
    {
        mwvector zero;
        SET_VECTOR(&zero,ZERO_REAL,ZERO_REAL,ZERO_REAL);
        return nbStepSystemPlain(ctx, st, &zero, &zero); 
    }
    else
    {
        return nbStepSystemPlain(ctx, st, &st->shiftByLMC[st->step], &st->shiftByLMC[st->step+1]);
    }
}

NBodyStatus nbRunSystem(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
  #if NBODY_OPENCL
    if (st->usesCL)
    {
        return nbRunSystemCL(ctx, st);
    }
  #endif

    return nbRunSystemPlain(ctx, st, nbf);
}

/* Output appropriate things depending on whether raw output, a
 * histogram, or just a likelihood is wanted.
 */
static NBodyStatus nbReportResults(const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf)
{
    MainStruct* data = NULL;
    MainStruct* histogram = NULL;
    real likelihood = mw_real_const(NAN);
    real likelihood_EMD = mw_real_const(NAN);
    real likelihood_Mass = mw_real_const(NAN);
    real likelihood_Beta = mw_real_const(NAN);
    real likelihood_Vel = mw_real_const(NAN);
    real likelihood_VelAvg = mw_real_const(NAN);
    real likelihood_BetaAvg = mw_real_const(NAN);
    real likelihood_Dist = mw_real_const(NAN);
    NBodyLikelihoodMethod method;

    real *likelihoodArray;

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

        histogram = (nbCreateHistogram(ctx, st, &hp));
        if (!histogram)
        {
            mw_printf("Failed to create histogram\n");
            return NBODY_LIKELIHOOD_ERROR;
        }
    }

    /* We want to write something whether or not the likelihood can be
     * calculated (i.e. given a histogram) so write this first */
    /*
     * this has been moved to the best likelihood calculation in nbody_plain.c
     */
    if (nbf->histoutFileName && !calculateLikelihood)
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
        
        
        likelihoodArray    = nbSystemLikelihood(st, data, histogram, method);
        likelihood         = likelihoodArray[0];
        likelihood_EMD     = likelihoodArray[1];
        likelihood_Mass    = likelihoodArray[2];
        likelihood_Beta    = likelihoodArray[3];
        likelihood_Vel     = likelihoodArray[4];
        likelihood_BetaAvg = likelihoodArray[5];
        likelihood_VelAvg  = likelihoodArray[6];
        likelihood_Dist    = likelihoodArray[7];

        /*
          Used to fix Windows platform issues.  Windows' infinity is expressed as:
          1.#INF00000, -1.#INF00000, or 0.#INF000000.  The server reads these as -1, 1, and 0
          respectively, accounting for the sign change.  Thus, I have changed overflow
          infinities (not errors) to be the worst case.  The worst case is now the actual
          worst thing that can happen.
        */
        /*
        * It previous returned the worse case when the likelihood == 0. 
        * Changed it to be best case, 1e-9 which has been added in nbody_defaults.h
        */
        if (showRealValue(&likelihood) > DEFAULT_WORST_CASE || showRealValue(&likelihood) < (-1 * DEFAULT_WORST_CASE) )
        {
            mw_printf("Poor likelihood.  Returning worst case.\n");
            setRealValue(&likelihood, DEFAULT_WORST_CASE);
        }
        else if(showRealValue(&likelihood) == 0.0)
        {
            setRealValue(&likelihood, DEFAULT_BEST_CASE);
        }
        
        
        /* if the end state likelihood is not better than the best likelihood then
         * replace it with the best likelihood 
         */
        /* only do the comparison if we are using the best likelihood code */
        if(mw_fabs_0(showRealValue(&likelihood)) > mw_fabs_0(showRealValue(&st->bestLikelihood)) && ctx->useBestLike)
        {
            likelihood         = st->bestLikelihood;
            likelihood_EMD     = st->bestLikelihood_EMD;
            likelihood_Mass    = st->bestLikelihood_Mass;
            likelihood_Beta    = st->bestLikelihood_Beta;
            likelihood_Vel     = st->bestLikelihood_Vel;
            likelihood_VelAvg  = st->bestLikelihood_VelAvg;
            likelihood_BetaAvg = st->bestLikelihood_BetaAvg;
            likelihood_Dist    = st->bestLikelihood_Dist;
        }
        else
        {
            /* If the end state likelihood is worse than best likelihood then replace it, and keep the hist that would have been written.
             * If there was never an improvement, i.e. the likelihood stayed worse case the entire time, then 
             * the likelihood each timestep and best like are always equal. Then the best like code would not have written a hist. 
             * if the best like is worse than the end state likelihood, or if they are equal,  rewrite the hist.
             */
            if(nbf->histoutFileName)
            {
                nbWriteHistogram(nbf->histoutFileName, ctx, st, histogram);
            }
        }        
    }
    
    free(histogram);
    free(data);

  if (calculateLikelihood)
    {
        /* Reported negated distance since the search maximizes this */
      if (isnan(showRealValue(&likelihood)))
        {
            setRealValue(&likelihood, DEFAULT_WORST_CASE);
            mw_printf("Likelihood was NAN. Returning worst case. \n");
            //mw_printf("<search_likelihood>%.15f</search_likelihood>\n", -showRealValue(&likelihood));
            //return NBODY_SUCCESS;
        }
        mw_printf("<search_likelihood>%.15f</search_likelihood>\n", -showRealValue(&likelihood));
        //printRealGradient(&likelihood, "likelihood");

        mw_printf("<search_likelihood_EMD>%.15f</search_likelihood_EMD>\n", -showRealValue(&likelihood_EMD));
        //printRealGradient(&likelihood_EMD, "likelihood_EMD");

        mw_printf("<search_likelihood_Mass>%.15f</search_likelihood_Mass>\n", -showRealValue(&likelihood_Mass));
        //printRealGradient(&likelihood_Mass, "likelihood_Mass");

	if (st->useBetaDisp)
        {
            mw_printf("<search_likelihood_Beta>%.15f</search_likelihood_Beta>\n", -showRealValue(&likelihood_Beta));
            //printRealGradient(&likelihood_Beta, "likelihood_Beta");
        }
	if (st->useVelDisp)
        {
            mw_printf("<search_likelihood_Vel>%.15f</search_likelihood_Vel>\n", -showRealValue(&likelihood_Vel));
            //printRealGradient(&likelihood_Vel, "likelihood_Vel");
        }
       if (st->useBetaComp)
       {
           mw_printf("<search_likelihood_BetaAvg>%.15f</search_likelihood_BetaAvg>\n", -showRealValue(&likelihood_BetaAvg));
           //printRealGradient(&likelihood_BetaAvg, "likelihood_BetaAvg");
       }
       if (st->useVlos)
       {
           mw_printf("<search_likelihood_VelAvg>%.15f</search_likelihood_VelAvg>\n", -showRealValue(&likelihood_VelAvg));
           //printRealGradient(&likelihood_VelAvg, "likelihood_VelAvg");
       }
       if (st->useDist)
       {
           mw_printf("<search_likelihood_Dist>%.15f</search_likelihood_Dist>\n", -showRealValue(&likelihood_Dist));
           //printRealGradient(&likelihood_Dist, "likelihood_Dist");
       }

       /* Write out the derivative information from AUTODIFF if file available */
       if(nbf->autoDiffFileName)
       {
           nbWriteAutoDiff(nbf->autoDiffFileName, &likelihood);
       }
    }


    return NBODY_SUCCESS;
}

int nbMain(const NBodyFlags* nbf)
{
    NBodyCtx* ctx = &_ctx;
    NBodyState* st = &_st;
    CLRequest clr;

    NBodyStatus rc = NBODY_SUCCESS;
    real_0 ts = 0.0, te = 0.0;

    if (!nbOutputIsUseful(nbf))
    {
        return NBODY_USER_ERROR;
    }

    nbSetCLRequestFromFlags(&clr, nbf);

    /* Find out what device we're using so we can tell the workunit
     * about it */
    if (NBODY_OPENCL && !nbf->noCL)
    {
        rc = nbInitCL(st, ctx, &clr);
        if (nbStatusIsFatal(rc))
        {
            destroyNBodyState(st);
            return rc;
        }
    }

    rc = nbResumeOrNewRun(ctx, st, nbf);
    if (nbStatusIsFatal(rc))
    {
        destroyNBodyState(st);
        return rc;
    }

    NBodyState initialState = EMPTY_NBODYSTATE;
    //for the first run, just assume the best likelihood timestep will occur in middle of best-likelihood window
    //convert eff_best_like_start to the original best like start
    real_0 ogBestLikeStart = (2*ctx->BestLikeStart)/(ctx->BestLikeStart + 1);
    if(ctx->useBestLike){
        //assume evolve time has been adjusted to be the end of the best-likelihood window
        st->previousForwardTime = ctx->timeEvolve/(2 - ogBestLikeStart);
    }else{
        st->previousForwardTime = ctx->timeEvolve;
    }
    if(ctx->pot.disk2.type != OrbitingBar){
        ctx->calibrationRuns = 0;
    }
    //Run forward evolution calibrationRuns + 1 times
    for(int i = 0; i <= ctx->calibrationRuns; i++){
        //these for checkpointing
        nbSetCtxFromFlags(ctx, nbf); /* Do this after setup to avoid the setup clobbering the flags */
        nbSetStateFromFlags(st, nbf); 

        if (NBODY_OPENCL && !nbf->noCL)
        {
            rc = nbInitNBodyStateCL(st, ctx);
            if (nbStatusIsFatal(rc))
            {
                destroyNBodyState(st);
                return rc;
            }
        }

        if (nbCreateSharedScene(st, ctx))
        {
            mw_printf("Failed to create shared scene\n");
        }

        if (nbf->visualizer && st->scene)
        {
            /* Make sure the first scene is available for the launched graphics */
            nbForceUpdateDisplayedBodies(ctx, st);

            /* Launch graphics and make sure we are sure the graphics is
            * attached in case we are using blocking mode */
            nbLaunchVisualizer(st, nbf->graphicsBin, nbf->visArgs);
        }

        if (nbf->reportProgress)
        {
            nbSetupCursesOutput();
        }

        ts = mwGetTime();

        st->useVelDisp = ctx->useVelDisp;
        st->useBetaDisp = ctx->useBetaDisp;
        st->useBetaComp = ctx->useBetaComp;
        st->useVlos = ctx->useVlos;
        st->useDist = ctx->useDist;

        //save the state if about to start first calibration run
        if(ctx->calibrationRuns > 0 && i == 0){
            cloneNBodyState(&initialState, st);
        }

        rc = nbRunSystem(ctx, st, nbf);

        //debug output for calibration runs
        /*real expectedForwardTime = st->timeEvolve;
        if(i == 0){
            expectedForwardTime = ctx->timeBack;
        }
        mw_printf("run: %d forwardTime: %f\n", i, st->bestLikelihood_time);
        mw_printf("expected forward time - real forward time = %f\n\n", expectedForwardTime - st->bestLikelihood_time);
        */

        if(i < ctx->calibrationRuns){
            //grab the best likelihood time
            real_0 forwardTime = st->bestLikelihood_time;
            //reset the state for the next run
            *st = (NBodyState)EMPTY_NBODYSTATE;
            cloneNBodyState(st, &initialState);
            //set previous forward time for the next run
            st->previousForwardTime = forwardTime;
        }
        nbResolveCheckpoint(st, nbf->checkpointFileName);
    }

    te = mwGetTime();
    //mw_printf("After end GetTime\n");

    if (nbf->reportProgress)
    {
        nbCleanupCursesOutput();
    }
    //mw_printf("After reportProgress\n");

    nbReportSimulationComplete(st);
    //mw_printf("After nbReportSimulationComplete\n");

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
    //mw_printf("After Status Check\n");

    rc = nbReportResults(ctx, st, nbf);
    //mw_printf("After nbReportResults\n");

    destroyNBodyState(st);
    //mw_printf("After destroyNBodyState\n");
    destroyNBodyState(&initialState);

    return rc;
}
