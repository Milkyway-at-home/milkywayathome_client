/*
Copyright (C) 2010  Matthew Arsenault

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

#include "milkyway_util.h"
#include "separation_types.h"
#include "calculated_constants.h"
#include "r_points.h"
#include "integrals.h"
#include "show_cal_types.h"
#include "separation_cal_run.h"
#include "separation_cal_setup.h"
#include "separation_cal_types.h"
#include "separation_cal_kernelgen.h"

static CALresult waitForKernel(MWCALInfo* ci, CALevent ev, CALuint initialWait, CALint pollingMode)
{
    CALresult err;

    mwMilliSleep(initialWait);  /* Sleep for initial estimate before polling */

    if (pollingMode > 0)
    {
        while (calCtxIsEventDone(ci->calctx, ev) == CAL_RESULT_PENDING)
            mwMilliSleep(pollingMode);
    }
    else if (pollingMode == 0)
    {
        err = mw_calCtxWaitForEvents(ci->calctx, &ev, 1, 0);
        if (err != CAL_RESULT_OK)
        {
            cal_warn("Error waiting for kernel", err);
            return err;
        }
    }
    else /* Busy wait */
    {
        while (calCtxIsEventDone(ci->calctx, ev) == CAL_RESULT_PENDING);
    }

    return CAL_RESULT_OK;
}

static CALresult runKernel(MWCALInfo* ci, const CALdomain* domain, CALint pollingMode, CALuint initialWait)
{
    CALresult err;
    CALevent ev = 0;

#if 1
    err = calCtxRunProgram(&ev, ci->calctx, ci->func, domain);
#else
    CALdomain3D global = { domain->width, domain->height, 1 };
    CALdomain3D local = { 64, 28, 1 };
    CALprogramGrid grid;

    grid.func = ci->func;
    grid.flags = 0;

    grid.gridBlock = local;

    grid.gridSize.width  = (global.width + local.width - 1) / local.width;
    grid.gridSize.height = (global.height + local.height - 1) / local.height;
    grid.gridSize.depth  = (global.depth + local.depth - 1) / local.depth;

    #if 0
    warn("arst %u %u %u -> { %u %u }\n", grid.gridSize.width, grid.gridSize.height, grid.gridSize.depth,

         grid.gridSize.width * local.width,
         grid.gridSize.height * local.height
        );
    #endif

    err = calCtxRunProgramGrid(&ev, ci->calctx, &grid);
#endif

    if (err != CAL_RESULT_OK)
    {
        cal_warn("Error running kernel", err);
        return err;
    }

    return waitForKernel(ci, ev, initialWait, pollingMode);
}
static CALresult setNuKernelArgs(SeparationCALMem* cm, const IntegralArea* ia, CALuint nuStep)
{
    CALresult err;
    CALdouble* nuBufPtr;
    CALfloat* nuStepPtr;
    CALuint pitch = 0;
    NuId nuid;

    err = mapMWMemRes(&cm->nuBuf, (CALvoid**) &nuBufPtr, &pitch);
    if (err != CAL_RESULT_OK)
        return err;

    nuid = calcNuStep(ia, nuStep);

    nuStepPtr = (CALfloat*) nuBufPtr;
    *nuStepPtr = (CALfloat) nuStep;
    nuBufPtr[1] = nuid.id;

    err = unmapMWMemRes(&cm->nuBuf);
    if (err != CAL_RESULT_OK)
        return err;

    return CAL_RESULT_OK;
}

static real sumResults(MWMemRes* mr, const IntegralArea* ia)
{
    CALuint i, j, pitch;
    Kahan* bufPtr;
    Kahan* tmp;
    Kahan ksum = ZERO_KAHAN;
    CALresult err = CAL_RESULT_OK;

    err = mapMWMemRes(mr, (CALvoid**) &bufPtr, &pitch);
    if (err != CAL_RESULT_OK)
        return NAN;

    for (i = 0; i < ia->mu_steps; ++i)
    {
        tmp = &bufPtr[i * pitch];
        for (j = 0; j < ia->r_steps; ++j)
        {
            KAHAN_ADD(ksum, tmp[j].sum);
        }
    }

    err = unmapMWMemRes(mr);
    if (err != CAL_RESULT_OK)
        return NAN;

    return ksum.sum + ksum.correction;
}

static void readResults(SeparationCALMem* cm, const IntegralArea* ia, EvaluationState* es)
{
    CALuint i;

    es->bgTmp = sumResults(&cm->outBg, ia);

    for (i = 0; i < es->numberStreams; ++i)
    {
        es->streamTmps[i] = sumResults(&cm->outStreams[i], ia);
    }
}

static CALresult chunkDivCheck(const char* name, CALuint dim, CALuint nChunk)
{
    if (!mwDivisible(dim, nChunk))
    {
        warn("%s (%u) not divisible by n chunks (%u)\n", name, dim, nChunk);
        return CAL_RESULT_ERROR;
    }

    return CAL_RESULT_OK;
}

static void printChunks(const IntegralArea* ia, const SeparationCALChunks* chunks)
{
    warn("Integration range: { nu_steps = %u, mu_steps = %u, r_steps = %u }\n"
         "Using { %u, %u } chunk(s) of size { %u, %u }\n",
         ia->nu_steps, ia->mu_steps, ia->r_steps,
         chunks->nChunkR, chunks->nChunkMu,
         chunks->chunkSizeR, chunks->chunkSizeMu);
}

static CALuint64 estimateWUFLOPsPerIter(const AstronomyParameters* ap, const IntegralArea* ia)
{
    CALuint64 perItem, perIter;

    CALuint64 tmp = 32 + ap->number_streams * 57;
    if (ap->aux_bg_profile)
        tmp += 8;

    perItem = tmp * ap->convolve + 1 + (ap->number_streams * 2);
    perIter = perItem * ia->mu_steps * ia->r_steps;

    return perIter;
}

static CALuint64 deviceFlopsEstimate(const CALdeviceattribs* d)
{
    CALuint64 vliw, flops, doubleFrac;

    switch (d->target)
    {
        case CAL_TARGET_600:
        case CAL_TARGET_610:
        case CAL_TARGET_630:
        case CAL_TARGET_670:
        case CAL_TARGET_7XX:
        case CAL_TARGET_770:
        case CAL_TARGET_710:
        case CAL_TARGET_730:
        case CAL_TARGET_CYPRESS:
        case CAL_TARGET_JUNIPER:
        case CAL_TARGET_REDWOOD:
        case CAL_TARGET_CEDAR:
            vliw = 5;
            doubleFrac = 5;
            break;

        case CAL_TARGET_WRESTLER:
        case CAL_TARGET_CAYMAN:
        case CAL_TARGET_BARTS:
            vliw = 4;
            doubleFrac = 4;
            break;

        case CAL_TARGET_RESERVED0:
        case CAL_TARGET_RESERVED1:
        case CAL_TARGET_RESERVED2:
        default:
            doubleFrac = 5;
            vliw = 5;
            warn("Unknown target type: %s (%d)\n", showCALtargetEnum(d->target), d->target);
    }

    flops = 2 * (d->numberOfSIMD * vliw * 16) * d->engineClock * 1000000;

  #if DOUBLEPREC
    flops /= doubleFrac;
  #endif

    return flops;
}

static CALuint deviceChunkEstimate(const AstronomyParameters* ap,
                                   const IntegralArea* ia,
                                   const CALdeviceattribs* devAttribs,
                                   const CLRequest* clr,
                                   CALuint* chunkWaitEstimate)
{
    CALuint64 flops, iterFlops;
    CALdouble effFlops, estIterTime, timePerIter;
    CALuint nChunk;

    if (clr->responsivenessFactor <= 0.01)
        return 1;

    flops = deviceFlopsEstimate(devAttribs);
    effFlops = 0.8 * (CALdouble) flops;
    iterFlops = estimateWUFLOPsPerIter(ap, ia);

    estIterTime = 1000.0 * (CALdouble) iterFlops / effFlops; /* milliseconds */
    estIterTime *= clr->responsivenessFactor;

    timePerIter = 1000.0 / clr->targetFrequency;

    nChunk = (CALuint) (estIterTime / timePerIter);
    if (nChunk <= 0)
        nChunk = 1;

    if (nChunk >= ia->mu_steps)
        return ia->mu_steps;

    /* Round up to next evenly divisible integer */
    while (!mwDivisible(ia->mu_steps, nChunk) && nChunk <= ia->mu_steps)
        ++nChunk;

    //*chunkWaitEstimate = (CALuint) (0.1 * estIterTime / nChunk); /* Sleep for 50% of estimated time before polling */
    *chunkWaitEstimate = 0;
    warn("Estimated iteration time %f ms\n"
         "Target frequency %f Hz, polling mode %d, using responsiveness factor of %f\n"
         "Dividing into %u chunks, initially sleeping for %u ms\n",
         estIterTime, clr->targetFrequency, clr->pollingMode, clr->responsivenessFactor, nChunk, *chunkWaitEstimate);

    return nChunk;
}

static CALresult findCALChunks(const AstronomyParameters* ap,
                               const MWCALInfo* ci,
                               const CLRequest* clr,
                               const IntegralArea* ia,
                               SeparationCALChunks* chunks)
{
    CALresult err = CAL_RESULT_OK;
    CALuint nChunk;

    nChunk = deviceChunkEstimate(ap, ia, &ci->devAttribs, clr, &chunks->chunkWaitTime);
    if (nChunk == 0)
    {
        warn("Invalid number of chunks: %u\n", nChunk);
        return CAL_RESULT_ERROR;
    }
    else if (nChunk == 1)
    {
        chunks->nChunkR = 1;
        chunks->nChunkMu = 1;
    }
    else
    {
        chunks->nChunkR = 1;
        chunks->nChunkMu = nChunk;
    }

    chunks->chunkSizeR = ia->r_steps / chunks->nChunkR;
    chunks->chunkSizeMu = ia->mu_steps / chunks->nChunkMu;

    printChunks(ia, chunks);
    err |= chunkDivCheck("r steps", ia->r_steps, chunks->nChunkR);
    err |= chunkDivCheck("mu steps", ia->mu_steps, chunks->nChunkMu);

    return err;
}

static inline CALuint runNuStep(MWCALInfo* ci,
                                SeparationCALMem* cm,
                                const IntegralArea* ia,
                                const SeparationCALChunks* chunks,
                                CALint pollingMode,
                                CALuint nuStep)
{
    CALdomain domain;
    CALresult err = CAL_RESULT_OK;

    /* It's much faster to do chunking in both dimensions when using
     * images in tiled format */

    domain.width = chunks->chunkSizeR;
    domain.height = chunks->chunkSizeMu;

    for (domain.x = 0; domain.x < ia->r_steps; domain.x += chunks->chunkSizeR)
    {
        for (domain.y = 0; domain.y < ia->mu_steps; domain.y += chunks->chunkSizeMu)
        {
            err = setNuKernelArgs(cm, ia, nuStep);
            if (err != CAL_RESULT_OK)
                break;

            err = runKernel(ci, &domain, pollingMode, chunks->chunkWaitTime);
            if (err != CAL_RESULT_OK)
                goto run_step_exit;
        }
    }

run_step_exit:

    return err;
}

static inline void reportProgress(const AstronomyParameters* ap,
                                  const IntegralArea* ia,
                                  EvaluationState* es,
                                  CALuint step,
                                  double dt)
{
  #if BOINC_APPLICATION
    CALuint prog;
    prog = es->current_calc_probs + ia->mu_steps * ia->r_steps * step;
    boinc_fraction_done((double) prog / ap->total_calc_probs);
  #else
    printf("Step %u: %fms\n", step, dt);
  #endif /* BOINC_APPLICATION */
}

static CALresult checkpointCAL(SeparationCALMem* cm, const IntegralArea* ia, EvaluationState* es)
{
    CALresult err;

    readResults(cm, ia, es);
    err = writeCheckpoint(es) ? CAL_RESULT_ERROR : CAL_RESULT_OK;

  #if BOINC_APPLICATION
    boinc_checkpoint_completed();
  #endif

    return err;
}


static CALresult runIntegral(const AstronomyParameters* ap,
                             const IntegralArea* ia,
                             EvaluationState* es,
                             const CLRequest* clr,
                             MWCALInfo* ci,
                             SeparationCALMem* cm)
{
    CALresult err;
    SeparationCALNames cn;
    double t1, t2, dt, tAcc = 0.0;
    SeparationCALChunks chunks;

    if (findCALChunks(ap, ci, clr, ia, &chunks) != CAL_RESULT_OK)
        return CAL_RESULT_ERROR;

    memset(&cn, 0, sizeof(SeparationCALNames));
    err = getModuleNames(ci, &cn, cm->numberStreams);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get module names", err);
        return err;
    }

    err = setKernelArguments(ci, cm, &cn);
    if (err != CAL_RESULT_OK)
    {
        destroyModuleNames(&cn);
        return err;
    }

    for (; es->nu_step < ia->nu_steps; es->nu_step++)
    {
        if (clr->enableCheckpointing && timeToCheckpointGPU(es, ia))
        {
            err = checkpointCAL(cm, ia, es);
            if (err != CAL_RESULT_OK)
                break;
        }

        t1 = mwGetTimeMilli();

        err = runNuStep(ci, cm, ia, &chunks, clr->pollingMode, es->nu_step);
        if (err != CAL_RESULT_OK)
            break;

        t2 = mwGetTimeMilli();
        dt = t2 - t1;
        tAcc += dt;

        reportProgress(ap, ia, es, es->nu_step + 1, dt);
    }

    es->nu_step = 0;

    warn("Integration time = %f s, average per iteration = %f ms\n", 1.0e-3 * tAcc, tAcc / ia->nu_steps);

    destroyModuleNames(&cn);

    if (err == CAL_RESULT_OK)
    {
        readResults(cm, ia, es);
        addTmpSums(es); /* Add final episode to running totals */
    }

    return err;
}

static void calculateCALSeparationSizes(CALSeparationSizes* sizes,
                                        const AstronomyParameters* ap,
                                        const IntegralArea* ia)
{
    sizes->outBg = sizeof(Kahan) * ia->mu_steps * ia->r_steps;
    sizes->outStreams = sizeof(Kahan) * ia->mu_steps * ia->r_steps * ap->number_streams;
    sizes->rPts = sizeof(RPoints) * ap->convolve * ia->r_steps;
    sizes->rc = sizeof(RConsts) * ia->r_steps;
    sizes->sg_dx = sizeof(real) * ap->convolve;
    sizes->lTrig = sizeof(LTrigPair) * ia->mu_steps * ia->nu_steps;
    sizes->bTrig = sizeof(real) * ia->mu_steps * ia->nu_steps;

    sizes->nuSteps = ia->nu_steps;
    sizes->muSteps = ia->mu_steps;
    sizes->rSteps = ia->r_steps;
}


CALresult integrateCAL(const AstronomyParameters* ap,
                       const IntegralArea* ia,
                       const StreamGauss sg,
                       EvaluationState* es,
                       const CLRequest* clr,
                       MWCALInfo* ci)
{
    CALresult err;
    SeparationCALMem cm;
    CALSeparationSizes sizes;

    calculateCALSeparationSizes(&sizes, ap, ia);

    memset(&cm, 0, sizeof(cm));
    err = createSeparationBuffers(ci, &cm, ap, ia, sg, &sizes);
    if (err != CAL_RESULT_OK)
        return err;

    err = runIntegral(ap, ia, es, clr, ci, &cm);

    err |= releaseSeparationBuffers(ci, &cm);
    if (err != CAL_RESULT_OK)
        return err;

    separationIntegralApplyCorrection(es);

    return CAL_RESULT_OK;
}

