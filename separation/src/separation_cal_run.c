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
#include "show_cal_types.h"
#include "separation_cal_run.h"
#include "separation_cal_setup.h"
#include "separation_cal_types.h"
#include "separation_cal_kernelgen.h"


static CALresult runKernel(MWCALInfo* ci, SeparationCALMem* cm, const CALdomain* domain)
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

    while (calCtxIsEventDone(ci->calctx, ev) == CAL_RESULT_PENDING);

    return CAL_RESULT_OK;
}
static CALresult setNuKernelArgs(MWCALInfo* ci,
                                 SeparationCALMem* cm,
                                 SeparationCALNames* cn,
                                 const IntegralArea* ia,
                                 CALuint nuStep)
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

static real readResults(MWCALInfo* ci,
                        SeparationCALMem* cm,
                        const IntegralArea* ia,
                        real* probs_results,
                        CALuint numberStreams)
{
    CALuint i;
    real result;

    result = sumResults(&cm->outBg, ia);

    for (i = 0; i < numberStreams; ++i)
        probs_results[i] = sumResults(&cm->outStreams[i], ia);

    return result;
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

static void printChunks(const SeparationCALChunks* chunks)
{
    warn("Using { %u, %u } chunk(s) of size { %u, %u }\n",
         chunks->nChunkR, chunks->nChunkMu,
         chunks->chunkSizeR, chunks->chunkSizeMu);
}

/* TODO: Actually do this */
static CALresult findCALChunks(const MWCALInfo* ci, const IntegralArea* ia, SeparationCALChunks* chunks)
{
    CALresult err = CAL_RESULT_OK;
    CALuint nChunk = 1;

    if (!mwEven(nChunk) && nChunk != 1)
    {
        warn("Number of chunks (%u) not even (or 1)\n", nChunk);
        return CAL_RESULT_ERROR;
    }

    if (nChunk == 1)
    {
        chunks->nChunkR = 1;
        chunks->nChunkMu = 1;
    }
    else
    {
        chunks->nChunkR = nChunk / 2;
        chunks->nChunkMu = nChunk / 2;
    }


    warn("HERPY: %u \n",
         8 * ci->devAttribs.wavefrontSize * ci->devAttribs.numberOfSIMD
        );

    warn("DERP %u\n", ia->mu_steps * ia->r_steps / (128 * ci->devAttribs.wavefrontSize * ci->devAttribs.numberOfSIMD));

    //chunks->nChunkR = 2;
    //chunks->nChunkMu = 8;

    //chunks->nChunkR = 10;
    //chunks->nChunkMu = 10;

    chunks->chunkSizeR = ia->r_steps / chunks->nChunkR;
    chunks->chunkSizeMu = ia->mu_steps / chunks->nChunkMu;

    printChunks(chunks);
    err |= chunkDivCheck("r steps", ia->r_steps, chunks->nChunkR);
    err |= chunkDivCheck("mu steps", ia->mu_steps, chunks->nChunkMu);

    return err;
}

static inline CALuint runNuStep(MWCALInfo* ci,
                                SeparationCALMem* cm,
                                SeparationCALNames* cn,
                                const IntegralArea* ia,
                                const SeparationCALChunks* chunks,
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
            err = setNuKernelArgs(ci, cm, cn, ia, nuStep);
            if (err != CAL_RESULT_OK)
                break;

            err = runKernel(ci, cm, &domain);
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

static real runIntegral(const AstronomyParameters* ap,
                        const IntegralArea* ia,
                        EvaluationState* es,
                        MWCALInfo* ci,
                        SeparationCALMem* cm,
                        real* probs_results)
{
    CALresult err;
    SeparationCALNames cn;
    double t1, t2;
    SeparationCALChunks chunks;

    if (findCALChunks(ci, ia, &chunks) != CAL_RESULT_OK)
        return NAN;

    memset(&cn, 0, sizeof(SeparationCALNames));
    err = getModuleNames(ci, &cn, cm->numberStreams);
    if (err != CAL_RESULT_OK)
    {
        cal_warn("Failed to get module names", err);
        return NAN;
    }

    err = setKernelArguments(ci, cm, &cn);
    if (err != CAL_RESULT_OK)
    {
        destroyModuleNames(&cn);
        return NAN;
    }

    for (es->nu_step = 0; es->nu_step < ia->nu_steps; es->nu_step++)
    {
        t1 = mwGetTime();

        err = runNuStep(ci, cm, &cn, ia, &chunks, es->nu_step);
        if (err != CAL_RESULT_OK)
            break;

        t2 = mwGetTime();
        reportProgress(ap, ia, es, es->nu_step + 1, 1000.0 * (t2 - t1));
    }

    destroyModuleNames(&cn);

    return (err != CAL_RESULT_OK) ? NAN : readResults(ci, cm, ia, probs_results, cm->numberStreams);
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


real integrateCAL(const AstronomyParameters* ap,
                  const IntegralArea* ia,
                  const StreamConstants* sc,
                  const StreamGauss sg,
                  real* st_probs,
                  EvaluationState* es,
                  const CLRequest* clr,
                  MWCALInfo* ci)
{
    CALresult err;
    SeparationCALMem cm;
    CALSeparationSizes sizes;
    real result = NAN;

    calculateCALSeparationSizes(&sizes, ap, ia);

    memset(&cm, 0, sizeof(SeparationCALMem));
    err = createSeparationBuffers(ci, &cm, ap, ia, sc, sg, &sizes);
    if (err != CAL_RESULT_OK)
        return NAN;

    result = runIntegral(ap, ia, es, ci, &cm, st_probs);

    err = releaseSeparationBuffers(ci, &cm);
    if (err != CAL_RESULT_OK)
        result = NAN;

    return result;
}

