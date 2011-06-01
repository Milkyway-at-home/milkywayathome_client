/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "separation_types.h"
#include "milkyway_util.h"
#include "setup_cl.h"
#include "milkyway_cl.h"
#include "mw_cl.h"
#include "separation_cl_buffers.h"
#include "separation_cl_defs.h"
#include "calculated_constants.h"
#include "run_cl.h"
#include "r_points.h"
#include "integrals.h"

static void sumStreamResults(real* streamResults,
                             const real* mappedResults,
                             const cl_uint mu_steps,
                             const cl_uint r_steps,
                             const cl_uint number_streams)
{
    cl_uint i, j, k, idx;
    Kahan* streamSums;

    streamSums = (Kahan*) mwCallocA(number_streams, sizeof(Kahan));

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            idx = (i * r_steps * number_streams) + (j * number_streams);
            for (k = 0; k < number_streams; ++k)
            {
                KAHAN_ADD(streamSums[k], mappedResults[idx + k]);
            }
        }
    }

    for (i = 0; i < number_streams; ++i)
        streamResults[i] = streamSums[i].sum + streamSums[i].correction;

    mwFreeA(streamSums);
}

static cl_int enqueueIntegralKernel(CLInfo* ci,
                                    cl_bool useDefault,
                                    const size_t offset[],
                                    const size_t global[],
                                    const size_t local[])
{
    cl_int err;

    err = clEnqueueNDRangeKernel(ci->queue,
                                 ci->kern,
                                 2,
                                 offset, global, useDefault ? NULL : local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error enqueueing integral kernel execution", err);
        return err;
    }

    return CL_SUCCESS;
}

/* The mu and r steps for a nu step were done in parallel. After that
 * we need to add the result to the running total + correction from
 * all the nu steps */
static real sumBgResults(const real* bgResults,
                         const cl_uint mu_steps,
                         const cl_uint r_steps)
{
    cl_uint i, j;
    Kahan bg_prob = ZERO_KAHAN;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            KAHAN_ADD(bg_prob, bgResults[r_steps * i + j]);
        }
    }

    return bg_prob.sum + bg_prob.correction;
}

static cl_int setNuKernelArgs(CLInfo* ci, const IntegralArea* ia, const cl_uint nu_step)
{
    cl_int err;
    NuId nuid;

    /* Avoid doing any trig in the broken ATI math. Also trig seems to
     * be more expensive there. Not doing the coordinate conversion
     * there also halves number of required registers, which prevents
     * enough threads to hide the horrible latency of the other
     * required reads. */
    nuid = calcNuStep(ia, nu_step);
    err = clSetKernelArg(ci->kern, 10, sizeof(real), &nuid.id);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error setting nu_id argument for step %u", err, nu_step);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int readKernelResults(CLInfo* ci,
                                SeparationCLMem* cm,
                                EvaluationState* es,
                                const IntegralArea* ia,
                                const cl_uint number_streams)
{
    cl_int err;
    real* bgResults;
    real* streamsTmp;
    size_t resultSize, streamsSize;

    resultSize = sizeof(real) * ia->mu_steps * ia->r_steps;
    bgResults = mapIntegralResults(ci, cm, resultSize);
    if (!bgResults)
    {
        warn("Failed to map integral results\n");
        return MW_CL_ERROR;
    }

    es->bgTmp = sumBgResults(bgResults, ia->mu_steps, ia->r_steps);

    err = clEnqueueUnmapMemObject(ci->queue, cm->outBg, bgResults, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to unmap results buffer", err);
        return err;
    }

    streamsSize = sizeof(real) * ia->mu_steps * ia->r_steps * number_streams;
    streamsTmp = mapStreamsResults(ci, cm, streamsSize);
    if (!streamsTmp)
    {
        warn("Failed to map stream results\n");
        return MW_CL_ERROR;
    }

    sumStreamResults(es->streamTmps, streamsTmp, ia->mu_steps, ia->r_steps, number_streams);

    err = clEnqueueUnmapMemObject(ci->queue, cm->outStreams, streamsTmp, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to unmap streams buffer", err);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int runNuStep(CLInfo* ci, const IntegralArea* ia, const RunSizes* runSizes, const cl_uint nu_step)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;
    size_t offset[2];
    size_t global[2];

    err = setNuKernelArgs(ci, ia, nu_step);
    if (err != CL_SUCCESS)
    {
        warn("Failed to set nu kernel argument\n");
        return err;
    }

    offset[1] = nu_step;
    global[1] = 1;
    for (i = 0; i < runSizes->nChunk && err == CL_SUCCESS; ++i)
    {
        offset[0] = runSizes->chunkBorders[i];
        global[0] = runSizes->chunkBorders[i + 1] - runSizes->chunkBorders[i];

        mw_begin_critical_section();
        err = enqueueIntegralKernel(ci, mwDivisible(global[0], runSizes->local[0]) ,
                                    offset, global, runSizes->local);
        if (err == CL_SUCCESS)
        {
            /* Give the screen a chance to redraw */
            err = clFinish(ci->queue);
            if (err != CL_SUCCESS)
                mwCLWarn("Failed to finish", err);
        }

        mw_end_critical_section();
    }

    return err;
}

static inline void reportProgress(const AstronomyParameters* ap,
                                  const IntegralArea* ia,
                                  EvaluationState* es,
                                  cl_uint step,
                                  double dt)
{
  #if BOINC_APPLICATION
    cl_uint prog;
    prog = es->current_calc_probs + ia->mu_steps * ia->r_steps * step;
    boinc_fraction_done((double) prog / ap->total_calc_probs);
  #else
    printf("Step %u: %fms\n", step, dt);
  #endif /* BOINC_APPLICATION */
}

static cl_int checkpointCL(CLInfo* ci, SeparationCLMem* cm, const IntegralArea* ia, EvaluationState* es)
{
    cl_int err;

    err = readKernelResults(ci, cm, es, ia, es->numberStreams);
    if (err != CL_SUCCESS)
        return err;

    err = writeCheckpoint(es) ? MW_CL_ERROR : CL_SUCCESS;

  #if BOINC_APPLICATION
    boinc_checkpoint_completed();
  #endif

    return err;
}

static cl_int runIntegral(CLInfo* ci,
                          SeparationCLMem* cm,
                          RunSizes* runSizes,
                          EvaluationState* es,
                          const CLRequest* clr,
                          const AstronomyParameters* ap,
                          const IntegralArea* ia)
{
    cl_int err = CL_SUCCESS;
    double t1, t2, dt;
    double tAcc = 0.0;

    for (; es->nu_step < ia->nu_steps; es->nu_step++)
    {
        if (clr->enableCheckpointing && timeToCheckpointGPU(es, ia))
        {
            err = checkpointCL(ci, cm, ia, es);
            if (err != CL_SUCCESS)
                break;
        }

        t1 = mwGetTimeMilli();
        err = runNuStep(ci, ia, runSizes, es->nu_step);
        if (err != CL_SUCCESS)
        {
            mwCLWarn("Failed to run nu step", err);
            return err;
        }
        t2 = mwGetTimeMilli();

        dt = t2 - t1;
        tAcc += dt;

        reportProgress(ap, ia, es, es->nu_step + 1, dt);
    }

    es->nu_step = 0;

    warn("Integration time: %f s. Average time per iteration = %f ms\n",
         tAcc / 1000.0, tAcc / (double) ia->nu_steps);

    if (err == CL_SUCCESS)
    {
        err = readKernelResults(ci, cm, es, ia, ap->number_streams);
        if (err != CL_SUCCESS)
            warn("Failed to read final kernel results\n");

        addTmpSums(es); /* Add final episode to running totals */
    }

    return err;
}

cl_int integrateCL(const AstronomyParameters* ap,
                   const IntegralArea* ia,
                   const StreamConstants* sc,
                   const StreamGauss sg,
                   EvaluationState* es,
                   const CLRequest* clr,
                   CLInfo* ci,
                   cl_bool useImages)
{
    cl_int err;
    RunSizes runSizes;
    SeparationSizes sizes;
    SeparationCLMem cm = EMPTY_SEPARATION_CL_MEM;

    /* Need to test sizes for each integral, since the area size can change */
    calculateSizes(&sizes, ap, ia);

    if (findRunSizes(&runSizes, ci, &ci->di, ia, clr))
    {
        warn("Failed to find good run sizes\n");
        return MW_CL_ERROR;
    }

    err = createSeparationBuffers(ci, &cm, ap, ia, sc, sg, &sizes, useImages);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to create CL buffers", err);
        return err;
    }

    err = separationSetKernelArgs(ci, &cm, &runSizes);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to set integral kernel arguments", err);
        return err;
    }

    err = runIntegral(ci, &cm, &runSizes, es, clr, ap, ia);

    releaseSeparationBuffers(&cm);

    separationIntegralApplyCorrection(es);

    return err;
}

