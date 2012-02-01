/*
 *  Copyright (c) 2010-2012 Matthew Arsenault
 *  Copyright (c) 2010-2012 Rensselaer Polytechnic Institute
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "separation_types.h"
#include "milkyway_util.h"
#include "setup_cl.h"
#include "milkyway_cl.h"
#include "separation_cl_buffers.h"
#include "calculated_constants.h"
#include "run_cl.h"
#include "r_points.h"
#include "integrals.h"

static void swapBuffers(cl_mem bufs[2])
{
    cl_mem tmp = bufs[0];
    bufs[0] = bufs[1];
    bufs[1] = tmp;
}

static cl_int runSummarization(CLInfo* ci,
                               SeparationCLMem* cm,
                               const IntegralArea* ia,
                               cl_uint which,
                               Kahan* resultOut)
{
    cl_int err = CL_SUCCESS;
    cl_mem buf;
    cl_uint offset;
    size_t global[1];
    size_t local[1];
    real result[2] = { -1.0, -1.0 };
    cl_uint nElements = ia->r_steps * ia->mu_steps;
    cl_mem sumBufs[2] = { cm->summarizationBufs[0], cm->summarizationBufs[1] };

    if (which == 0)
    {
        buf = cm->outBg;
        offset = 0;
    }
    else
    {
        buf = cm->outStreams;
        offset = (which - 1) * nElements;
    }


    /* First call reads from an offset into one of the output buffers */
    err |= clSetKernelArg(_summarizationKernel, 0, sizeof(cl_mem), &sumBufs[0]);
    err |= clSetKernelArg(_summarizationKernel, 1, sizeof(cl_mem), &buf);
    err |= clSetKernelArg(_summarizationKernel, 2, sizeof(cl_uint), &nElements);
    err |= clSetKernelArg(_summarizationKernel, 3, sizeof(cl_uint), &offset);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error setting summarization kernel arguments");
        return err;
    }

    local[0] = _summarizationWorkgroupSize;
    global[0] = mwNextMultiple(local[0], nElements);

    err = clEnqueueNDRangeKernel(ci->queue, _summarizationKernel, 1,
                                 NULL, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error enqueuing summarization kernel");
        return err;
    }

    /* Why is this necessary? It seems to frequently break on the 7970 and nowhere else without it */
    err = clFinish(ci->queue);
    //err = clFlush(ci->queue);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error finishing summarization kernel");
        return err;
    }

    /* Later calls swap between summarization buffers without an offset */
    nElements = (cl_uint) mwDivRoundup(global[0], local[0]);
    offset = 0;
    err |= clSetKernelArg(_summarizationKernel, 3, sizeof(cl_uint), &offset);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error setting summarization kernel offset argument");
        return err;
    }

    while (nElements > 1)
    {
        /* Swap old summarization buffer to the input and shrink the range */
        swapBuffers(sumBufs);

        global[0] = mwNextMultiple(local[0], nElements);

        err |= clSetKernelArg(_summarizationKernel, 0, sizeof(cl_mem), &sumBufs[0]);
        err |= clSetKernelArg(_summarizationKernel, 1, sizeof(cl_mem), &sumBufs[1]);
        err |= clSetKernelArg(_summarizationKernel, 2, sizeof(cl_uint), &nElements);
        if (err != CL_SUCCESS)
        {
            mwPerrorCL(err, "Error setting summarization kernel arguments");
            return err;
        }

        /*
        err = clEnqueueBarrier(ci->queue);
        if (err != CL_SUCCESS)
        {
            mwPerrorCL(err, "Error enqueuing summarization barrier");
            return err;
        }
        */

        err = clEnqueueNDRangeKernel(ci->queue, _summarizationKernel, 1,
                                     NULL, global, local,
                                     0, NULL, NULL);
        if (err != CL_SUCCESS)
        {
            mwPerrorCL(err, "Error enqueuing summarization kernel");
            return err;
        }

        err = clFinish(ci->queue);
        if (err != CL_SUCCESS)
        {
            mwPerrorCL(err, "Error finishing summarization kernel");
            return err;
        }

        nElements = mwDivRoundup(global[0], local[0]);
    }


    err = clEnqueueBarrier(ci->queue);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error enqueuing summarization barrier");
        return err;
    }

    err = clEnqueueReadBuffer(ci->queue, sumBufs[0], CL_TRUE,
                              0, 2 * sizeof(real), result,
                              0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error reading summarization result buffer");
        return err;
    }

    resultOut->sum = result[0];
    resultOut->correction = result[1];

    return CL_SUCCESS;
}


static cl_int runIntegralKernel(CLInfo* ci, const RunSizes* runSizes, const size_t offset[1])
{
    cl_int err;
    cl_event ev;

    err = clEnqueueNDRangeKernel(ci->queue,
                                 _separationKernel,
                                 1,
                                 offset, runSizes->global, runSizes->local,
                                 0, NULL, &ev);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error enqueueing integral kernel execution");
        return err;
    }

    /* Give the screen a chance to redraw */
    err = mwCLWaitForEvent(ci, ev, runSizes->initialWait);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to wait for integral event");
        return err;
    }

    return CL_SUCCESS;
}

static cl_int setNuKernelArgs(const IntegralArea* ia, cl_uint nu_step)
{
    cl_int err;
    NuId nuid;

    nuid = calcNuStep(ia, nu_step);
    err = clSetKernelArg(_separationKernel, 13, sizeof(real), &nuid.id);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error setting nu_id argument for step %u", nu_step);
        return err;
    }

    err = clSetKernelArg(_separationKernel, 14, sizeof(cl_uint), &nu_step);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error setting nu_id argument for step %u", nu_step);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int readKernelResults(CLInfo* ci, SeparationCLMem* cm, EvaluationState* es, const IntegralArea* ia)
{
    cl_int err = CL_SUCCESS;
    cl_int i;

    err = runSummarization(ci, cm, ia, 0, &es->bgSumCheckpoint);
    for (i = 0; err == CL_SUCCESS && i < es->numberStreams; ++i)
    {
        err = runSummarization(ci, cm, ia, i + 1, &es->streamSumsCheckpoint[i]);
    }

    return err;
}

static cl_int runNuStep(CLInfo* ci, const IntegralArea* ia, const RunSizes* runSizes, cl_uint nu_step)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;
    size_t offset[1];

    err = setNuKernelArgs(ia, nu_step);
    if (err != CL_SUCCESS)
    {
        mw_printf("Failed to set nu kernel argument\n");
        return err;
    }

    offset[0] = 0;
    for (i = 0; i < runSizes->nChunk && err == CL_SUCCESS; ++i)
    {
        mw_begin_critical_section();
        err = runIntegralKernel(ci, runSizes, offset);
        mw_end_critical_section();

        offset[0] += runSizes->global[0];
    }

    return err;
}

static inline void reportProgress(const AstronomyParameters* ap,
                                  const IntegralArea* ia,
                                  EvaluationState* es,
                                  cl_uint step,
                                  double dt)
{
    cl_long prog;
    prog = es->current_calc_probs + (cl_ulong) ia->mu_steps * ia->r_steps * step;
    mw_fraction_done((cl_double) prog / ap->total_calc_probs);

    if (!BOINC_APPLICATION)
    {
        printf("Step %u: %fms\n", step, dt);
    }
}

static cl_int checkpointCL(CLInfo* ci, SeparationCLMem* cm, const IntegralArea* ia, EvaluationState* es)
{
    cl_int err;

    err = readKernelResults(ci, cm, es, ia);
    if (err != CL_SUCCESS)
        return err;

    err = writeCheckpoint(es) ? MW_CL_ERROR : CL_SUCCESS;
    mw_checkpoint_completed();

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
            mwPerrorCL(err, "Failed to run nu step");
            return err;
        }
        t2 = mwGetTimeMilli();

        dt = t2 - t1;
        tAcc += dt;

        reportProgress(ap, ia, es, es->nu_step + 1, dt);
    }

    es->nu_step = 0;

    mw_printf("Integration time: %f s. Average time per iteration = %f ms\n",
              tAcc / 1000.0, tAcc / (double) ia->nu_steps);

    if (err == CL_SUCCESS)
    {
        err = readKernelResults(ci, cm, es, ia);
        if (err != CL_SUCCESS)
            mw_printf("Failed to read final kernel results\n");

        /* Add final episode to running totals */
        addTmpCheckpointSums(es);
    }

    return err;
}

cl_int integrateCL(const AstronomyParameters* ap,
                   const IntegralArea* ia,
                   const StreamConstants* sc,
                   const StreamGauss sg,
                   EvaluationState* es,
                   const CLRequest* clr,
                   CLInfo* ci)
{
    cl_int err;
    RunSizes runSizes;
    SeparationSizes sizes;
    SeparationCLMem cm = EMPTY_SEPARATION_CL_MEM;

    /* Need to test sizes for each integral, since the area size can change */
    calculateSizes(&sizes, ap, ia);

    if (findRunSizes(&runSizes, ci, &ci->di, ap, ia, clr))
    {
        mw_printf("Failed to find good run sizes\n");
        return MW_CL_ERROR;
    }

    err = createSeparationBuffers(ci, &cm, ap, ia, sc, sg, &sizes);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to create CL buffers");
        return err;
    }

    err = separationSetKernelArgs(&cm, &runSizes);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to set integral kernel arguments");
        return err;
    }

    err = runIntegral(ci, &cm, &runSizes, es, clr, ap, ia);

    releaseSeparationBuffers(&cm);

    separationIntegralGetSums(es);

    return err;
}

