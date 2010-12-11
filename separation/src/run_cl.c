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
#include "integrals_common.h"


static void sumStreamResults(Kahan* probs_results,
                             const real* probs_V_reff_xr_rp3,
                             const cl_uint number_streams)
{
    cl_uint i;

    for (i = 0; i < number_streams; ++i)
        KAHAN_ADD(probs_results[i], probs_V_reff_xr_rp3[i]);
}

static void sumProbsResults(real* probs_results,
                            const real* st_probs_V_reff_xr_rp3_mu_r,
                            const cl_uint mu_steps,
                            const cl_uint r_steps,
                            const cl_uint number_streams)
{
    cl_uint i, j, idx;
    Kahan* probs_sum;

    probs_sum = (Kahan*) mwCalloc(number_streams, sizeof(Kahan));

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            idx = (i * r_steps * number_streams) + (j * number_streams);
            sumStreamResults(probs_sum, &st_probs_V_reff_xr_rp3_mu_r[idx], number_streams);
        }
    }

    for (i = 0; i < number_streams; ++i)
        probs_results[i] = probs_sum[i].sum + probs_sum[i].correction;

    free(probs_sum);
}

static cl_int enqueueIntegralKernel(CLInfo* ci,
                                    const size_t offset[],
                                    const size_t global[],
                                    const size_t local[])
{
    cl_int err;

    err = clEnqueueNDRangeKernel(ci->queue,
                                 ci->kern,
                                 2,
                                 offset, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Error enqueueing integral kernel execution: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

/* The mu and r steps for a nu step were done in parallel. After that
 * we need to add the result to the running total + correction from
 * all the nu steps */
static real sumMuResults(const real* mu_results,
                         const cl_uint mu_steps,
                         const cl_uint r_steps)
{
    cl_uint i, j;
    Kahan bg_prob = ZERO_KAHAN;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            KAHAN_ADD(bg_prob, mu_results[r_steps * i + j]);
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
        warn("Error setting nu_id argument for step %u: %s\n", nu_step, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline real readKernelResults(CLInfo* ci,
                                     SeparationCLMem* cm,
                                     real* probs_results,
                                     const cl_uint mu_steps,
                                     const cl_uint r_steps,
                                     const cl_uint number_streams)
{
    cl_int err;
    real* mu_results;
    real* probs_tmp;
    real result;
    size_t resultSize, probsSize;

    resultSize = sizeof(real) * mu_steps * r_steps;
    mu_results = mapIntegralResults(ci, cm, resultSize);
    if (!mu_results)
    {
        warn("Failed to map integral results\n");
        return NAN;
    }

    result = sumMuResults(mu_results, mu_steps, r_steps);

    err = clEnqueueUnmapMemObject(ci->queue, cm->outMu, mu_results, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to unmap results buffer: %s\n", showCLInt(err));
        return NAN;
    }

    probsSize = sizeof(real) * mu_steps * r_steps * number_streams;
    probs_tmp = mapProbsResults(ci, cm, probsSize);
    if (!probs_tmp)
    {
        warn("Failed to map probs results\n");
        return NAN;
    }

    sumProbsResults(probs_results, probs_tmp, mu_steps, r_steps, number_streams);

    err = clEnqueueUnmapMemObject(ci->queue, cm->outProbs, probs_tmp, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to unmap probs buffer: %s\n", showCLInt(err));
        return NAN;
    }

    return result;
}

static cl_int runNuStep(CLInfo* ci,
                        SeparationCLMem* cm,

                        const IntegralArea* ia,
                        const RunSizes* runSizes,
                        const cl_uint number_streams,
                        const cl_uint nu_step)
{
    cl_int err;
    size_t offset[2];
    unsigned int i;

    err = setNuKernelArgs(ci, ia, nu_step);
    if (err != CL_SUCCESS)
    {
        warn("Failed to set nu kernel argument\n");
        return err;
    }

    offset[1] = nu_step;
    for (i = 0; i < runSizes->numChunks; ++i)
    {
        offset[0] = i * runSizes->chunkSize;

        err = enqueueIntegralKernel(ci, offset, runSizes->global, runSizes->local);
        if (err != CL_SUCCESS)
        {
            warn("Failed to enqueue integral kernel: %s\n", showCLInt(err));
            return err;
        }

        /* Give the screen a chance to redraw */
        err = clFinish(ci->queue);
        if (err != CL_SUCCESS)
        {
            warn("Failed to finish: %s\n", showCLInt(err));
            return err;
        }
    }

    return CL_SUCCESS;
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

static real runIntegral(CLInfo* ci,
                        SeparationCLMem* cm,
                        RunSizes* runSizes,
                        real* probs_results,
                        EvaluationState* es,
                        const CLRequest* clr,
                        const AstronomyParameters* ap,
                        const IntegralArea* ia)
{
    cl_int err;
    real result;
    double t1, t2, dt;
    double tAcc = 0.0;

    for (es->nu_step = 0; es->nu_step < ia->nu_steps; es->nu_step++)
    {
        t1 = mwGetTimeMilli();
        err = runNuStep(ci, cm, ia, runSizes, ap->number_streams, es->nu_step);
        if (err != CL_SUCCESS)
        {
            warn("Failed to run nu step: %s\n", showCLInt(err));
            return NAN;
        }
        t2 = mwGetTimeMilli();

        dt = t2 - t1;
        tAcc += dt;

        reportProgress(ap, ia, es, es->nu_step + 1, dt);
    }

    warn("Integration time: %f s. Average time per iteration = %f ms\n",
         tAcc / 1000.0, tAcc / (double) ia->nu_steps);

    /* Read results from final step */
    result = readKernelResults(ci, cm, probs_results, ia->mu_steps, ia->r_steps, ap->number_streams);
    if (isnan(result))
        warn("Failed to read final kernel results\n");

    return result;
}

real integrateCL(const AstronomyParameters* ap,
                 const IntegralArea* ia,
                 const StreamConstants* sc,
                 const StreamGauss sg,
                 real* probs_results,
                 EvaluationState* es,
                 const CLRequest* clr,
                 CLInfo* ci,
                 DevInfo* di,
                 cl_bool useImages)
{
    real result;
    cl_int err;
    SeparationSizes sizes;
    RunSizes runSizes;
    SeparationCLMem cm = EMPTY_SEPARATION_CL_MEM;

    /* Need to test sizes for each integral, since the area size can change */
    calculateSizes(&sizes, ap, ia);
    if (findGoodRunSizes(&runSizes, ci, di, ia, clr))
    {
        warn("Failed to find good run sizes\n");
        return NAN;
    }

    if (!separationCheckDevCapabilities(di, &sizes))
    {
        warn("Device failed capability check\n");
        return NAN;
    }

    err = createSeparationBuffers(ci, &cm, ap, ia, sc, sg, &sizes, useImages);
    if (err != CL_SUCCESS)
    {
        warn("Failed to create CL buffers: %s\n", showCLInt(err));
        return NAN;
    }

    err = separationSetKernelArgs(ci, &cm, &runSizes);
    if (err != CL_SUCCESS)
    {
        warn("Failed to set integral kernel arguments: %s\n", showCLInt(err));
        return NAN;
    }

    result = runIntegral(ci, &cm, &runSizes, probs_results, es, clr, ap, ia);

    releaseSeparationBuffers(&cm);

    return result;
}

