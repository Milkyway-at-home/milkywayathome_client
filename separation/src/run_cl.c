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
#include "show_cl_types.h"
#include "setup_cl.h"
#include "milkyway_cl.h"
#include "separation_cl_buffers.h"
#include "separation_cl_defs.h"
#include "calculated_constants.h"
#include "run_cl.h"
#include "r_points.h"


static cl_int readIntegralResults(CLInfo* ci,
                                  SeparationCLMem* cm,
                                  real* mu_results,
                                  size_t resultsSize)
{
    cl_int err;

    err = clEnqueueReadBuffer(ci->queue,
                              cm->outMu,
                              CL_TRUE,
                              0, resultsSize, mu_results,
                              0, NULL, NULL);

    if (err != CL_SUCCESS)
    {
        warn("Error reading integral result buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline void sumStreamResults(ST_PROBS* probs_results,
                                    const real* probs_V_reff_xr_rp3,
                                    const unsigned int number_streams)
{
    unsigned int i;

    for (i = 0; i < number_streams; ++i)
        KAHAN_ADD(probs_results[i].st_prob_int, probs_V_reff_xr_rp3[i], probs_results[i].st_prob_int_c);
}

static inline void sumProbsResults(ST_PROBS* probs_results,
                                   const real* st_probs_V_reff_xr_rp3_mu_r,
                                   const unsigned int mu_steps,
                                   const unsigned int r_steps,
                                   const unsigned int number_streams)
{
    unsigned int i, j, idx;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            idx = (i * r_steps * number_streams) + (j * number_streams);
            sumStreamResults(probs_results, &st_probs_V_reff_xr_rp3_mu_r[idx], number_streams);
        }
    }
}

static cl_int readProbsResults(CLInfo* ci,
                               SeparationCLMem* cm,
                               real* probs_tmp,
                               size_t size)
{
    cl_int err;

    err = clEnqueueReadBuffer(ci->queue,
                              cm->outProbs,
                              CL_TRUE,
                              0, size, probs_tmp,
                              0, NULL, NULL);

    if (err != CL_SUCCESS)
        warn("Error reading probs result buffer for stream: %s\n", showCLInt(err));

    return err;
}

static cl_int enqueueIntegralKernel(CLInfo* ci,
                                    const unsigned int mu_steps,
                                    const unsigned int r_steps)
{
    cl_int err;
    const size_t global[] = { mu_steps, r_steps };

    err = clEnqueueNDRangeKernel(ci->queue,
                                 ci->kern,
                                 2,
                                 NULL, global, NULL,
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
static inline void sumMuResults(BG_PROB* bg_prob,
                                const real* mu_results,
                                const unsigned int mu_steps,
                                const unsigned int r_steps)
{
    unsigned int i, j;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            KAHAN_ADD(bg_prob->bg_int, mu_results[r_steps * i + j], bg_prob->correction);
        }
    }
}

static cl_int setNuKernelArg(CLInfo* ci, SeparationCLMem* cm, const unsigned int nu_step)
{
    cl_int err;

    err = clSetKernelArg(ci->kern, 6, sizeof(unsigned int), &nu_step);
    if (err != CL_SUCCESS)
    {
        warn("Error setting nu step argument for step %u: %s\n", nu_step, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static cl_int runNuStep(CLInfo* ci,
                        SeparationCLMem* cm,

                        BG_PROB* bg_progress,    /* Accumulating results over nu steps */
                        ST_PROBS* probs_results,

                        real* mu_results,        /* Scratch buffer for reading from cl buffer */
                        real* probs_tmp,

                        const unsigned int mu_steps,
                        const unsigned int r_steps,
                        const unsigned int number_streams,
                        const unsigned int nu_step)
{
    cl_int err;
    size_t resultSize = sizeof(real) * mu_steps * r_steps;
    size_t probsSize = sizeof(real) * mu_steps * r_steps * number_streams;

    err = setNuKernelArg(ci, cm, nu_step);
    if (err != CL_SUCCESS)
    {
        warn("Failed to set nu kernel argument\n");
        return err;
    }

    err = enqueueIntegralKernel(ci, mu_steps, r_steps);
    if (err != CL_SUCCESS)
    {
        warn("Failed to enqueue integral kernel: %s\n", showCLInt(err));
        return err;
    }

    err = readIntegralResults(ci, cm, mu_results, resultSize);
    if (err != CL_SUCCESS)
    {
        warn("Failed to read integral results: %s\n", showCLInt(err));
        return err;
    }

    err = readProbsResults(ci, cm, probs_tmp, probsSize);
    if (err != CL_SUCCESS)
    {
        warn("Failed to read probs results: %s\n", showCLInt(err));
        return err;
    }

    sumMuResults(bg_progress, mu_results, mu_steps, r_steps);
    sumProbsResults(probs_results, probs_tmp, mu_steps, r_steps, number_streams);

    return CL_SUCCESS;
}

static real runIntegral(CLInfo* ci,
                        SeparationCLMem* cm,
                        ST_PROBS* probs_results,
                        const ASTRONOMY_PARAMETERS* ap,
                        const INTEGRAL_AREA* ia)
{
    unsigned int i;
    real* mu_results;
    real* probs_tmp;
    cl_int err;
    BG_PROB bg_sum = ZERO_BG_PROB;
    size_t resultSize = sizeof(real) * ia->mu_steps * ia->r_steps;
    size_t probsTmpSize = sizeof(real) * ia->mu_steps * ia->r_steps * ap->number_streams;

    mu_results = (real*) mwMallocAligned(resultSize, sizeof(real));
    probs_tmp = (real*) mwMallocAligned(probsTmpSize, sizeof(real));

    for (i = 0; i < ia->nu_steps; ++i)
    {
        err = runNuStep(ci, cm,
                        &bg_sum, probs_results,
                        mu_results, probs_tmp,
                        ia->mu_steps, ia->r_steps, ap->number_streams, i);

        if (err != CL_SUCCESS)
        {
            warn("Failed to run nu step: %s\n", showCLInt(err));
            bg_sum.bg_int = bg_sum.correction = NAN;
            break;
        }
    }

    mwAlignedFree(mu_results);
    mwAlignedFree(probs_tmp);

    return bg_sum.bg_int + bg_sum.correction;
}

/* FIXME: This can only work right now for 1 integral */
real integrateCL(const ASTRONOMY_PARAMETERS* ap,
                 const INTEGRAL_AREA* ia,
                 const STREAM_CONSTANTS* sc,
                 const STREAM_GAUSS* sg,
                 ST_PROBS* probs_results)
{
    real result = NAN;
    CLInfo ci = EMPTY_CL_INFO;
    SeparationCLMem cm = EMPTY_SEPARATION_CL_MEM;

    if (setupSeparationCL(ap, ia, sc, sg, &ci, &cm) != CL_SUCCESS)
        warn("Failed to setup up CL\n");
    else
        result = runIntegral(&ci, &cm, probs_results, ap, ia);

    destroyCLInfo(&ci);
    releaseSeparationBuffers(&cm);

    return result;
}

