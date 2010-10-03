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

typedef struct
{
    cl_event ndRange;  /* kernel execution event */
    cl_event outRead;  /* reading main output buffer */
    cl_event probRead; /* reading stream probs buffer */

    cl_event probCalc;  /* custom events for summing temporary output buffers */
    cl_event stProbCalc;
} _SeparationCLEvents;

#define N_SEP_CL_EVENTS (sizeof(_SeparationCLEvents) / sizeof(cl_event))

typedef union
{
    _SeparationCLEvents events;
    cl_event eventv[N_SEP_CL_EVENTS];
} SeparationCLEvents;

static void printSeparationEventTimes(SeparationCLEvents* evs)
{
    printf("NDrange:    %.15g s\n"
           "Read out:   %.15g s\n"
           "Read probs: %.15g s\n",
           mwEventTime(evs->events.ndRange),
           mwEventTime(evs->events.outRead),
           mwEventTime(evs->events.probRead));
}

static real* mapIntegralResults(CLInfo* ci,
                                SeparationCLMem* cm,
                                SeparationCLEvents* evs,
                                size_t resultsSize)
{
    cl_int err;
    real* mapOutMu;

    mapOutMu = (real*) clEnqueueMapBuffer(ci->queue,
                                          cm->outMu,
                                          CL_TRUE, CL_MAP_READ,
                                          0, resultsSize,
                                          1, &evs->events.ndRange, /* Wait for calculation */
                                          &evs->events.outRead,
                                          &err);
    if (err != CL_SUCCESS)
        warn("Error mapping integral result buffer: %s\n", showCLInt(err));

    return mapOutMu;
}

static real* mapProbsResults(CLInfo* ci,
                             SeparationCLMem* cm,
                             SeparationCLEvents* evs,
                             size_t probsResultsSize)
{
    cl_int err;
    real* mapOutProbs;

    mapOutProbs = (real*) clEnqueueMapBuffer(ci->queue,
                                             cm->outProbs,
                                             CL_TRUE, CL_MAP_READ,
                                             0, probsResultsSize,
                                             1, &evs->events.ndRange, /* Wait for calculation */
                                             &evs->events.probRead,
                                             &err);
    if (err != CL_SUCCESS)
        warn("Error mapping probs result buffer: %s\n", showCLInt(err));

    return mapOutProbs;
}

static inline void sumStreamResults(ST_PROBS* probs_results,
                                    const real* probs_V_reff_xr_rp3,
                                    const cl_uint number_streams)
{
    cl_uint i;

    for (i = 0; i < number_streams; ++i)
        KAHAN_ADD(probs_results[i].st_prob_int, probs_V_reff_xr_rp3[i], probs_results[i].st_prob_int_c);
}

static inline void sumProbsResults(ST_PROBS* probs_results,
                                   const real* st_probs_V_reff_xr_rp3_mu_r,
                                   const cl_uint mu_steps,
                                   const cl_uint r_steps,
                                   const cl_uint number_streams)
{
    cl_uint i, j, idx;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            idx = (i * r_steps * number_streams) + (j * number_streams);
            sumStreamResults(probs_results, &st_probs_V_reff_xr_rp3_mu_r[idx], number_streams);
        }
    }
}

static cl_int enqueueIntegralKernel(CLInfo* ci,
                                    SeparationCLEvents* evs,
                                    const cl_uint mu_steps,
                                    const cl_uint r_steps)
{
    cl_int err;
    const size_t global[] = { mu_steps, r_steps };

    err = clEnqueueNDRangeKernel(ci->queue,
                                 ci->kern,
                                 2,
                                 NULL, global, NULL,
                                 0, NULL, &evs->events.ndRange);
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
                                const cl_uint mu_steps,
                                const cl_uint r_steps)
{
    cl_uint i, j;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            KAHAN_ADD(bg_prob->bg_int, mu_results[r_steps * i + j], bg_prob->correction);
        }
    }
}

static cl_int setNuKernelArg(CLInfo* ci, const cl_uint nu_step)
{
    cl_int err;

    err = clSetKernelArg(ci->kern, 7, sizeof(cl_uint), &nu_step);
    if (err != CL_SUCCESS)
    {
        warn("Error setting nu step argument for step %u: %s\n", nu_step, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int readKernelResults(CLInfo* ci,
                                       SeparationCLMem* cm,
                                       SeparationCLEvents* evs,
                                       BG_PROB* bg_progress,
                                       ST_PROBS* probs_results,
                                       const cl_uint mu_steps,
                                       const cl_uint r_steps,
                                       const cl_uint number_streams)
{
    cl_int err;
    real* mu_results;
    real* probs_tmp;

    /* We must map/unmap the buffer on each step. A kernel writing to
     * a mapped buffer is undefined. */

    size_t resultSize = sizeof(real) * mu_steps * r_steps;
    mu_results = mapIntegralResults(ci, cm, evs, resultSize);
    if (!mu_results)
    {
        warn("Failed to map integral results\n");
        return -1;
    }

    size_t probsSize = sizeof(real) * mu_steps * r_steps * number_streams;
    probs_tmp = mapProbsResults(ci, cm, evs, probsSize);
    if (!probs_tmp)
    {
        warn("Failed to map probs results\n");
        return -1;
    }

    //cl_event sumEvent = mwCreateEvent(ci);

    sumMuResults(bg_progress, mu_results, mu_steps, r_steps);
    sumProbsResults(probs_results, probs_tmp, mu_steps, r_steps, number_streams);

    //mwFinishEvent(sumEvent);


    err = clEnqueueUnmapMemObject(ci->queue, cm->outMu, mu_results, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to unmap results buffer: %s\n", showCLInt(err));
        return err;
    }

    err = clEnqueueUnmapMemObject(ci->queue, cm->outProbs, probs_tmp, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to unmap probs buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static cl_int runNuStep(CLInfo* ci,
                        SeparationCLMem* cm,
                        SeparationCLEvents* evs,

                        BG_PROB* bg_progress,    /* Accumulating results over nu steps */
                        ST_PROBS* probs_results,

                        const cl_uint mu_steps,
                        const cl_uint r_steps,
                        const cl_uint number_streams,
                        const cl_uint nu_step)
{
    cl_int err;

    err = setNuKernelArg(ci, nu_step);
    if (err != CL_SUCCESS)
    {
        warn("Failed to set nu kernel argument\n");
        return err;
    }

    err = separationSwapOutputBuffers(ci, cm);
    if (err != CL_SUCCESS)
    {
        warn("Failed to set output buffer arguments on step %u: %s\n", nu_step, showCLInt(err));
        return err;
    }

    err = enqueueIntegralKernel(ci, evs, mu_steps, r_steps);
    if (err != CL_SUCCESS)
    {
        warn("Failed to enqueue integral kernel: %s\n", showCLInt(err));
        return err;
    }

    err = readKernelResults(ci, cm, evs, bg_progress, probs_results, mu_steps, r_steps, number_streams);
    if (err != CL_SUCCESS)
    {
        warn("Failed to read kernel results: %s\n", showCLInt(err));
        return err;
    }

    printf("Nu step %u:\n", nu_step);
    //printSeparationEventTimes(evs);

    return CL_SUCCESS;
}

static real runIntegral(CLInfo* ci,
                        SeparationCLMem* cm,
                        ST_PROBS* probs_results,
                        const ASTRONOMY_PARAMETERS* ap,
                        const INTEGRAL_AREA* ia)
{
    cl_uint i;
    cl_int err;
    BG_PROB bg_sum = ZERO_BG_PROB;
    SeparationCLEvents evs;

    //mwEnableProfiling(ci);

    for (i = 0; i < ia->nu_steps; ++i)
    {
        err = runNuStep(ci, cm, &evs,
                        &bg_sum, probs_results,
                        ia->mu_steps, ia->r_steps, ap->number_streams, i);

        if (err != CL_SUCCESS)
        {
            warn("Failed to run nu step: %s\n", showCLInt(err));
            bg_sum.bg_int = bg_sum.correction = NAN;
            break;
        }
    }

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

