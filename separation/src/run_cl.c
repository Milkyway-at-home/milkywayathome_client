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
                                  BG_PROB* mu_results,
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

inline static void sumStreamResults(ST_PROBS* probs_results,
                                    ST_PROBS* probs,
                                    const unsigned int number_streams)
{
    unsigned int i;
    for (i = 0; i < number_streams; ++i)
        KAHAN_ADD(probs_results[i].st_prob_int, probs[i].st_prob_int, probs[i].st_prob_int_c);
}

inline static void sumProbsResults(ST_PROBS* probs_results,
                                   ST_PROBS* probs_r_nu,
                                   const unsigned int r_steps,
                                   const unsigned int nu_steps,
                                   const unsigned int number_streams)
{
    unsigned int i, j, idx;

    for (i = 0; i < r_steps; ++i)
    {
        for (j = 0; j < nu_steps; ++j)
        {
            idx = (i * nu_steps * number_streams) + (j * number_streams);
            sumStreamResults(probs_results, &probs_r_nu[idx], number_streams);
        }
    }
}

static cl_int readProbsResults(CLInfo* ci,
                               SeparationCLMem* cm,
                               ST_PROBS* probs_results,
                               const unsigned int r_steps,
                               const unsigned int nu_steps,
                               const unsigned int number_streams)
{
    ST_PROBS* probs_tmp;
    cl_int err = CL_SUCCESS;

    size_t size = sizeof(ST_PROBS) * r_steps * nu_steps * number_streams;
    probs_tmp = mallocSafe(size);

    err = clEnqueueReadBuffer(ci->queue,
                              cm->outProbs,
                              CL_TRUE,
                              0, size, probs_tmp,
                              0, NULL, NULL);

    if (err != CL_SUCCESS)
        warn("Error reading probs result buffer for stream: %s\n", showCLInt(err));
    else
        sumProbsResults(probs_results, probs_tmp, r_steps, nu_steps, number_streams);

    free(probs_tmp);
    return err;
}

static cl_int enqueueIntegralKernel(CLInfo* ci,
                                    const unsigned int r_steps,
                                    const unsigned int nu_steps)
{
    cl_int err;
    const size_t global[] = { r_steps, nu_steps };
    //const size_t local[] = { 0, nu_steps };
    //printf("local = %u\n", nu_steps);

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

/* The nu steps were done in parallel. After that we need to sum the results */
inline static real sumMuResults(BG_PROB* mu_results,
                                const unsigned int r_steps,
                                const unsigned int nu_steps)
{
    unsigned int i, j;
    BG_PROB bg_prob = ZERO_BG_PROB;

    for (i = 0; i < r_steps; ++i)
    {
        for (j = 0; j < nu_steps; ++j)
            INCADD_BG_PROB(bg_prob, mu_results[nu_steps * i + j]);
    }

    return bg_prob.bg_int + bg_prob.correction;
}

static R_POINTS* prepare_r_pts(const ASTRONOMY_PARAMETERS* ap,
                               const INTEGRAL_AREA* ia,
                               const STREAM_CONSTANTS* sc,
                               const STREAM_GAUSS* sg,
                               const unsigned int r_steps)
{
    unsigned int i;
    R_PRIME rp;
    const unsigned int nconvolve = ap->convolve;
    R_POINTS* r_pts_all = mallocSafe(sizeof(R_POINTS) * r_steps * nconvolve);

    for (i = 0; i < r_steps; ++i)
    {
        rp = calcRPrime(ia, i);
        set_r_points(ap, sg, nconvolve, rp.rPrime, &r_pts_all[i * nconvolve]);
    }

    return r_pts_all;
}

static real runIntegral(CLInfo* ci,
                        SeparationCLMem* cm,
                        ST_PROBS* probs_results,
                        const unsigned int r_steps,
                        const unsigned int nu_steps,
                        const unsigned int number_streams)
{
    BG_PROB* mu_results;
    real bg_result;
    size_t resultSize = sizeof(BG_PROB) * r_steps * nu_steps;

    enqueueIntegralKernel(ci, r_steps, nu_steps);

    mu_results = mallocSafe(resultSize);
    readIntegralResults(ci, cm, mu_results, resultSize);

    printf("Read integral results\n");
    bg_result = sumMuResults(mu_results, r_steps, nu_steps);
    free(mu_results);

    readProbsResults(ci, cm, probs_results, r_steps, nu_steps, number_streams);

    return bg_result;
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
    R_POINTS* r_pts_all;

    r_pts_all = prepare_r_pts(ap, ia, sc, sg, ia->r_steps);
    if (setupSeparationCL(ap, ia, sc, r_pts_all, &ci, &cm) != CL_SUCCESS)
        warn("Failed to setup up CL\n");
    else
        result = runIntegral(&ci, &cm, probs_results, ia->r_steps, ia->nu_steps, ap->number_streams);

    destroyCLInfo(&ci);
    releaseSeparationBuffers(&cm);
    free(r_pts_all);

    return result;
}

