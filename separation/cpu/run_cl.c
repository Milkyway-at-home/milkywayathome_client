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
#include <assert.h>
#include <OpenCL/cl.h>
#include <OpenCL/cl_ext.h>

#include "milkyway_util.h"
#include "show_cl_types.h"
#include "separation.h"
#include "setup_cl.h"
#include "separation_cl_buffers.h"
#include "separation_cl_defs.h"
#include "run_cl.h"

static cl_int readIntegralResults(CLInfo* ci,
                                  SeparationCLMem* cm,
                                  BG_PROB* nu_results,
                                  const unsigned int r_steps)
{
    cl_int err;
    err = clEnqueueReadBuffer(ci->queue,
                              cm->outNu,
                              CL_TRUE,
                              0, sizeof(BG_PROB) * r_steps, nu_results,
                              0, NULL, NULL);

    if (err != CL_SUCCESS)
    {
        warn("Error reading integral result buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}


static cl_int enqueueIntegralKernel(CLInfo* ci, const unsigned int r_steps)
{
    cl_int err;
    const size_t global[] = { r_steps };

    err = clEnqueueNDRangeKernel(ci->queue,
                                 ci->kern,
                                 1,
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
inline static double sumNuResults(BG_PROB* nu_results, const unsigned int r_steps)
{
    unsigned int i;
    BG_PROB bg_prob = ZERO_BG_PROB;

    for (i = 0; i < r_steps; ++i)
        INCADD_BG_PROB(bg_prob, nu_results[i]);

    return bg_prob.bg_int + bg_prob.correction;
}


static double runIntegral(CLInfo* ci,
                          SeparationCLMem* cm,
                          const unsigned int r_steps)
{
    BG_PROB* nu_results;
    double bg_result;

    nu_results = mallocSafe(sizeof(BG_PROB) * r_steps);

    enqueueIntegralKernel(ci, r_steps);
    readIntegralResults(ci, cm, nu_results, r_steps);
    bg_result = sumNuResults(nu_results, r_steps);

    free(nu_results);

    return bg_result;
}

void separationCL(const ASTRONOMY_PARAMETERS* ap,
                  const INTEGRAL_AREA* ia,
                  const STREAM_CONSTANTS* sc)
{
    CLInfo ci;
    SeparationCLMem cm;
    NU_CONSTANTS* nu_consts;
    double result;

    nu_consts = prepare_nu_constants(ia->nu_steps, ia->nu_step_size, ia->nu_min);
    setupSeparationCL(ap, ia, sc, nu_consts, &ci, &cm);
    free(nu_consts);

    result = runIntegral(&ci, &cm, ia->r_steps);
    releaseSeparationBuffers(&cm);
    destroyCLInfo(&ci);


    printf("Result = %g\n", result);

    mw_finish(EXIT_SUCCESS);
}


