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

#define BUFSIZE 4096

inline static void releaseSeparationCLMem(SeparationCLMem* cm)
{
    clReleaseMemObject(cm->outNu);
    clReleaseMemObject(cm->sc);
    clReleaseMemObject(cm->nuConsts);
}

inline static cl_int separationSetKernelArgs(const ASTRONOMY_PARAMETERS* ap,
                                             const INTEGRAL_AREA* ia,
                                             const CLInfo* ci,
                                             SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;

    /* Output buffer */
    err |= clSetKernelArg(ci->kern, 0, sizeof(cl_mem), &cm->outNu);

    /* The constant arguments */
    err |= clSetKernelArg(ci->kern, 1, sizeof(ASTRONOMY_PARAMETERS), ap);
    err |= clSetKernelArg(ci->kern, 2, sizeof(INTEGRAL_AREA), ia);
    err |= clSetKernelArg(ci->kern, 3, sizeof(cl_mem), &cm->sc);
    err |= clSetKernelArg(ci->kern, 4, sizeof(cl_mem), &cm->nuConsts);


    /* Local workspaces */
    err |= clSetKernelArg(ci->kern, 5, sizeof(ST_PROBS) * ap->number_streams, NULL); /* st_probs */
    err |= clSetKernelArg(ci->kern, 6, sizeof(vector) * ap->convolve, NULL);         /* xyz */
    err |= clSetKernelArg(ci->kern, 7, sizeof(R_POINTS) * ap->convolve, NULL);       /* r_pts */

    if (err != CL_SUCCESS)
    {
        warn("Error setting kernel arguments: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

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


static cl_int enqueueIntegralKernel(CLInfo* ci, SeparationCLMem* cm, const unsigned int r_steps)
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
inline static BG_PROB sumNuResults(BG_PROB* nu_results, const unsigned int r_steps)
{
    unsigned int i;
    BG_PROB bg_prob = ZERO_BG_PROB;

    for (i = 0; i < r_steps; ++i)
        INCADD_BG_PROB(bg_prob, nu_results[i]);

    return bg_prob;
}

int setupSeparationCL(const ASTRONOMY_PARAMETERS* ap,
                      const INTEGRAL_AREA* ia,
                      const STREAM_CONSTANTS* sc,
                      const NU_CONSTANTS* nu_consts)
{
    CLInfo ci;
    SeparationCLMem cm;
    char* compileDefs;
    char* testSrc;

    testSrc = mwReadFile("/Users/matt/src/milkywayathome_client/separation/kernels/test_kernel.cl");

    compileDefs = separationCLDefs(ap,
                                   "-I/Users/matt/src/milkywayathome_client/separation/include "
                                   "-I/Users/matt/src/milkywayathome_client/milkyway/include "
                                   "-DDOUBLEPREC=1 ");

    if (getCLInfo(&ci, CL_DEVICE_TYPE_CPU, "testKernel", &testSrc, compileDefs))
        fail("Failed to setup OpenCL device\n");

    free(testSrc);
    free(compileDefs);

    if (createSeparationBuffers(ap, ia, sc, nu_consts, &ci, &cm) != CL_SUCCESS)
        fail("Failed to create CL buffers\n");

    separationSetKernelArgs(ap, ia, &ci, &cm);
    enqueueIntegralKernel(&ci, &cm, ia->r_steps);

    BG_PROB* nu_results = mallocSafe(sizeof(BG_PROB) * ia->r_steps);

    readIntegralResults(&ci, &cm, nu_results, ia->r_steps);

    printf("arstarstarst\n");

    BG_PROB result = sumNuResults(nu_results, ia->r_steps);
    printf("Result = %g, %g\n", result.bg_int, result.correction);


    mw_finish(EXIT_SUCCESS);

    return 0;
}

