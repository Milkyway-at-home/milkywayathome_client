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

#define BUFSIZE 4096

static const char* cl_integrate_src = "I am not actually a kernel.";

inline static void releaseSeparationCLMem(SeparationCLMem* cm)
{
    clReleaseMemObject(cm->outNu);
    clReleaseMemObject(cm->ap);
    clReleaseMemObject(cm->sc);
    clReleaseMemObject(cm->rConsts);
    clReleaseMemObject(cm->rPoints);
    clReleaseMemObject(cm->nuConsts);
    clReleaseMemObject(cm->ia);
}

inline static cl_int separationSetKernelArgs(const ASTRONOMY_PARAMETERS* ap,
                                             const CLInfo* ci,
                                             SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;

    /* The constant buffers */
    err |= clSetKernelArg(ci->kern, 0, sizeof(cl_mem), &cm->outNu);
    err |= clSetKernelArg(ci->kern, 1, sizeof(cl_mem), &cm->ap);
    err |= clSetKernelArg(ci->kern, 2, sizeof(cl_mem), &cm->sc);
    err |= clSetKernelArg(ci->kern, 3, sizeof(cl_mem), &cm->rConsts);
    err |= clSetKernelArg(ci->kern, 4, sizeof(cl_mem), &cm->rPoints);
    err |= clSetKernelArg(ci->kern, 5, sizeof(cl_mem), &cm->nuConsts);
    err |= clSetKernelArg(ci->kern, 6, sizeof(cl_mem), &cm->ia);

    /* Local workspace used by nu_sum and r_sum */
    err |= clSetKernelArg(ci->kern, 7, sizeof(ST_PROBS) * ap->number_streams, NULL); /* st_probs */
    err |= clSetKernelArg(ci->kern, 8, sizeof(vector) * ap->convolve, NULL);         /* xyz */

    if (err != CL_SUCCESS)
    {
        warn("Error setting kernel arguments: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

const char* src = "\n" \
"__kernel void sampleKernel(                                            \n" \
"   __global float* input,                                              \n" \
"   __global float* output,                                             \n" \
"   const unsigned int count)                                           \n" \
"{                                                                      \n" \
"   int i = get_global_id(0);                                           \n" \
"   if(i < count)                                                       \n" \
"       output[i] = input[i] * input[i];                                \n" \
"}                                                                      \n" \
"\n";

int setupSeparationCL(const ASTRONOMY_PARAMETERS* ap,
                      const STREAM_CONSTANTS* sc,
                      const R_CONSTANTS* r_consts,
                      const R_POINTS* r_points,
                      const NU_CONSTANTS* nu_consts,
                      const INTEGRAL_AREA* ia)
{
    CLInfo ci;
    SeparationCLMem cm;
    char* compileDefs;

    asprintf(&compileDefs,
             "-I/Users/matt/src/milkywayathome_client/separation/include "
             "-DDOUBLEPREC=0 ");

    if (getCLInfo(&ci, CL_DEVICE_TYPE_CPU, "sampleKernel", &src, compileDefs))
        fail("Failed to setup OpenCL device\n");

    free(compileDefs);

    printf("arstarstarst\n");
    mw_finish(EXIT_SUCCESS);


    if (createSeparationBuffers(ap, sc, r_consts, r_points, nu_consts, ia, &ci, &cm) != CL_SUCCESS)
        fail("Failed to create CL buffers\n");

    return 0;
}

