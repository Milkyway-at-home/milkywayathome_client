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

#include "milkyway_util.h"
#include "show_cl_types.h"
#include "milkyway_cl.h"
#include "setup_cl.h"
#include "separation_cl_buffers.h"
#include "separation_cl_defs.h"

#define BUFSIZE 4096

inline static cl_int separationSetKernelArgs(const ASTRONOMY_PARAMETERS* ap,
                                             const INTEGRAL_AREA* ia,
                                             const CLInfo* ci,
                                             SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;

    /* Output buffer */
    err |= clSetKernelArg(ci->kern, 0, sizeof(cl_mem), &cm->outNu);
    err |= clSetKernelArg(ci->kern, 1, sizeof(cl_mem), &cm->outProbs);

    /* The constant arguments */
    err |= clSetKernelArg(ci->kern, 2, sizeof(ASTRONOMY_PARAMETERS), ap);
    err |= clSetKernelArg(ci->kern, 3, sizeof(INTEGRAL_AREA), ia);
    err |= clSetKernelArg(ci->kern, 4, sizeof(cl_mem), &cm->sc);
    err |= clSetKernelArg(ci->kern, 5, sizeof(cl_mem), &cm->sg);
    err |= clSetKernelArg(ci->kern, 6, sizeof(cl_mem), &cm->nuConsts);

    /* Local workspaces */
    err |= clSetKernelArg(ci->kern, 7, sizeof(ST_PROBS) * ap->number_streams, NULL); /* st_probs */
    err |= clSetKernelArg(ci->kern, 8, sizeof(vector) * ap->convolve, NULL);         /* xyz */
    err |= clSetKernelArg(ci->kern, 9, sizeof(R_POINTS) * ap->convolve, NULL);       /* r_pts */

    if (err != CL_SUCCESS)
    {
        warn("Error setting kernel arguments: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

#if DOUBLEPREC
  #define DOUBLEPREC_DEF_STRING "-D DOUBLEPREC=1"
#else
  #define DOUBLEPREC_DEF_STRING "-D DOUBLEPREC=0 -cl-single-precision-constant "
#endif /* DOUBLEPREC */

cl_int setupSeparationCL(const ASTRONOMY_PARAMETERS* ap,
                         const INTEGRAL_AREA* ia,
                         const STREAM_CONSTANTS* sc,
                         const STREAM_GAUSS* sg,
                         const NU_CONSTANTS* nu_consts,
                         CLInfo* ci,
                         SeparationCLMem* cm)
{
    cl_int err;
    //char* compileDefs;
    char* kernelSrc;
    char* rPointsSrc;

    static const char* extraDefs = DOUBLEPREC_DEF_STRING
                                   "-cl-strict-aliasing "
                                   "-cl-finite-math-only "
                                   "-I../src "
                                   "-I../include "
                                   "-I../../milkyway/include ";

    kernelSrc = mwReadFile("../kernels/integrals.cl");
    if (!kernelSrc)
    {
        warn("Failed to read kernel file\n");
        return -1;
    }

    rPointsSrc = mwReadFile("../src/r_points.c");
    if (!rPointsSrc)
    {
        warn("Failed to read r_points file\n");
        return -1;
    }

    char* allSrc[] = { rPointsSrc, kernelSrc };

    //compileDefs = separationCLDefs(ap, extraDefs);
    err = getCLInfo(ci, CL_DEVICE_TYPE_CPU, "r_sum_kernel", allSrc, 2, extraDefs);

    free(kernelSrc);
    free(rPointsSrc);
    //free(compileDefs);

    if (err != CL_SUCCESS)
    {
        fail("Failed to setup OpenCL device: %s\n", showCLInt(err));
        return err;
    }

    err = createSeparationBuffers(ap, ia, sc, sg, nu_consts, ci, cm);
    if (err != CL_SUCCESS)
    {
        fail("Failed to create CL buffers: %s\n", showCLInt(err));
        return err;
    }

    err = separationSetKernelArgs(ap, ia, ci, cm);
    if (err != CL_SUCCESS)
    {
        fail("Failed to set integral kernel arguments: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

