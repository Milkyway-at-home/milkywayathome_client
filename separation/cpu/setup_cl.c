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

#define BUFSIZE 4096

static const char* cl_integrate_src = "I am not actually a kernel.";

inline static void destroySeparationCLInfo(SeparationCLInfo* ci)
{
    cl_int err = CL_SUCCESS;
    err |= clReleaseCommandQueue(ci->queue);
    err |= clReleaseProgram(ci->prog);
    err |= clReleaseKernel(ci->kern);
    err |= clReleaseContext(ci->clctx);

    /* TODO: or'ing the err and showing = useless */
    if (err)
        warn("Error cleaning up SeparationCLInfo: %s\n", showCLInt(err));
}

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
                                             const SeparationCLInfo* ci,
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

static cl_int buildProgram(SeparationCLInfo* ci)
{
    cl_int err;
    char* compileDefinitions;
    char buildLog[BUFSIZE] = "";
    cl_int infoErr;
    cl_build_status stat;
    size_t failSize;

    asprintf(&compileDefinitions,
             "-I/Users/matt/src/milkywayathome_client/separation/include "
             "-DDOUBLEPREC=0 ");

    err = clBuildProgram(ci->prog,
                         1,
                         &ci->dev,
                         compileDefinitions,
                         NULL,
                         NULL);

    free(compileDefinitions);

    if (err == CL_SUCCESS)
        return CL_SUCCESS;

    infoErr = clGetProgramBuildInfo(ci->prog,
                                    ci->dev,
                                    CL_PROGRAM_BUILD_STATUS,
                                    sizeof(stat),
                                    &stat,
                                    NULL);

    if (infoErr != CL_SUCCESS)
        warn("Get build status failed: %s\n", showCLInt(infoErr));
    else
        printf("Build status: %s\n", showCLBuildStatus(stat));

    clGetProgramBuildInfo(ci->prog,
                          ci->dev,
                          CL_PROGRAM_BUILD_LOG,
                          sizeof(buildLog),
                          buildLog,
                          &failSize);

    if (failSize > BUFSIZE)
    {
        char* bigBuf = callocSafe(sizeof(char), failSize + 1);

        clGetProgramBuildInfo(ci->prog,
                              ci->dev,
                              CL_PROGRAM_BUILD_LOG,
                              failSize,
                              bigBuf,
                              NULL);

        printf("Large build message: \n%s\n", bigBuf);
        free(bigBuf);
    }

    warn("Build failure: %s: log = %s\n", showCLInt(err), buildLog);

    return err;
}


/* Query one device specified by type, create a context, command
 * queue, and compile the kernel for it
 * Returns: non-zero on failure
 *  TODO: Multiple device support
 *  TODO: Caching of compiled binaries
 */
static int getCLInfo(SeparationCLInfo* ci, cl_device_type type)
{
    cl_int err;
    cl_uint maxComputeUnits, clockFreq;
    cl_ulong memSize;

    err = clGetDeviceIDs(NULL, type, 1, &ci->dev, &ci->devCount);
    if (err != CL_SUCCESS)
    {
        warn("Error getting device: %s\n", showCLInt(err));
        return 1;
    }

    if (ci->devCount == 0)
    {
        warn("Didn't find any %s devices\n", showCLDeviceType(type));
        return 1;
    }

    ci->devType = type;

    ci->clctx = clCreateContext(NULL, 1, &ci->dev, clLogMessagesToStdoutAPPLE, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating context: %s\n", showCLInt(err));
        return 1;
    }

    ci->queue = clCreateCommandQueue(ci->clctx, ci->dev, 0, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating command Queue: %s\n", showCLInt(err));
        return 1;
    }

    /* Print some device information */
    clGetDeviceInfo(ci->dev, CL_DEVICE_MAX_COMPUTE_UNITS,   sizeof(cl_uint),  &maxComputeUnits, NULL);
    clGetDeviceInfo(ci->dev, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint),  &clockFreq, NULL);
    clGetDeviceInfo(ci->dev, CL_DEVICE_GLOBAL_MEM_SIZE,     sizeof(cl_ulong), &memSize, NULL);

    printf("arst device %s: %u %u %lu\n",
           showCLDeviceType(type), maxComputeUnits, clockFreq, (unsigned long) memSize);

    ci->prog = clCreateProgramWithSource(ci->clctx, 1, &cl_integrate_src, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating program: %s\n", showCLInt(err));
        return 1;
    }

    buildProgram(ci);

    ci->kern = clCreateKernel(ci->prog, "I don't exist yet", &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating kernel: %s\n", showCLInt(err));
        return 1;
    }

    clUnloadCompiler();

    return 0;
}


int setupSeparationCL(const ASTRONOMY_PARAMETERS* ap,
                      const STREAM_CONSTANTS* sc,
                      const R_CONSTANTS* r_consts,
                      const R_POINTS* r_points,
                      const NU_CONSTANTS* nu_st,
                      const INTEGRAL_AREA* ia)
{
    SeparationCLMem cm;
    SeparationCLInfo ci;

    if (getCLInfo(&ci, CL_DEVICE_TYPE_CPU))
        fail("Failed to setup OpenCL device\n");

    if (createSeparationBuffers(ap, sc, r_consts, r_points, nu_st, ia, &ci, &cm) != CL_SUCCESS)
        fail("Failed to create CL buffers\n");

    return 0;
}

