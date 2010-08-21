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
#include "build_cl.h"

#define BUFSIZE 4096

cl_int printCLExtensions(cl_device_id dev)
{
    cl_int err;
    char exts[1024];

    err = clGetDeviceInfo(dev, CL_DEVICE_EXTENSIONS, sizeof(exts), exts, NULL);

    if (err != CL_SUCCESS)
    {
        warn("Failed to print CL extensions for device: %s", exts);
        return err;
    }
    else
        printf("Extensions: %s\n", exts);

    return CL_SUCCESS;
}

cl_int destroyCLInfo(CLInfo* ci)
{
    cl_int err = CL_SUCCESS;
    err |= clReleaseCommandQueue(ci->queue);
    err |= clReleaseProgram(ci->prog);
    err |= clReleaseKernel(ci->kern);
    err |= clReleaseContext(ci->clctx);

    /* TODO: or'ing the err and showing = useless */
    if (err)
        warn("Error cleaning up CLInfo: %s\n", showCLInt(err));

    return err;
}

void milkywayBuildCB(cl_program prog, void* user_data)
{
    char buildLog[BUFSIZE] = "";
    cl_int infoErr;
    cl_build_status stat;
    size_t failSize;

    CLInfo* ci = (CLInfo*) user_data;

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
}

static cl_int milkywayBuildProgram(CLInfo* ci, const char** src, cl_uint srcCount, const char* compileDefs)
{
    cl_int err = CL_SUCCESS;

    ci->prog = clCreateProgramWithSource(ci->clctx, 2, src, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating program: %s\n", showCLInt(err));
        return err;
    }

    err = clBuildProgram(ci->prog, 1, &ci->dev, compileDefs, NULL, NULL);
    if (err != CL_SUCCESS)
        warn("Build failure: %s\n", showCLInt(err));

    return err;
}

static cl_int createCtxQueue(CLInfo* ci)
{
    cl_int err = CL_SUCCESS;

    ci->clctx = clCreateContext(NULL, 1, &ci->dev, clLogMessagesToStdoutAPPLE, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating context: %s\n", showCLInt(err));
        return err;
    }

    ci->queue = clCreateCommandQueue(ci->clctx, ci->dev, 0, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating command queue: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static cl_int getDevInfo(CLInfo* ci, cl_device_type type)
{
    cl_int err;
    cl_uint maxComputeUnits, clockFreq;
    cl_ulong memSize;

    err = clGetDeviceIDs(NULL, type, 1, &ci->dev, &ci->devCount);
    if (err != CL_SUCCESS)
    {
        warn("Error getting device: %s\n", showCLInt(err));
        return err;
    }

    if (ci->devCount == 0)
    {
        warn("Didn't find any %s devices\n", showCLDeviceType(type));
        return -1; /* FIXME: Meaningful error? */
    }

    ci->devType = type;

    /* Print some device information */
    clGetDeviceInfo(ci->dev, CL_DEVICE_MAX_COMPUTE_UNITS,   sizeof(cl_uint),  &maxComputeUnits, NULL);
    clGetDeviceInfo(ci->dev, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint),  &clockFreq, NULL);
    clGetDeviceInfo(ci->dev, CL_DEVICE_GLOBAL_MEM_SIZE,     sizeof(cl_ulong), &memSize, NULL);

    printf("arst device %s: %u %u %lu\n",
           showCLDeviceType(type), maxComputeUnits, clockFreq, (unsigned long) memSize);
    printCLExtensions(ci->dev);

    return CL_SUCCESS;
}

/* Query one device specified by type, create a context, command
 * queue, and compile the kernel for it
 * Returns: non-zero on failure
 *  TODO: Multiple device support
 *  TODO: Caching of compiled binaries
 */
cl_int getCLInfo(CLInfo* ci,
                 cl_device_type type,
                 const char* kernName,
                 const char** src,
                 const cl_uint srcCount,
                 const char* compileDefs)
{
    cl_int err;

    err = getDevInfo(ci, type);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get information about device\n");
        return err;
    }

    err = createCtxQueue(ci);
    if (err != CL_SUCCESS)
    {
        warn("Error creating CL context and command queue: %s\n", showCLInt(err));
        return err;
    }

    err = milkywayBuildProgram(ci, src, srcCount, compileDefs);
    if (err != CL_SUCCESS)
    {
        warn("Error building program: %s\n", showCLInt(err));
        return err;
    }

    ci->kern = clCreateKernel(ci->prog, kernName, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating kernel '%s': %s\n", kernName, showCLInt(err));
        return err;
    }

    clUnloadCompiler();

    return CL_SUCCESS;
}

