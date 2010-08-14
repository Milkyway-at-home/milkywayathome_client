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

#include "nbody_util.h"
#include "nbody_priv.h"
#include "gravmap_opencl.h"
#include "ckernels/cl_gravmap.h"
#include "show_cl_types.h"
//#include "ckernels/cl_nbody_types.h"

inline static void destroyNBodyCLInfo(NBodyCLInfo* ci)
{
    cl_int err = CL_SUCCESS;
    err |= clReleaseCommandQueue(ci->queue);
    err |= clReleaseProgram(ci->prog);
    err |= clReleaseKernel(ci->kern);
    err |= clReleaseContext(ci->clctx);

    /* TODO: or'ing the err and showing = useless */
    if (err)
        warn("Error cleaning up NBodyCLInfo: %s\n", showCLInt(err));
}

inline static void releaseNBodyCLMem(NBodyCLMem* cm)
{
    clReleaseMemObject(cm->acc);
    clReleaseMemObject(cm->bodies);
    clReleaseMemObject(cm->root);
    clReleaseMemObject(cm->nbctx);
}

inline static int nbodySetKernelArgs(const NBodyCLInfo* ci, NBodyCLMem* cm, const size_t nbody)
{
    cl_int err = CL_SUCCESS;
    err |= clSetKernelArg(ci->kern, 0, sizeof(cl_mem), &cm->nbctx);
    err |= clSetKernelArg(ci->kern, 1, sizeof(cl_mem), &cm->root);
    err |= clSetKernelArg(ci->kern, 2, sizeof(cl_mem), &cm->bodies);
    err |= clSetKernelArg(ci->kern, 3, sizeof(cl_mem), &cm->acc);
    err |= clSetKernelArg(ci->kern, 4, sizeof(size_t), &nbody);
    if (err != CL_SUCCESS)
    {
        warn("Error setting kernel arguments: %s\n", showCLInt(err));
        return 1;
    }

    return 0;
}

static int createBuffers(NBodyCLInfo* ci,
                         NBodyCLMem* cm,
                         NBodyCtx* ctx,
                         NBodyState* st)
{
    const size_t nbody    = ctx->model.nbody;
    const size_t bodysize = nbody * sizeof(body);
    const size_t accSize  = sizeof(vector) * nbody;
    cl_int err;

    cm->nbctx = clCreateBuffer(ci->clctx,
                               CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               sizeof(NBodyCtx),
                               ctx,
                               &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating context buffer: %s\n", showCLInt(err));
        return 1;
    }

    /* FIXME: I know this is completey wrong, and only "sort of" works
     * for the cpu, in that it doesn't crash. */
    cm->root = clCreateBuffer(ci->clctx,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                              sizeof(nodeptr),
                              &st->tree.root,
                              &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating root buffer: %s\n", showCLInt(err));
        return 1;
    }

    cm->bodies = clCreateBuffer(ci->clctx,
                                CL_MEM_READ_ONLY,
                                bodysize,
                                NULL,
                                &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating body buffer: %s\n", showCLInt(err));
        return 1;
    }

    cm->acc = clCreateBuffer(ci->clctx, CL_MEM_WRITE_ONLY, accSize, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating acceleration buffer: %s\n", showCLInt(err));
        return 1;
    }

    return 0;
}

/* Query one device specified by type, create a context, command
 * queue, and compile the kernel for it
 * Returns: non-zero on failure
 *  TODO: Multiple device support
 *  TODO: Caching of compiled binaries
 */
static int getCLInfo(NBodyCLInfo* ci, NBodyCtx* ctx, cl_device_type type)
{
    cl_int err;
    cl_uint maxComputeUnits, clockFreq;
    cl_ulong memSize;
    char* compileDefinitions;

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

    ci->prog = clCreateProgramWithSource(ci->clctx, 1, &cl_gravmap_src, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating program: %s\n", showCLInt(err));
        return 1;
    }

    asprintf(&compileDefinitions,
             "-I/Users/matt/src/milkywayathome_client/nbody/include "
             "-DNBODY_OPENCL=1 "
             "-DDOUBLEPREC=0 "
             "-DSPHERICALTYPE=%d -DDISKTYPE=%d -DHALOTYPE=%d ",
             ctx->pot.sphere[0].type,
             ctx->pot.disk.type,
             ctx->pot.halo.type);

    err = clBuildProgram(ci->prog,
                         1,
                         &ci->dev,
                         compileDefinitions,
                         NULL,
                         NULL);
    free(compileDefinitions);
    if (err != CL_SUCCESS)
    {
        char buildLog[BUFSIZE] = "";
        cl_int infoErr;
        cl_build_status stat;
        size_t failSize;

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
        return 1;
    }

    ci->kern = clCreateKernel(ci->prog, "gravMap", &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating kernel: %s\n", showCLInt(err));
        return 1;
    }

    clUnloadCompiler();

    return 0;
}

inline static int enqueueGravMap(NBodyCLInfo* ci, const size_t nbody)
{
    const size_t global[] = { nbody };
    cl_int err;
    err = clEnqueueNDRangeKernel(ci->queue,
                                 ci->kern,
                                 1,
                                 NULL, global, NULL,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Error enqueueing kernel execution: %s\n", showCLInt(err));
        return 1;
    }

    return 0;
}

/* Debugging */
static void printResults(vector* vs, const char* name, const size_t n)
{
    size_t i;
    char* vstr;
    for (i = 0; i < n; ++i)
    {
        vstr = showVector(vs[i]);
        printf("%s[%zu]: %s\n", name, i, vstr);
        free(vstr);
    }
}

inline static int prepareExec(NBodyCLInfo* ci, NBodyCLMem* cm, const NBodyCtx* ctx, NBodyState* st)
{
    int rc = 0;
    rc |= createBuffers(ci, cm, ctx, st);
    rc |= nbodySetKernelArgs(ci, cm, ctx->model.nbody);
    return rc;
}

/* Write bodies into the CL buffer */
inline static void writeBodies(NBodyState* st, const size_t bodySize)
{
    cl_int err;
    err = clEnqueueWriteBuffer(st->ci.queue,
                               st->cm.bodies,
                               CL_TRUE,
                               0, bodySize, st->bodytab,
                               0, NULL, NULL);

    if (err != CL_SUCCESS)
        fail("Error writing CL body buffer\n");
}

/* Read resulting acceleration vectors from the CL buffer */
static inline void readAccels(NBodyState* st, const size_t accSize)
{
    cl_int err;
    err = clEnqueueReadBuffer(st->ci.queue,
                              st->cm.acc,
                               CL_TRUE,
                               0, accSize, st->acctab,
                               0, NULL, NULL);

    if (err != CL_SUCCESS)
        fail("Error reading CL acceleration buffer\n");
}

int setupNBodyCL(NBodyCtx* ctx, NBodyState* st)
{
    if (getCLInfo(&st->ci, ctx, CL_DEVICE_TYPE_CPU))
        fail("Failed to setup OpenCL device\n");

    /* Kernel arguments only need to be set once */
    if (prepareExec(&st->ci, &st->cm, ctx, st))
        fail("Error running prepareExec\n");

    return 0;
}

void cleanupNBodyCL(NBodyState* st)
{
    destroyNBodyCLInfo(&st->ci);
    releaseNBodyCLMem(&st->cm);
}

void gravMapCL(const NBodyCtx* ctx, NBodyState* st)
{
    const size_t nbody = ctx->model.nbody;

    makeTree(ctx, st);

    writeBodies(st, nbody * sizeof(body));
    enqueueGravMap(&st->ci, nbody);
    readAccels(st, sizeof(vector) * nbody);
}

