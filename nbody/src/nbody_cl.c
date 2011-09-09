/*
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "milkyway_cl.h"
#include "milkyway_util.h"
#include "nbody_cl.h"
#include "nbody_show.h"
#include "nbody_util.h"

#define MAXDEPTH 26

typedef struct
{
    cl_mem pos[3];
    cl_mem vel[3];
    cl_mem acc[3];
    cl_mem max[3];
    cl_mem min[3];
    cl_mem masses;
    cl_mem treeStatus;

    cl_mem start; /* TODO: We can reuse other buffers with this later to save memory */
    cl_mem count;
    cl_mem child;
    cl_mem sort;

    cl_mem critRadii; /* Used by the alternative cell opening criterion.
                         Unnecessary for BH86.
                         BH86 will be the fastest option since it won't need to load from this
                       */

    cl_mem debug;
} NBodyBuffers;

/* CHECKME: Padding between these fields might be a good idea */
typedef struct NBODY_ALIGN
{
    real radius;
    int bottom;
    int maxDepth;
    int errorCode;
    unsigned int blkCnt;
} TreeStatus;

typedef struct
{
    real f[32];
    int i[64];
} Debug;

typedef struct
{
    size_t factors[6];
    size_t threads[6];
    double timings[6];        /* In a single iteration */
    double chunkTimings[6];   /* Average time per chunk */
    double kernelTimings[6];  /* Running totals */

    size_t global[6];
    size_t local[6];
} NBodyWorkSizes;

static NBodyWorkSizes _workSizes;

static void printNBodyWorkSizes(const NBodyWorkSizes* ws)
{
    mw_printf("\n"
              "Kernel launch sizes:\n"
              "  Bounding box kernel:  "ZU", "ZU"\n"
              "  Tree build kernel:    "ZU", "ZU"\n"
              "  Summarization kernel: "ZU", "ZU"\n"
              "  Sort kernel:          "ZU", "ZU"\n"
              "  Force kernel:         "ZU", "ZU"\n"
              "  Integration kernel:   "ZU", "ZU"\n"
              "\n",
              ws->global[0], ws->local[0],
              ws->global[1], ws->local[1],
              ws->global[2], ws->local[2],
              ws->global[3], ws->local[3],
              ws->global[4], ws->local[4],
              ws->global[5], ws->local[5]
        );
}

static cl_bool setWorkSizes(NBodyWorkSizes* ws, const DevInfo* di)
{
    cl_uint blocks = di->maxCompUnits;

    ws->global[0] = ws->threads[0] * ws->factors[0] * blocks;
    ws->local[0] = ws->threads[0];

    ws->global[1] = ws->threads[1] * ws->factors[1] * blocks;
    ws->local[1] = ws->threads[1];

    ws->global[2] = ws->threads[2] * ws->factors[2] * blocks;
    ws->local[2] = ws->threads[2];

    ws->global[3] = ws->threads[3] * ws->factors[3] * blocks;
    ws->local[3] = ws->threads[3];

    ws->global[4] = ws->threads[4] * ws->factors[4] * blocks;
    ws->local[4] = ws->threads[4];

    ws->global[5] = ws->threads[5] * ws->factors[5] * blocks;
    ws->local[5] = ws->threads[5];

    return CL_FALSE;
}

/* Return CL_TRUE if some error */
static cl_bool setThreadCounts(NBodyWorkSizes* ws, const DevInfo* di)
{
    ws->factors[0] = 1;
    ws->factors[1] = 2;
    ws->factors[2] = 1;
    ws->factors[3] = 1;
    ws->factors[4] = 1;
    ws->factors[5] = 1;

    if (di->devType == CL_DEVICE_TYPE_CPU)
    {
        ws->threads[0] = 1;
        ws->threads[1] = 1;
        ws->threads[2] = 1;
        ws->threads[3] = 1;
        ws->threads[4] = 1;
        ws->threads[5] = 1;
    }
    else if (computeCapabilityIs(di, 1, 3))
    {
        ws->threads[0] = 512;
        ws->threads[1] = 288;
        ws->threads[2] = 256;
        ws->threads[3] = 512;
        ws->threads[4] = 384;
        ws->threads[5] = 512;
    }
    else if (minComputeCapabilityCheck(di, 2, 0))
    {
        ws->threads[0] = 512;
        ws->threads[1] = 1024;
        ws->threads[2] = 1024;
        ws->threads[3] = 256;
        ws->threads[4] = 256;
        ws->threads[5] = 512;
    }
    else
    {
        ws->threads[0] = 256;
        ws->threads[1] = 256;
        ws->threads[2] = 256;
        ws->threads[3] = 256;
        ws->threads[4] = 256;
        ws->threads[5] = 256;
    }

    return CL_FALSE;
}


static void* mapBuffer(CLInfo* ci, cl_mem mem, cl_map_flags flags, size_t size)
{
    return clEnqueueMapBuffer(ci->queue, mem, CL_TRUE, flags, 0,
                              size,
                              0, NULL, NULL, NULL);

}

static void printDebug(const Debug* d)
{
    int i;
    for (i = 0; i < 32; ++i)
    {
        mw_printf("Debug.int[%d] = %d\n", i, d->i[i]);
    }

    for (i = 0; i < 32; ++i)
    {
        mw_printf("Debug.float[%d] = %.15f\n", i, d->f[i]);
    }
}

static void printTreeStatus(const TreeStatus* ts)
{
    mw_printf("TreeStatus = {\n"
              "  radius    = %.15f\n"
              "  bottom    = %d\n"
              "  maxDepth  = %d\n"
              "  errorCode = %d\n"
              "  blckCnt   = %u\n"
              "}\n",
              ts->radius,
              ts->bottom,
              ts->maxDepth,
              ts->errorCode,
              ts->blkCnt);
}


static struct
{
    cl_kernel boundingBox;
    cl_kernel buildTree;
    cl_kernel summarization;
    cl_kernel sort;
    cl_kernel forceCalculation;
    cl_kernel integration;
} kernels;

#define NKERNELS 6

/* Set arguments (idx, idx + 1, idx + 2) to the buffers in mem[3] */
static cl_int setMemArrayArgs(CLInfo* ci, cl_kernel kern, cl_mem mem[3], cl_uint idx)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;

    for (i = 0; i < 3; ++i)
    {
        err |= clSetKernelArg(kern, idx + i, sizeof(cl_mem), &mem[i]);
    }

    return err;
}


static cl_int setKernelArguments(cl_kernel kern, CLInfo* ci, NBodyBuffers* nbb)
{
    cl_int err = CL_SUCCESS;
    cl_int step = 0;

    err |= setMemArrayArgs(ci, kern, nbb->pos, 0);
    err |= setMemArrayArgs(ci, kern, nbb->vel, 3);
    err |= setMemArrayArgs(ci, kern, nbb->acc, 6);
    err |= setMemArrayArgs(ci, kern, nbb->max, 9);
    err |= setMemArrayArgs(ci, kern, nbb->min, 12);

    err |= clSetKernelArg(kern, 15, sizeof(cl_mem), &nbb->masses);

    err |= clSetKernelArg(kern, 16, sizeof(cl_mem), &nbb->start);
    err |= clSetKernelArg(kern, 17, sizeof(cl_mem), &nbb->count);
    err |= clSetKernelArg(kern, 18, sizeof(cl_mem), &nbb->child);
    err |= clSetKernelArg(kern, 19, sizeof(cl_mem), &nbb->sort);
    err |= clSetKernelArg(kern, 20, sizeof(cl_mem), &nbb->treeStatus);


    /* Set the step to be 0 here. Only the force calculation kernel
     * actually cares. We'll set it before then. */
    err |= clSetKernelArg(kern, 21, sizeof(cl_int), &step);

    /* Value also doesn't matter */
    err |= clSetKernelArg(kern, 22, sizeof(cl_int), &step);


    if (nbb->critRadii)
    {
        err |= clSetKernelArg(kern, 23, sizeof(cl_mem), &nbb->critRadii);
    }
    else
    {
        /* Just be anything. Won't be used, it just needs to be not null to not error */
        err |= clSetKernelArg(kern, 23, sizeof(cl_mem), &nbb->treeStatus);
    }

    err |= clSetKernelArg(kern, 24, sizeof(cl_mem), &nbb->debug);

    return err;
}

static cl_int setAllKernelArguments(CLInfo* ci, NBodyBuffers* nbb)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;
    cl_kernel* kerns = (cl_kernel*) &kernels;

    for (i = 0; i < NKERNELS; ++i)
    {
        err |= setKernelArguments(kerns[i], ci, nbb);
    }

    return err;
}


static cl_int createKernel(cl_kernel* kern, CLInfo* ci, const char* name)
{
    cl_int err = CL_SUCCESS;

    *kern = clCreateKernel(ci->prog, name, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to create kernel '%s'", err, name);
        return err;
    }

    return err;
}


static cl_int createKernels(CLInfo* ci)
{
    cl_int err = CL_SUCCESS;

    err |= createKernel(&kernels.boundingBox, ci, "boundingBox");
    err |= createKernel(&kernels.buildTree, ci, "buildTree");
    err |= createKernel(&kernels.summarization, ci, "summarization");
    err |= createKernel(&kernels.sort, ci, "sort");
    err |= createKernel(&kernels.forceCalculation, ci, "forceCalculation");
    err |= createKernel(&kernels.integration, ci, "integration");

    return err;
}

static cl_int clReleaseKernel_quiet(cl_kernel kern)
{
    if (!kern)
        return CL_SUCCESS;

    return clReleaseKernel(kern);
}

static cl_int releaseKernels()
{
    cl_int err = CL_SUCCESS;

    err |= clReleaseKernel_quiet(kernels.boundingBox);
    err |= clReleaseKernel_quiet(kernels.buildTree);
    err |= clReleaseKernel_quiet(kernels.summarization);
    err |= clReleaseKernel_quiet(kernels.sort);
    err |= clReleaseKernel_quiet(kernels.forceCalculation);
    err |= clReleaseKernel_quiet(kernels.integration);

    if (err != CL_SUCCESS)
        mwCLWarn("Error releasing kernels", err);

    return err;
}

static cl_uint findNNode(const DevInfo* di, cl_int nbody)
{
    cl_uint nNode = 2 * nbody;

    if (nNode < 1024 * di->maxCompUnits)
        nNode = 1024 * di->maxCompUnits;
    while ((nNode & (di->warpSize - 1)) != 0)
        ++nNode;

    return nNode - 1;
}

static char* getCompileFlags(const NBodyCtx* ctx, const NBodyState* st, const DevInfo* di)
{
    char* buf;
    static NBodyWorkSizes* ws = &_workSizes;

    /* Put a space between the -D. if the thing begins with D, there's
     * an Apple OpenCL compiler bug where the D will be
     * stripped. -DDOUBLEPREC=1 will actually define OUBLEPREC */
    if (asprintf(&buf,
                 "-D DOUBLEPREC=%d "
               #if !DOUBLEPREC
                 "-cl-single-precision-constant "
               #endif

                 "-D NBODY=%d "
                 "-D NNODE=%u "
                 "-D WARPSIZE=%u "

                 "-D NOSORT=%d "

                 "-D THREADS1="ZU" "
                 "-D THREADS2="ZU" "
                 "-D THREADS3="ZU" "
                 "-D THREADS4="ZU" "
                 "-D THREADS5="ZU" "
                 "-D THREADS6="ZU" "
                 "-D MAXDEPTH=%u "

                 "-D TIMESTEP=%a "
                 "-D EPS2=%a "
                 "-D THETA=%a "

                 "-D NEWCRITERION=%d "
                 "-D SW93=%d "
                 "-D BH86=%d "
                 "-D EXACT=%d "

                 "%s ",
                 DOUBLEPREC,

                 st->nbody,
                 findNNode(di, st->nbody),
                 di->warpSize,

                 (di->devType == CL_DEVICE_TYPE_CPU),

                 ws->threads[0],
                 ws->threads[1],
                 ws->threads[2],
                 ws->threads[3],
                 ws->threads[4],
                 ws->threads[5],
                 MAXDEPTH,

                 ctx->timestep,
                 ctx->eps2,
                 ctx->theta,

                 ctx->criterion == NewCriterion,
                 ctx->criterion == SW93,
                 ctx->criterion == BH86,
                 ctx->criterion == Exact,
                 hasNvidiaCompilerFlags(di) ? "-cl-nv-verbose" : ""
            ) < 1)
    {
        mw_printf("Error getting compile flags\n");
        return NULL;
    }

    return buf;
}

static cl_int loadKernels(CLInfo* ci, const NBodyCtx* ctx, const NBodyState* st)
{
    cl_int err = CL_SUCCESS;
    char* src = NULL;
    char* compileFlags = NULL;

    src = mwReadFile("kernels/nbody_kernels.cl");
    if (!src)
        return MW_CL_ERROR;

    compileFlags = getCompileFlags(ctx, st, &ci->di);
    assert(compileFlags);
    err = mwSetProgramFromSrc(ci, (const char**) &src, 1, compileFlags);
    if (err != CL_SUCCESS)
    {
        mw_printf("Failed build flags: %s\n", compileFlags);
    }

    free(src);
    free(compileFlags);

    return err;
}

/* Return CL_FALSE if device isn't capable of running this */
static cl_bool nbodyCheckDevCapabilities(const DevInfo* di, const NBodyCtx* ctx, NBodyState* st)
{
    return CL_TRUE;
}

static cl_int debug(CLInfo* ci, NBodyBuffers* nbb)
{
    cl_int err;
    Debug d;

    memset(&d, 0, sizeof(d));

    err = clEnqueueReadBuffer(ci->queue,
                              nbb->debug,
                              CL_TRUE,
                              0, sizeof(d), &d,
                              0, NULL, NULL);

    mw_printf("Debug:\n");
    printDebug(&d);
    return err;
}

static cl_int readTreeStatus(TreeStatus* tc, CLInfo* ci, NBodyBuffers* nbb)
{
    assert(tc);

    return clEnqueueReadBuffer(ci->queue,
                               nbb->treeStatus,
                               CL_TRUE,
                               0, sizeof(*tc), tc,
                               0, NULL, NULL);
}

static cl_int printBuffer(CLInfo* ci, cl_mem mem, size_t n, const char* name, int type)
{
    size_t i;
    void* p;

    p = mapBuffer(ci, mem, CL_MAP_READ, n * (type == 0 ? sizeof(real) : sizeof(int)));
    if (!p)
        return MW_CL_ERROR;

    if (type == 0)
    {
        const real* pr = (const real*) p;
        for (i = 0; i < n; ++i)
        {
            mw_printf("%s["ZU"] = %.15f\n", name, i, pr[i]);
        }
    }
    else
    {
        const int* pi = (const int*) p;
        for (i = 0; i < n; ++i)
        {
            mw_printf("%s["ZU"] = %d\n", name, i, pi[i]);
        }
    }

    return clEnqueueUnmapMemObject(ci->queue, mem, p, 0, NULL, NULL);
}

static void stdDebugPrint(CLInfo* ci, NBodyState* st, NBodyBuffers* nbb)
{
    mw_printf("--------------------------------------------------------------------------------\n");
    cl_int err;
    cl_uint nNode = findNNode(&ci->di, st->nbody);
    mw_printf("BEGIN CHILD\n");
    printBuffer(ci, nbb->child, NSUB * (nNode + 1), "child", 1);
    mw_printf("BEGIN START\n");
    printBuffer(ci, nbb->start, nNode, "start", 1);

    mw_printf("BEGIN MASS\n");
    printBuffer(ci, nbb->masses, nNode + 1, "mass", 0);

    {
        TreeStatus tc;
        memset(&tc, 0, sizeof(tc));
        err = readTreeStatus(&tc, ci, nbb);
        if (err != CL_SUCCESS)
            mwCLWarn("Reading tree status failed\n", err);
        else
            printTreeStatus(&tc);
    }

    debug(ci, nbb);
    mw_printf("--------------------------------------------------------------------------------\n");
}

/* Check the error code */
static cl_bool checkKernelErrorCode(CLInfo* ci, NBodyBuffers* nbb)
{
    cl_int err;
    TreeStatus ts;

    err = readTreeStatus(&ts, ci, nbb);
    if (err != CL_SUCCESS)
        return CL_TRUE;

    if (ts.errorCode != 0)
    {
        mw_printf("Kernel reported error: %d\n", ts.errorCode);
        return CL_TRUE;
    }

    return CL_FALSE;
}

static cl_double waitReleaseEventWithTime(cl_event ev)
{
    cl_double t;
    cl_int err;

    err = clWaitForEvents(1, &ev);
    if (err != CL_SUCCESS)
        return 0.0;

    t = mwEventTimeMS(ev);

    err = clReleaseEvent(ev);
    if (err != CL_SUCCESS)
        return 0.0;

    return t;
}

/* If kernels are taking too long, try running with fewer threads if
 * possible. Return CL_TRUE if anything changed
 */
static cl_bool reevaluateWorkSizes(NBodyWorkSizes* ws)
{
    cl_double buildTreeThreshold = 33.0;
    cl_double forceThreshold = 33.0;

    if (ws->chunkTimings[1] > buildTreeThreshold)
    {
        /* Try fewer threads if possible */
    }

    if (ws->chunkTimings[4] > forceThreshold)
    {
        /* Try fewer threads if possible */
    }

    return CL_FALSE;
}

static cl_int stepSystemCL(CLInfo* ci, const NBodyCtx* ctx, NBodyState* st)
{
    cl_int err;
    cl_uint i;
    int upperBound;
    size_t offset[1];
    size_t chunk, nChunk;
    cl_bool ignoreResponsive = CL_FALSE;
    cl_event boxEv, sumEv, sortEv, integrateEv;
    NBodyWorkSizes* ws = &_workSizes;

    memset(ws->timings, 0, sizeof(ws->timings));

    err = clSetKernelArg(kernels.forceCalculation, 21, sizeof(int), &st->step);
    if (err != CL_SUCCESS)
        return err;

    err = clEnqueueNDRangeKernel(ci->queue, kernels.boundingBox, 1,
                                 NULL, &ws->global[0], &ws->local[0],
                                 0, NULL, &boxEv);
    if (err != CL_SUCCESS)
        return err;

    ws->timings[0] += waitReleaseEventWithTime(boxEv);

    nChunk     = ignoreResponsive ?         1 : mwDivRoundup((size_t) st->nbody, ws->global[0]);
    upperBound = ignoreResponsive ? st->nbody : (int) ws->global[0];
    offset[0] = 0;
    for (chunk = 0; chunk < nChunk; ++chunk)
    {
        cl_event ev;

        if (upperBound > st->nbody)
            upperBound = st->nbody;

        err = clSetKernelArg(kernels.buildTree, 22, sizeof(int), &upperBound);
        if (err != CL_SUCCESS)
            return err;

        err = clEnqueueNDRangeKernel(ci->queue, kernels.buildTree, 1,
                                     offset, &ws->global[1], &ws->local[1],
                                     0, NULL, &ev);
        if (err != CL_SUCCESS)
            return err;

        ws->timings[1] += waitReleaseEventWithTime(ev);

        upperBound += (int) ws->global[1];
        offset[0] += ws->global[1];
    }
    ws->chunkTimings[1] = ws->timings[1] / (double) nChunk;

    err = clEnqueueNDRangeKernel(ci->queue, kernels.summarization, 1,
                                 NULL, &ws->global[2], &ws->local[2],
                                 0, NULL, &sumEv);
    if (err != CL_SUCCESS)
        return err;

    /* FIXME: This does not work unless ALL of the threads are
     * launched at once. This may be bad when we need
     * responsiveness. This also means it will always hang with
     * CPUs */
    err = clEnqueueNDRangeKernel(ci->queue, kernels.sort, 1,
                                 NULL, &ws->global[3], &ws->local[3],
                                 0, NULL, &sortEv);
    if (err != CL_SUCCESS)
        return err;

    ws->timings[3] += waitReleaseEventWithTime(sortEv);

    nChunk     = ignoreResponsive ?         1 : mwDivRoundup((size_t) st->nbody, ws->global[4]);
    upperBound = ignoreResponsive ? st->nbody : (int) ws->global[4];
    offset[0] = 0;
    for (chunk = 0; chunk < nChunk; ++chunk)
    {
        cl_event ev;
        double t;

        if (upperBound > st->nbody)
            upperBound = st->nbody;

        err = clSetKernelArg(kernels.forceCalculation, 22, sizeof(int), &upperBound);
        if (err != CL_SUCCESS)
            return err;
        err = clEnqueueNDRangeKernel(ci->queue, kernels.forceCalculation, 1,
                                     offset, &ws->global[4], &ws->local[4],
                                     0, NULL, &ev);
        if (err != CL_SUCCESS)
            return err;

        t = waitReleaseEventWithTime(ev);

        ws->timings[4] += t;
        upperBound += (int) ws->global[4];
        offset[0] += ws->global[4];
    }
    ws->chunkTimings[4] = ws->timings[4] / (double) nChunk;

    err = clEnqueueNDRangeKernel(ci->queue, kernels.integration, 1,
                                 NULL, &ws->global[5], &ws->local[5],
                                 0, NULL, &integrateEv);
    if (err != CL_SUCCESS)
        return err;

    ws->timings[5] += waitReleaseEventWithTime(integrateEv);

    err = clFinish(ci->queue);
    if (err != CL_SUCCESS)
        return err;

    mw_printf("Step %d:\n"
              "  boundingBox:      %15f ms\n"
              "  buildTree:        %15f ms%15f ms\n"
              "  summarization:    %15f ms\n"
              "  sort:             %15f ms\n"
              "  forceCalculation: %15f ms%15f ms\n"
              "  integration:      %15f ms\n"
              "\n",
              st->step,
              ws->timings[0],
              ws->timings[1], ws->chunkTimings[1],
              ws->timings[2],
              ws->timings[3],
              ws->timings[4], ws->chunkTimings[4],
              ws->timings[5]);

    for (i = 0; i < 6; ++i) /* Add timings to running totals */
    {
        ws->kernelTimings[i] += ws->timings[i];
    }

    return err;
}

static cl_int nbodyMainLoop(CLInfo* ci, const NBodyCtx* ctx, NBodyState* st, NBodyBuffers* nbb)
{
    cl_int err = CL_SUCCESS;
    const real tstop = ctx->timeEvolve - ctx->timestep / 1024.0;

    st->tnow = 0;
    st->step = 0;
    while (err == CL_SUCCESS && st->tnow < tstop)
    {
        if (checkKernelErrorCode(ci, nbb))
        {
            err = MW_CL_ERROR;
            break;
        }

        mw_printf("Running step %d (%f%%)\n",
                  st->step,
                  100.0 * st->tnow / tstop);
        err = stepSystemCL(ci, ctx, st);
        st->tnow += ctx->timestep;
        st->step++;
    }
    mw_printf("Broke on step %d, %s\n", st->step - 1, showCLInt(err));

    return err;
}

/* This is dumb and errors if mem isn't set */
static cl_int clReleaseMemObject_quiet(cl_mem mem)
{
    if (!mem)
        return CL_SUCCESS;

    return clReleaseMemObject(mem);
}

static cl_int releaseBuffers(NBodyBuffers* nbb)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;

    for (i = 0; i < 3; ++i)
    {
        err |= clReleaseMemObject_quiet(nbb->pos[i]);
        err |= clReleaseMemObject_quiet(nbb->vel[i]);
        err |= clReleaseMemObject_quiet(nbb->acc[i]);
        err |= clReleaseMemObject_quiet(nbb->max[i]);
        err |= clReleaseMemObject_quiet(nbb->min[i]);
    }

    err |= clReleaseMemObject_quiet(nbb->masses);
    err |= clReleaseMemObject_quiet(nbb->treeStatus);
    err |= clReleaseMemObject_quiet(nbb->start);
    err |= clReleaseMemObject_quiet(nbb->count);
    err |= clReleaseMemObject_quiet(nbb->child);
    err |= clReleaseMemObject_quiet(nbb->sort);

    if (err != CL_SUCCESS)
        mwCLWarn("Error releasing buffers", err);

    return err;
}

static cl_int setInitialTreeStatus(CLInfo* ci, NBodyBuffers* nbb)
{
    TreeStatus iniTreeStatus;

    memset(&iniTreeStatus, 0, sizeof(iniTreeStatus));

    iniTreeStatus.radius = 0;
    iniTreeStatus.bottom = 0;
    iniTreeStatus.maxDepth = 1;
    iniTreeStatus.errorCode = 0;
    iniTreeStatus.blkCnt = 0;

    printTreeStatus(&iniTreeStatus);

    return clEnqueueWriteBuffer(ci->queue,
                                nbb->treeStatus,
                                CL_TRUE,
                                0, sizeof(TreeStatus), &iniTreeStatus,
                                0, NULL, NULL);
}

static cl_uint findInc(cl_uint warpSize, cl_uint nbody)
{
    return (nbody + warpSize - 1) & (-warpSize);
}

static cl_int createBuffers(const NBodyCtx* ctx, NBodyState* st, CLInfo* ci, NBodyBuffers* nbb)
{
    cl_uint i;
    cl_uint nNode = findNNode(&ci->di, st->nbody);
    cl_uint inc = findInc(ci->di.warpSize, st->nbody);

    mw_printf("NNODE = %u, nbody = %d\n"
              "(NNODE + 1) * NSUB = %u\n"
              "inc %u\n",
              nNode, st->nbody, NSUB * (nNode + 1), inc);

    for (i = 0; i < 3; ++i)
    {
        nbb->pos[i] = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
        nbb->vel[i] = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(real));
        nbb->acc[i] = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(real));
        nbb->min[i] = mwCreateZeroReadWriteBuffer(ci, ci->di.maxCompUnits * sizeof(real));
        nbb->max[i] = mwCreateZeroReadWriteBuffer(ci, ci->di.maxCompUnits * sizeof(real));

        if (!nbb->pos[i] || !nbb->vel[i] || !nbb->acc[i] || !nbb->min[i] || !nbb->max[i])
        {
            return MW_CL_ERROR;
        }
    }

    nbb->masses = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
    nbb->treeStatus = mwCreateZeroReadWriteBuffer(ci, sizeof(TreeStatus));

    nbb->start = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(cl_int));
    nbb->count = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(cl_int));
    nbb->sort = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(cl_int));
    nbb->child = mwCreateZeroReadWriteBuffer(ci, NSUB * (nNode + 1) * sizeof(cl_int));

    nbb->debug = mwCreateZeroReadWriteBuffer(ci, sizeof(Debug));

    if (!nbb->masses || !nbb->treeStatus || !nbb->start || !nbb->count || !nbb->sort || !nbb->child)
        return MW_CL_ERROR;

    if (ctx->criterion == SW93 || ctx->criterion == NewCriterion)
    {
        /* This only is for cells, so we could subtract nbody if we wanted */
        nbb->critRadii = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
        if (!nbb->critRadii)
            return MW_CL_ERROR;
    }

    return CL_SUCCESS;
}

static cl_int mapBodies(real* pos[3], real* vel[3], real** mass, NBodyBuffers* nbb, CLInfo* ci, cl_map_flags flags, NBodyState* st)
{
    cl_uint i;

    for (i = 0; i < 3; ++i)
    {
        pos[i] = (real*) mapBuffer(ci, nbb->pos[i], flags, st->nbody * sizeof(real));
        vel[i] = (real*) mapBuffer(ci, nbb->vel[i], flags, st->nbody * sizeof(real));
        if (!pos[i] || !vel[i])
        {
            return MW_CL_ERROR;
        }
    }

    *mass = (real*) mapBuffer(ci, nbb->masses, flags, st->nbody * sizeof(real));

    return CL_SUCCESS;
}

static cl_int unmapBodies(real* pos[3], real* vel[3], real* mass, NBodyBuffers* nbb, CLInfo* ci)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;

    for (i = 0; i < 3; ++i)
    {
        if (pos[i])
            err |= clEnqueueUnmapMemObject(ci->queue, nbb->pos[i], pos[i], 0, NULL, NULL);

        if (vel[i])
            err |= clEnqueueUnmapMemObject(ci->queue, nbb->vel[i], vel[i], 0, NULL, NULL);
    }

    if (mass)
        err |= clEnqueueUnmapMemObject(ci->queue, nbb->masses, mass, 0, NULL, NULL);

    return err;
}

/* If last parameter is true, copy to the buffers. if false, copy from the buffers */
static cl_int marshalBodies(NBodyBuffers* nbb, CLInfo* ci, NBodyState* st, cl_bool marshalIn)
{
    cl_int i;
    cl_int err = CL_SUCCESS;
    const Body* b;
    real* pos[3] = { NULL, NULL, NULL };
    real* vel[3] = { NULL, NULL, NULL };
    real* mass = NULL;
    cl_map_flags flags = marshalIn ? CL_MAP_WRITE : CL_MAP_READ;

    err = mapBodies(pos, vel, &mass, nbb, ci, flags, st);
    if (err != CL_SUCCESS)
    {
        unmapBodies(pos, vel, mass, nbb, ci);
        return err;
    }

    if (marshalIn)
    {
        for (i = 0, b = st->bodytab; b < st->bodytab + st->nbody; ++i, ++b)
        {
            pos[0][i] = X(Pos(b));
            pos[1][i] = Y(Pos(b));
            pos[2][i] = Z(Pos(b));

            vel[0][i] = X(Vel(b));
            vel[1][i] = Y(Vel(b));
            vel[2][i] = Z(Vel(b));

            mass[i] = Mass(b);
        }
    }
    else
    {
        for (i = 0, b = st->bodytab; b < st->bodytab + st->nbody; ++i, ++b)
        {
            X(Pos(b)) = pos[0][i];
            Y(Pos(b)) = pos[1][i];
            Z(Pos(b)) = pos[2][i];

            X(Vel(b)) = vel[0][i];
            Y(Vel(b)) = vel[1][i];
            Z(Vel(b)) = vel[2][i];

            Mass(b) = mass[i];
        }
    }

    return unmapBodies(pos, vel, mass, nbb, ci);
}

static void printKernelTimings(const DevInfo* di, const NBodyState* st)
{

    double totalTime = 0.0;
    cl_uint i;
    double nStep = (double) st->step;
    const double* kernelTimings = _workSizes.kernelTimings;

    for (i = 0; i < 6; ++i)
    {
        totalTime += kernelTimings[i];
    }

    mw_printf("\n--------------------------------------------------------------------------------\n"
              "Total timing over %d steps:\n"
              "                         Average             Total            Fraction\n"
              "                    ----------------   ----------------   ----------------\n"
              "  boundingBox:      %16f   %16f   %15.4f%%\n"
              "  buildTree:        %16f   %16f   %15.4f%%\n"
              "  summarization:    %16f   %16f   %15.4f%%\n"
              "  sort:             %16f   %16f   %15.4f%%\n"
              "  forceCalculation: %16f   %16f   %15.4f%%\n"
              "  integration:      %16f   %16f   %15.4f%%\n"
              "\n--------------------------------------------------------------------------------\n"
              "\n",
              st->step,
              kernelTimings[0] / nStep, kernelTimings[0], 100.0 * kernelTimings[0] / totalTime,
              kernelTimings[1] / nStep, kernelTimings[1], 100.0 * kernelTimings[1] / totalTime,
              kernelTimings[2] / nStep, kernelTimings[2], 100.0 * kernelTimings[2] / totalTime,
              kernelTimings[3] / nStep, kernelTimings[3], 100.0 * kernelTimings[3] / totalTime,
              kernelTimings[4] / nStep, kernelTimings[4], 100.0 * kernelTimings[4] / totalTime,
              kernelTimings[5] / nStep, kernelTimings[5], 100.0 * kernelTimings[5] / totalTime
        );
}

static void setCLRequestFromFlags(CLRequest* clr, const NBodyFlags* nbf)
{
    clr->platform = nbf->platform;
    clr->devNum = nbf->devNum;
    clr->verbose = TRUE;
    clr->enableCheckpointing = FALSE;
}

/* Setup kernels and buffers */
static cl_int setupExec(CLInfo* ci, const NBodyCtx* ctx, NBodyState* st, NBodyBuffers* nbb)
{
    cl_int err;

    err = loadKernels(ci, ctx, st);
    if (err != CL_SUCCESS)
        return err;

    err = createKernels(ci);
    if (err != CL_SUCCESS)
        return err;

    err = createBuffers(ctx, st, ci, nbb);
    if (err != CL_SUCCESS)
        return err;

    err = setInitialTreeStatus(ci, nbb);
    if (err != CL_SUCCESS)
        return err;

    err = setAllKernelArguments(ci, nbb);
    if (err != CL_SUCCESS)
        return err;

    return CL_SUCCESS;
}

NBodyStatus runSystemCL(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    cl_int err = CL_SUCCESS;
    CLInfo ci;
    CLRequest clr;
    NBodyBuffers nbb;

    memset(&ci, 0, sizeof(ci));
    memset(&clr, 0, sizeof(clr));
    memset(&nbb, 0, sizeof(nbb));

    setCLRequestFromFlags(&clr, nbf);
    clr.enableProfiling = TRUE;

    err = mwSetupCL(&ci, &clr);
    if (err != CL_SUCCESS)
        return NBODY_ERROR;

    if (!nbodyCheckDevCapabilities(&ci.di, ctx, st))
        return NBODY_ERROR;

    if (setThreadCounts(&_workSizes, &ci.di) || setWorkSizes(&_workSizes, &ci.di))
        return NBODY_ERROR;

    printNBodyWorkSizes(&_workSizes);

    err = setupExec(&ci, ctx, st, &nbb);
    if (err != CL_SUCCESS)
        goto fail;

    err = marshalBodies(&nbb, &ci, st, CL_TRUE);
    if (err != CL_SUCCESS)
        goto fail;

    err = nbodyMainLoop(&ci, ctx, st, &nbb);
    if (err != CL_SUCCESS)
        goto fail;

    err = marshalBodies(&nbb, &ci, st, CL_FALSE);
    if (err != CL_SUCCESS)
        goto fail;

    //printBodies(st->bodytab, st->nbody);

    printKernelTimings(&ci.di, st);

fail:
    debug(&ci, &nbb);

    mwDestroyCLInfo(&ci);
    releaseKernels();
    releaseBuffers(&nbb);

    return (err == CL_SUCCESS) ? NBODY_SUCCESS : NBODY_ERROR;
}

