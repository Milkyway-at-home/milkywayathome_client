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

#define FACTOR1 1
#define FACTOR2 2
#define FACTOR3 1
#define FACTOR4 1
#define FACTOR5 2
#define FACTOR6 1

#define THREADS1 256
#define THREADS2 288
#define THREADS3 256
#define THREADS4 512
#define THREADS5 256
#define THREADS6 512

#define MAXDEPTH 26

/* FIXME: Need a way to track which bodies are ignored when they are
 * shuffled around by sorting
 */
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
    real f[16];
    int i[16];
} Debug;

static NBodyBuffers* _nbb = NULL;
static cl_uint _nNode = 0;

static void* mapBuffer(CLInfo* ci, cl_mem mem, cl_map_flags flags, size_t size)
{
    return clEnqueueMapBuffer(ci->queue, mem, CL_TRUE, flags, 0,
                              size,
                              0, NULL, NULL, NULL);

}

static void printDebug(const Debug* d)
{
    int i;
    for (i = 0; i < 16; ++i)
    {
        warn("Debug.int[%d] = %d\n", i, d->i[i]);
    }

    for (i = 0; i < 16; ++i)
    {
        warn("Debug.float[%d] = %.15f\n", i, d->f[i]);
    }
}

static void printTreeStatus(const TreeStatus* ts)
{
    warn("TreeStatus = {\n"
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

    err |= clSetKernelArg(kern, 22, sizeof(cl_mem), &nbb->debug);

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
        return err;

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
    cl_uint warpSize = 32;
    //cl_uint warpSize = 1;

    if (nNode < 1024 * di->maxCompUnits)
        nNode = 1024 * di->maxCompUnits;
    while ((nNode & (warpSize - 1)) != 0)
        ++nNode;

    return nNode - 1;
}

static char* getCompileFlags(const NBodyCtx* ctx, const NBodyState* st, const DevInfo* di)
{
    char* buf;
    cl_uint warpSize = 32;

    /* Put a space between the -D. if the thing begins with D, there's
     * an Apple OpenCL compiler bug where the D will be
     * stripped. -DDOUBLEPREC=1 will actually define OUBLEPREC */
    if (asprintf(&buf,
                 "-D DOUBLEPREC=%d "
               #if !DOUBLEPREC
                 "-cl-single-precision-constant"
               #endif

                 "-D NBODY=%d "
                 "-D NNODE=%u "
                 "-D WARPSIZE=%u "

                 "-D THREADS1=%u "
                 "-D THREADS2=%u "
                 "-D THREADS3=%u "
                 "-D THREADS4=%u "
                 "-D THREADS5=%u "
                 "-D THREADS6=%u "
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
                 warpSize,

                 THREADS1,
                 THREADS2,
                 THREADS3,
                 THREADS4,
                 THREADS5,
                 THREADS6,
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
        warn("Error getting compile flags\n");
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

    warn("Debug:\n");
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
            warn("%s["ZU"] = %.15f\n", name, i, pr[i]);
        }
    }
    else
    {
        const int* pi = (const int*) p;
        for (i = 0; i < n; ++i)
        {
            warn("%s["ZU"] = %d\n", name, i, pi[i]);
        }
    }

    return clEnqueueUnmapMemObject(ci->queue, mem, p, 0, NULL, NULL);
}

static void stdDebugPrint(CLInfo* ci, NBodyState* st)
{
    warn("--------------------------------------------------------------------------------\n");
    cl_int err;
    cl_uint nNode = findNNode(&ci->di, st->nbody);
    warn("BEGIN CHILD\n");
    printBuffer(ci, _nbb->child, NSUB * (nNode + 1), "child", 1);
    warn("BEGIN START\n");
    printBuffer(ci, _nbb->start, nNode, "start", 1);

    warn("BEGIN MASS\n");
    printBuffer(ci, _nbb->masses, nNode + 1, "mass", 0);

    {
        TreeStatus tc;
        memset(&tc, 0, sizeof(tc));
        err = readTreeStatus(&tc, ci, _nbb);
        if (err != CL_SUCCESS)
            mwCLWarn("Reading tree status failed\n", err);
        else
            printTreeStatus(&tc);
    }

    debug(ci, _nbb);
    warn("--------------------------------------------------------------------------------\n");
}

/* Check the error code and reset the block count */
static cl_bool checkKernelErrorCode(CLInfo* ci, NBodyBuffers* nbb)
{
    cl_int err;
    TreeStatus* ts;
    cl_bool rc = CL_FALSE;

    warn("Checking tree status\n");

    ts = clEnqueueMapBuffer(ci->queue,
                            nbb->treeStatus,
                            CL_TRUE,
                            CL_MAP_READ | CL_MAP_WRITE,
                            0,
                            sizeof(TreeStatus),
                            0, NULL, NULL,
                            &err);
    if (!ts)
    {
        mwCLWarn("Failed to map tree status", err);
        return CL_TRUE;
    }

    ts->blkCnt = 0;
    if (ts->errorCode != 0)
    {
        warn("Kernel reported error: %d\n", ts->errorCode);
        rc = CL_TRUE;
    }

    err = clEnqueueUnmapMemObject(ci->queue, nbb->treeStatus, ts, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to unmap tree status", err);
        return CL_TRUE;
    }

    return rc;
}


static cl_int stepSystemCL(CLInfo* ci, const NBodyCtx* ctx, NBodyState* st)
{
    cl_int err;
    size_t global[1];
    size_t local[1];
    size_t warpSize = 64;
    cl_uint blocks = ci->di.maxCompUnits;


    err = clSetKernelArg(kernels.forceCalculation, 21, sizeof(int), &st->step);
    if (err != CL_SUCCESS)
        return err;


    global[0] = THREADS1 * FACTOR1 * blocks;
    local[0] = THREADS1;
    warn("Bounding box kernel: %zu, %zu\n", global[0], local[0]);
    err = clEnqueueNDRangeKernel(ci->queue, kernels.boundingBox, 1,
                                 NULL, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
        return err;

    st->tnow += ctx->timestep;
    st->step++;

    global[0] = THREADS2 * FACTOR2 * blocks;
    local[0] = THREADS2;
    warn("Tree build kernel: %zu, %zu\n", global[0], local[0]);
    err = clEnqueueNDRangeKernel(ci->queue, kernels.buildTree, 1,
                                 NULL, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
        return err;

    err = clFinish(ci->queue);
    if (err != CL_SUCCESS)
    {
        warn("Failure on build tree kernel\n");
        return err;
    }


    global[0] = THREADS3 * FACTOR3 * blocks;
    local[0] = THREADS3;
    warn("Summarization kernel: %zu, %zu\n", global[0], local[0]);
    err = clEnqueueNDRangeKernel(ci->queue, kernels.summarization, 1,
                                 NULL, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
        return err;

    err = clFinish(ci->queue);
    if (err != CL_SUCCESS)
    {
        warn("Failure on summarization kernel\n");
        return err;
    }

    //return err;

#if 1
    global[0] = THREADS4 * FACTOR4 * blocks;
    local[0] = THREADS4;
    warn("Sort kernel: %zu, %zu\n", global[0], local[0]);
    err = clEnqueueNDRangeKernel(ci->queue, kernels.sort, 1,
                                 NULL, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
        return err;
    err = clFinish(ci->queue);
    if (err != CL_SUCCESS)
        return err;
#endif

    global[0] = THREADS5 * FACTOR5 * blocks;
    local[0] = THREADS5;
    warn("Force kernel: %zu, %zu\n", global[0], local[0]);
    err = clEnqueueNDRangeKernel(ci->queue, kernels.forceCalculation, 1,
                                 NULL, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
        return err;

    err = clFinish(ci->queue);
    if (err != CL_SUCCESS)
        return err;

    global[0] = THREADS6 * FACTOR6 * blocks;
    local[0] = THREADS6;
    warn("Integration kernel: %zu, %zu\n", global[0], local[0]);
    err = clEnqueueNDRangeKernel(ci->queue, kernels.integration, 1,
                                 NULL, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
        return err;
    clFinish(ci->queue);

    return clFinish(ci->queue);
}

static cl_int nbodyMainLoop(CLInfo* ci, const NBodyCtx* ctx, NBodyState* st, NBodyBuffers* nbb)
{
    cl_int err = CL_SUCCESS;
    const real tstop = ctx->timeEvolve - ctx->timestep / 1024.0;

    st->step = -1;

    //while (err == CL_SUCCESS && st->tnow < tstop)
    for (int i = 0; i < 20; ++i)
    {
        if (checkKernelErrorCode(ci, nbb))
        {
            err = MW_CL_ERROR;
            break;
        }

        warn("Running step %d (%f%%)\n",
             st->step,
             100.0 * st->tnow / tstop);
        err = stepSystemCL(ci, ctx, st);

        //stdDebugPrint(ci, st);
    }
    warn("Broke on step %d, %s\n", st->step, showCLInt(err));

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

static cl_int createBuffers(NBodyBuffers* nbb, CLInfo* ci, NBodyState* st)
{
    cl_uint i;
    cl_uint nNode;
    cl_int err = CL_SUCCESS;

    nNode = findNNode(&ci->di, st->nbody);
    _nNode = nNode;
    warn("NNODE = %u, nbody = %d\n", nNode + 1, st->nbody);

    warn("(NNODE + 1) * NSUB = %u\n", (nNode + 1) * NSUB);

    int warpsize = 64;
    warn("inc %d\n", (st->nbody + warpsize - 1) & (-warpsize));
    for (i = 0; i < 3; ++i)
    {
        nbb->pos[i] = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
        nbb->vel[i] = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(real));
        nbb->acc[i] = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(real));
        nbb->min[i] = mwCreateZeroReadWriteBuffer(ci, ci->di.maxCompUnits * sizeof(real));
        nbb->max[i] = mwCreateZeroReadWriteBuffer(ci, ci->di.maxCompUnits * sizeof(real));

        if (!nbb->pos[i] || !nbb->vel[i] || !nbb->acc[i] || !nbb->min[i] || !nbb->max[i])
        {
            err = MW_CL_ERROR;
            break;
        }
    }

    nbb->masses = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
    nbb->treeStatus = mwCreateZeroReadWriteBuffer(ci, sizeof(TreeStatus));

    nbb->start = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(int));
    nbb->count = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(int));
    nbb->sort = mwCreateZeroReadWriteBuffer(ci, 2 * st->nbody * sizeof(int)); /* CHECKME */
    nbb->child = mwCreateZeroReadWriteBuffer(ci, NSUB * (nNode + 1) * sizeof(int));

    nbb->debug = mwCreateZeroReadWriteBuffer(ci, sizeof(Debug));


    if (!nbb->masses || !nbb->treeStatus || !nbb->start || !nbb->count || !nbb->sort || !nbb->child)
        err = MW_CL_ERROR;

    if (err != CL_SUCCESS)
    {
        warn("Error creating NBody Buffers\n");
        releaseBuffers(nbb);
        return err;
    }

    return setInitialTreeStatus(ci, nbb);
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

static void setCLRequestFromFlags(CLRequest* clr, const NBodyFlags* nbf)
{
    clr->platform = nbf->platform;
    clr->devNum = nbf->devNum;
    clr->verbose = TRUE;
    clr->enableCheckpointing = FALSE;
}

cl_int runSystemCL(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    cl_int err = CL_SUCCESS;
    CLInfo ci;
    CLRequest clr;
    NBodyBuffers nbb;

    memset(&ci, 0, sizeof(ci));
    memset(&clr, 0, sizeof(clr));
    memset(&nbb, 0, sizeof(nbb));

    setCLRequestFromFlags(&clr, nbf);

    err = mwSetupCL(&ci, &clr);
    if (err != CL_SUCCESS)
        return err;

    if (!nbodyCheckDevCapabilities(&ci.di, ctx, st))
        return MW_CL_ERROR;

    err = loadKernels(&ci, ctx, st);
    if (err != CL_SUCCESS)
        goto fail;

    err = createKernels(&ci);
    if (err != CL_SUCCESS)
        goto fail;

    err = createBuffers(&nbb, &ci, st);
    if (err != CL_SUCCESS)
        goto fail;

    err = marshalBodies(&nbb, &ci, st, CL_TRUE);
    if (err != CL_SUCCESS)
        goto fail;

    err = setAllKernelArguments(&ci, &nbb);
    if (err != CL_SUCCESS)
        goto fail;

    err = nbodyMainLoop(&ci, ctx, st, &nbb);
    if (err != CL_SUCCESS)
        goto fail;



    err = marshalBodies(&nbb, &ci, st, CL_FALSE);
    if (err != CL_SUCCESS)
        goto fail;

    //printBodies(st->bodytab, st->nbody);

fail:
    debug(&ci, &nbb);

    mwDestroyCLInfo(&ci);
    releaseKernels();
    releaseBuffers(&nbb);

    return err;
}

