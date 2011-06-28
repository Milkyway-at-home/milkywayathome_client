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
    cl_mem treeControl;

    cl_mem start; /* TODO: We can reuse other buffers with this later to save memory */
    cl_mem count;
    cl_mem child;
    cl_mem sort;
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

static void printTreeStatus(const TreeStatus* tc)
{
    warn("TreeStatus = {\n"
         "  radius    = %.15f\n"
         "  bottom    = %d\n"
         "  maxDepth  = %d\n"
         "  errorCode = %d\n"
         "  blckCnt   = %u\n"
         "}\n",
         tc->radius,
         tc->bottom,
         tc->maxDepth,
         tc->errorCode,
         tc->blkCnt);
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
    err |= clSetKernelArg(kern, 20, sizeof(cl_mem), &nbb->treeControl);


    /* Set the step to be 0 here. Only the force calculation kernel
     * actually cares. We'll set it before then. */
    err |= clSetKernelArg(kern, 21, sizeof(cl_int), &step);

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
    //cl_uint warpSize = 32;
    cl_uint warpSize = 1;

    if (nNode < 1024 * di->maxCompUnits)
        nNode = 1024 * di->maxCompUnits;
    while ((nNode & (warpSize - 1)) != 0)
        nNode++;

    return nNode - 1;
}

static char* getCompileFlags(const NBodyCtx* ctx, const NBodyState* st, const DevInfo* di)
{
    char* buf;
    cl_uint warpSize = 1;

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

                 "-D MAXDEPTH=%u "
                 "-D THREADS1=%u "
                 "-D THREADS2=%u "
                 "-D THREADS3=%u "
                 "-D THREADS4=%u "
                 "-D THREADS5=%u "
                 "-D THREADS6=%u "

                 "-D TIMESTEP=%a "
                 "-D EPS2=%a "
                 "-D THETA=%a "

                 "-D NEWCRITERION=%d "
                 "-D SW93=%d "
                 "-D BH86=%d "
                 "-D EXACT=%d ",
                 DOUBLEPREC,

                 st->nbody,
                 findNNode(di, st->nbody),
                 warpSize,

                 512, /* THREADS1 */
                 288,
                 256,
                 512,
                 384,
                 512,
                 26,  /* MAXDEPTH */

                 ctx->timestep,
                 ctx->eps2,
                 ctx->theta,

                 ctx->criterion == NewCriterion,
                 ctx->criterion == SW93,
                 ctx->criterion == BH86,
                 ctx->criterion == Exact
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

static cl_int readTreeStatus(TreeStatus* tc, CLInfo* ci, NBodyBuffers* nbb)
{
    assert(tc);

    return clEnqueueReadBuffer(ci->queue,
                               nbb->treeControl,
                               CL_TRUE,
                               0, sizeof(*tc), tc,
                               0, NULL, NULL);
}

static cl_int stepSystemCL(CLInfo* ci, const NBodyCtx* ctx, NBodyState* st)
{
    cl_int err;
    size_t global[1] = { st->nbody };
    size_t* local = NULL;

    err = clSetKernelArg(kernels.forceCalculation, 21, sizeof(int), &st->step);
    if (err != CL_SUCCESS)
        return err;

    err = clEnqueueNDRangeKernel(ci->queue, kernels.boundingBox, 1,
                                 NULL, global, local,
                                 0, NULL, NULL);
    if (err != CL_SUCCESS)
        return err;

    st->tnow += ctx->timestep;
    st->step++;

    return clFinish(ci->queue);
}

static cl_int nbodyMainLoop(CLInfo* ci, const NBodyCtx* ctx, NBodyState* st)
{
    cl_int err = CL_SUCCESS;
    const real tstop = ctx->timeEvolve - ctx->timestep / 1024.0;

    st->step = -1;

    //while (err == CL_SUCCESS && st->tnow < tstop)
    {
        err = stepSystemCL(ci, ctx, st);
    }

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
    err |= clReleaseMemObject_quiet(nbb->treeControl);
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
                                nbb->treeControl,
                                CL_TRUE,
                                0, sizeof(TreeStatus), &iniTreeStatus,
                                0, NULL, NULL);
}

static cl_int createBuffers(NBodyBuffers* nbb, CLInfo* ci, NBodyState* st)
{
    cl_uint i;
    cl_uint nNode;
    cl_int err = CL_SUCCESS;

    for (i = 0; i < 3; ++i)
    {
        nbb->pos[i] = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(real));
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

    nbb->masses = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(real));
    nbb->treeControl = mwCreateZeroReadWriteBuffer(ci, sizeof(TreeStatus));

    nNode = findNNode(&ci->di, st->nbody);
    nbb->start = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(int));
    nbb->count = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(int));
    nbb->sort = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(int));
    nbb->child = mwCreateZeroReadWriteBuffer(ci, (NSUB * nNode + 1) * sizeof(int));


    if (!nbb->masses || !nbb->treeControl || !nbb->start)
        err = MW_CL_ERROR;

    if (err != CL_SUCCESS)
    {
        warn("Error creating NBody Buffers\n");
        releaseBuffers(nbb);
        return err;
    }

    return setInitialTreeStatus(ci, nbb);
}

static real* mapRealBuffer(CLInfo* ci, cl_mem mem, size_t nElement)
{
    return (real*) clEnqueueMapBuffer(ci->queue, mem, CL_TRUE, CL_MAP_READ, 0,
                                      nElement * sizeof(real),
                                      0, NULL, NULL, NULL);

}

static cl_int mapBodies(real* pos[3], real* vel[3], NBodyBuffers* nbb, CLInfo* ci, NBodyState* st)
{
    cl_uint i;

    for (i = 0; i < 3; ++i)
    {
        pos[i] = mapRealBuffer(ci, nbb->pos[i], (size_t) st->nbody);
        vel[i] = mapRealBuffer(ci, nbb->vel[i], (size_t) st->nbody);
        if (!pos[i] || !vel[i])
        {
            return MW_CL_ERROR;
        }
    }

    return CL_SUCCESS;
}

static cl_int unmapBodies(real* pos[3], real* vel[3], NBodyBuffers* nbb, CLInfo* ci)
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

    err = mapBodies(pos, vel, nbb, ci, st);
    if (err != CL_SUCCESS)
    {
        unmapBodies(pos, vel, nbb, ci);
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
        }
    }

    return unmapBodies(pos, vel, nbb, ci);
}

static cl_int printRealBuffer(CLInfo* ci, cl_mem mem, size_t n, const char* name)
{
    size_t i;
    real* p;

    p = mapRealBuffer(ci, mem, n);
    if (!p)
        return MW_CL_ERROR;

    for (i = 0; i < n; ++i)
    {
        warn("%s["ZU"] = %.15f\n", name, i, p[i]);
    }

    return clEnqueueUnmapMemObject(ci->queue, mem, p, 0, NULL, NULL);
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

    err = nbodyMainLoop(&ci, ctx, st);
    if (err != CL_SUCCESS)
        goto fail;

    {
        TreeStatus tc;
        memset(&tc, 0, sizeof(tc));
        err = readTreeStatus(&tc, &ci, &nbb);
        if (err != CL_SUCCESS)
            goto fail;
        printTreeStatus(&tc);
    }

    err = marshalBodies(&nbb, &ci, st, CL_FALSE);
    if (err != CL_SUCCESS)
        goto fail;

fail:
    mwDestroyCLInfo(&ci);
    releaseKernels();
    releaseBuffers(&nbb);

    return err;
}

