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
#include "nbody_curses.h"


typedef struct MW_ALIGN_TYPE_V(64)
{
    real radius;
    cl_int bottom;
    cl_int maxDepth;
    cl_uint blkCnt;

    cl_int errorCode;
    cl_int assertionLine;

    char _pad[64 - (sizeof(real) - 5 * sizeof(cl_int))];

    struct
    {
        real f[32];
        cl_int i[64];
    } debug;
} TreeStatus;


extern const unsigned char nbody_kernels_cl[];
extern const size_t nbody_kernels_cl_len;


static cl_ulong nbCalculateDepthLimitationFromCalculatedForceKernelLocalMemoryUsage(const DevInfo* di, const NBodyWorkSizes* ws, cl_bool useQuad)
{
    cl_ulong estMaxDepth;
    cl_ulong wgSize = ws->threads[5];
    cl_ulong warpPerWG = wgSize / di->warpSize;

    /* I don't trust the device parameters reporting anymore */
    cl_ulong localMemSize = (di->localMemSize > 0) ? di->localMemSize : 16384;

    /* Pieces which are not part of the "stack" */
    cl_ulong maxDepth = sizeof(cl_int);
    cl_ulong rootCritRadius = sizeof(real);
    cl_ulong allBlock = wgSize * sizeof(cl_int);

    cl_ulong ch = warpPerWG * sizeof(cl_int);
    cl_ulong nx = warpPerWG * sizeof(cl_int);
    cl_ulong ny = warpPerWG * sizeof(cl_int);
    cl_ulong nz = warpPerWG * sizeof(cl_int);
    cl_ulong nm = warpPerWG * sizeof(cl_int);

    cl_ulong constantPieces = maxDepth + rootCritRadius + allBlock + ch + nx + ny + nz + nm;

    /* Individual sizes of elements on the cell stack. */
    cl_ulong pos = sizeof(cl_int);
    cl_ulong node = sizeof(cl_int);
    cl_ulong dq = sizeof(real);
    cl_ulong quadPieces = 6 * sizeof(real);

    cl_ulong stackItemCount = pos + node + dq;

    if (useQuad)
    {
        stackItemCount += quadPieces;
    }

    /* We now have the size requirement as:
       d * warpPerWG * stackItemCount + constantPieces <= localMemSize
       Solve for d.
    */
    estMaxDepth = (localMemSize - constantPieces) / (warpPerWG * stackItemCount);

    return estMaxDepth - 1;  /* A bit extra will be used. Some kind of rounding up */
}

static cl_uint nbFindMaxDepthForDevice(const DevInfo* di, const NBodyWorkSizes* ws, cl_bool useQuad)
{
    cl_ulong d;

    d = nbCalculateDepthLimitationFromCalculatedForceKernelLocalMemoryUsage(di, ws, useQuad);

    return (cl_uint) d;
}


static void nbPrintNBodyWorkSizes(const NBodyWorkSizes* ws)
{
    mw_printf("\n"
              "Kernel launch sizes:\n"
              "  Bounding box kernel:  "ZU", "ZU"\n"
              "  Tree build kernel:    "ZU", "ZU"\n"
              "  Summarization kernel: "ZU", "ZU"\n"
              "  Sort kernel:          "ZU", "ZU"\n"
              "  Quadrupole kernel:    "ZU", "ZU"\n"
              "  Force kernel:         "ZU", "ZU"\n"
              "  Integration kernel:   "ZU", "ZU"\n"
              "\n",
              ws->global[0], ws->local[0],
              ws->global[1], ws->local[1],
              ws->global[2], ws->local[2],
              ws->global[3], ws->local[3],
              ws->global[4], ws->local[4],
              ws->global[5], ws->local[5],
              ws->global[6], ws->local[6]
        );
}

cl_bool nbSetWorkSizes(NBodyWorkSizes* ws, const DevInfo* di)
{
    cl_uint i;
    cl_uint blocks = di->maxCompUnits;

    for (i = 0; i < 8; ++i)
    {
        ws->global[i] = ws->threads[i] * ws->factors[i] * blocks;
        ws->local[i] = ws->threads[i];
    }

    return CL_FALSE;
}

/* CHECKME: May not apply on GT200? */
static cl_bool nbShouldForceLargeGroup(const DevInfo* di, const NBodyCtx* ctx)
{
    return !ctx->useQuad && mwIsNvidiaGPUDevice(di) && mwHasNvidiaCompilerFlags(di);
}

static const char* nbMaybeNvMaxRegCount(const DevInfo* di, const NBodyCtx* ctx)
{
    return nbShouldForceLargeGroup(di, ctx) ? "-cl-nv-maxrregcount=32 " : "";
}

/* Return CL_TRUE if some error */
cl_bool nbSetThreadCounts(NBodyWorkSizes* ws, const DevInfo* di, const NBodyCtx* ctx)
{
    /* Numbers need playing for float and different opening criteria */

    ws->factors[0] = 1;
    ws->factors[1] = 1;
    ws->factors[2] = 1;  /* Must be 1. All workitems must be resident */
    ws->factors[3] = 1;  /* Also must be 1 for the same reason */
    ws->factors[4] = 1;  /* Also must be 1 for the same reason */
    ws->factors[5] = 1;
    ws->factors[6] = 1;
    ws->factors[7] = 1;

    ws->threads[0] = 64;
    ws->threads[1] = 64;
    ws->threads[2] = 64;
    ws->threads[3] = 64;
    ws->threads[4] = 64;
    ws->threads[5] = 64;
    ws->threads[6] = 64;
    ws->threads[7] = 64;

    if (di->devType == CL_DEVICE_TYPE_CPU)
    {
        ws->threads[0] = 1;
        ws->threads[1] = 1;
        ws->threads[2] = 1;
        ws->threads[3] = 1;
        ws->threads[4] = 1;
        ws->threads[5] = 1;
        ws->threads[6] = 1;
        ws->threads[7] = 1;
    }
    else if (mwComputeCapabilityIs(di, 1, 3))
    {
        ws->threads[0] = 256;
        ws->threads[1] = 288;
        ws->threads[2] = 256;
        ws->threads[3] = 512;
        ws->threads[4] = 256;
        ws->threads[5] = 384;
        ws->threads[6] = 512;
        ws->threads[7] = 448;
    }
    else if (mwMinComputeCapabilityCheck(di, 2, 0))
    {
        ws->factors[0] = 1;
        ws->factors[1] = 1;
        ws->factors[2] = 1;
        ws->factors[3] = 1;
        ws->factors[4] = 1;
        ws->factors[5] = 1;
        ws->factors[6] = 4;
        ws->factors[7] = 4;

        ws->threads[0] = 1024;
        ws->threads[1] = 1024;
        ws->threads[2] = 1024;
        ws->threads[3] = 1024;
        ws->threads[4] = 1024;

        /* It's faster to restrain the used number of registers and
         * get a larger workgroup size, but when using quadrupole
         * moments this gives a very small constraining maximum depth */
        ws->threads[5] = nbShouldForceLargeGroup(di, ctx) ? 1024 : 512;

        ws->threads[6] = 1024;
        ws->threads[7] = 1024;
    }
    else
    {
        ws->factors[0] = 1;
        ws->factors[1] = 1;
        ws->factors[2] = 1;
        ws->factors[3] = 1;
        ws->factors[4] = 1;
        ws->factors[5] = 1;
        ws->factors[6] = 2;
        ws->factors[7] = 1;

        ws->threads[0] = 256;
        ws->threads[1] = 256;
        ws->threads[2] = 256;
        ws->threads[3] = 256;
        ws->threads[4] = 256;
        ws->threads[5] = 256;
        ws->threads[6] = 256;
        ws->threads[7] = 256;
    }

    return CL_FALSE;
}


static void* mapBuffer(CLInfo* ci, cl_mem mem, cl_map_flags flags, size_t size)
{
    return clEnqueueMapBuffer(ci->queue, mem, CL_TRUE, flags, 0,
                              size,
                              0, NULL, NULL, NULL);
}

static void nbPrintDebug(const TreeStatus* ts)
{
    int i;
    for (i = 0; i < 32; ++i)
    {
        mw_printf("Debug.int[%d] = %d\n", i, ts->debug.i[i]);
    }

    for (i = 0; i < 32; ++i)
    {
        mw_printf("Debug.float[%d] = %.15f\n", i, ts->debug.f[i]);
    }
}

static void nbPrintTreeStatus(const TreeStatus* ts)
{
    mw_printf("TreeStatus = {\n"
              "  radius        = %.15f\n"
              "  bottom        = %d\n"
              "  maxDepth      = %d\n"
              "  blckCnt       = %u\n"
              "  errorCode     = %d\n"
              "  assertionLine = %d\n"
              "}\n",
              ts->radius,
              ts->bottom,
              ts->maxDepth,
              ts->blkCnt,
              ts->errorCode,
              ts->assertionLine
        );
}

/* Set arguments (idx, idx + 1, idx + 2) to the buffers in mem[3] */
static cl_int nbSetMemArrayArgs(cl_kernel kern, cl_mem mem[3], cl_uint idx)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;

    for (i = 0; i < 3; ++i)
    {
        err |= clSetKernelArg(kern, idx + i, sizeof(cl_mem), &mem[i]);
    }

    return err;
}

static cl_int nbSetKernelArguments(cl_kernel kern, NBodyBuffers* nbb, cl_bool exact)
{
    cl_int err = CL_SUCCESS;
    cl_int zeroVal = 0;

    /* Just any valid buffer that we'll use for arguments we don't need */
    cl_mem* anything = &nbb->masses;

    if (!exact)
    {
        err |= nbSetMemArrayArgs(kern, nbb->pos, 0);
        err |= nbSetMemArrayArgs(kern, nbb->vel, 3);
        err |= nbSetMemArrayArgs(kern, nbb->acc, 6);

        err |= clSetKernelArg(kern, 9, sizeof(cl_mem), &nbb->masses);

        err |= nbSetMemArrayArgs(kern, nbb->max, 10);
        err |= nbSetMemArrayArgs(kern, nbb->min, 13);

        err |= clSetKernelArg(kern, 16, sizeof(cl_mem), &nbb->start);
        err |= clSetKernelArg(kern, 17, sizeof(cl_mem), &nbb->count);
        err |= clSetKernelArg(kern, 18, sizeof(cl_mem), &nbb->child);
        err |= clSetKernelArg(kern, 19, sizeof(cl_mem), &nbb->sort);
        err |= clSetKernelArg(kern, 20, sizeof(cl_mem), nbb->critRadii ? &nbb->critRadii : anything);


        if (nbb->quad.xx) /* If we're using quadrupole moments */
        {
            err |= clSetKernelArg(kern, 21, sizeof(cl_mem), &nbb->quad.xx);
            err |= clSetKernelArg(kern, 22, sizeof(cl_mem), &nbb->quad.xy);
            err |= clSetKernelArg(kern, 23, sizeof(cl_mem), &nbb->quad.xz);

            err |= clSetKernelArg(kern, 24, sizeof(cl_mem), &nbb->quad.yy);
            err |= clSetKernelArg(kern, 25, sizeof(cl_mem), &nbb->quad.yz);

            err |= clSetKernelArg(kern, 26, sizeof(cl_mem), &nbb->quad.zz);
        }
        else
        {
            err |= clSetKernelArg(kern, 21, sizeof(cl_mem), anything);
            err |= clSetKernelArg(kern, 22, sizeof(cl_mem), anything);
            err |= clSetKernelArg(kern, 23, sizeof(cl_mem), anything);

            err |= clSetKernelArg(kern, 24, sizeof(cl_mem), anything);
            err |= clSetKernelArg(kern, 25, sizeof(cl_mem), anything);

            err |= clSetKernelArg(kern, 26, sizeof(cl_mem), anything);
        }

        err |= clSetKernelArg(kern, 27, sizeof(cl_mem), &nbb->treeStatus);
        err |= clSetKernelArg(kern, 28, sizeof(cl_int), &zeroVal);
        err |= clSetKernelArg(kern, 29, sizeof(cl_int), &zeroVal);
    }
    else
    {
        cl_int i;

        err |= nbSetMemArrayArgs(kern, nbb->pos, 0);
        err |= nbSetMemArrayArgs(kern, nbb->vel, 3);
        err |= nbSetMemArrayArgs(kern, nbb->acc, 6);

        err |= clSetKernelArg(kern, 9, sizeof(cl_mem), &nbb->masses);

        for (i = 10; i < 27; ++i)
        {
            err |= clSetKernelArg(kern, i, sizeof(cl_mem), anything);
        }

        err |= clSetKernelArg(kern, 27, sizeof(cl_mem), &nbb->treeStatus);
        err |= clSetKernelArg(kern, 28, sizeof(cl_int), &zeroVal);
        err |= clSetKernelArg(kern, 29, sizeof(cl_int), &zeroVal);
    }

    return err;
}

cl_int nbSetAllKernelArguments(NBodyState* st)
{
    cl_int err = CL_SUCCESS;
    NBodyKernels* k = st->kernels;
    cl_bool exact = st->usesExact;

    if (!exact)
    {
        err |= nbSetKernelArguments(k->boundingBox, st->nbb, exact);
        err |= nbSetKernelArguments(k->buildTree, st->nbb, exact);
        err |= nbSetKernelArguments(k->summarization, st->nbb, exact);
        err |= nbSetKernelArguments(k->sort, st->nbb, exact);
        err |= nbSetKernelArguments(k->quadMoments, st->nbb, exact);
        err |= nbSetKernelArguments(k->forceCalculation, st->nbb, exact);
        err |= nbSetKernelArguments(k->integration, st->nbb, exact);
    }
    else
    {
        err |= nbSetKernelArguments(k->forceCalculation_Exact, st->nbb, exact);
        err |= nbSetKernelArguments(k->integration, st->nbb, exact);
    }

    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error setting kernel arguments");
    }

    return err;
}

static cl_int clReleaseKernel_quiet(cl_kernel kern)
{
    return kern ? clReleaseKernel(kern) : CL_SUCCESS;
}

cl_int nbReleaseKernels(NBodyState* st)
{
    cl_int err = CL_SUCCESS;
    NBodyKernels* kernels = st->kernels;

    err |= clReleaseKernel_quiet(kernels->boundingBox);
    err |= clReleaseKernel_quiet(kernels->buildTree);
    err |= clReleaseKernel_quiet(kernels->summarization);
    err |= clReleaseKernel_quiet(kernels->quadMoments);
    err |= clReleaseKernel_quiet(kernels->sort);
    err |= clReleaseKernel_quiet(kernels->forceCalculation);
    err |= clReleaseKernel_quiet(kernels->integration);

    if (err != CL_SUCCESS)
        mwPerrorCL(err, "Error releasing kernels");

    return err;
}

static cl_uint nbFindNNode(const DevInfo* di, cl_int nbody)
{
    cl_uint nNode = 2 * nbody;

    if (nNode < 1024 * di->maxCompUnits)
        nNode = 1024 * di->maxCompUnits;
    while ((nNode & (di->warpSize - 1)) != 0)
        ++nNode;

    return nNode - 1;
}

static char* nbGetCompileFlags(const NBodyCtx* ctx, const NBodyState* st, const DevInfo* di)
{
    char* buf;
    const NBodyWorkSizes* ws = st->workSizes;
    const Potential* p = &ctx->pot;

    /* Put a space between the -D. if the thing begins with D, there's
     * an Apple OpenCL compiler bug where the D will be
     * stripped. -DDOUBLEPREC=1 will actually define OUBLEPREC */
    if (asprintf(&buf,
                 "-D DOUBLEPREC=%d "
               #if !DOUBLEPREC
                 "-cl-single-precision-constant "
               #endif

                 "-cl-mad-enable "

               #ifndef NDEBUG
                 "-D DEBUG=1 "
               #else
                 "-D DEBUG=0 "
               #endif

                 "-D NBODY=%d "
                 "-D EFFNBODY=%d "
                 "-D NNODE=%u "
                 "-D WARPSIZE=%u "

                 "-D NOSORT=%d "

                 "-D THREADS1="ZU" "
                 "-D THREADS2="ZU" "
                 "-D THREADS3="ZU" "
                 "-D THREADS4="ZU" "
                 "-D THREADS5="ZU" "
                 "-D THREADS6="ZU" "
                 "-D THREADS7="ZU" "
                 "-D THREADS8="ZU" "

                 "-D MAXDEPTH=%u "

                 "-D TIMESTEP=%a "
                 "-D EPS2=%a "
                 "-D THETA=%a "
                 "-D USE_QUAD=%d "

                 "-D NEWCRITERION=%d "
                 "-D SW93=%d "
                 "-D BH86=%d "
                 "-D EXACT=%d "

                 /* Potential */
                 "-D USE_EXTERNAL_POTENTIAL=%d "
                 "-D MIYAMOTO_NAGAI_DISK=%d "
                 "-D EXPONENTIAL_DISK=%d "
                 "-D LOG_HALO=%d "
                 "-D NFW_HALO=%d "
                 "-D TRIAXIAL_HALO=%d "

                 /* Spherical constants */
                 "-D SPHERICAL_MASS=%a "
                 "-D SPHERICAL_SCALE=%a "

                 /* Disk constants */
                 "-D DISK_MASS=%a "
                 "-D DISK_SCALE_LENGTH=%a "
                 "-D DISK_SCALE_HEIGHT=%a "

                 /* Halo constants */
                 "-D HALO_VHALO=%a "
                 "-D HALO_SCALE_LENGTH=%a "
                 "-D HALO_FLATTEN_Z=%a "
                 "-D HALO_FLATTEN_Y=%a "
                 "-D HALO_FLATTEN_X=%a "
                 "-D HALO_TRIAX_ANGLE=%a "
                 "-D HALO_C1=%a "
                 "-D HALO_C2=%a "
                 "-D HALO_C3=%a "

                 "%s "
                 "%s "
                 "-D HAVE_INLINE_PTX=%d ",
                 DOUBLEPREC,

                 st->nbody,
                 st->effNBody,
                 nbFindNNode(di, st->nbody),
                 di->warpSize,

                 (di->devType == CL_DEVICE_TYPE_CPU),

                 ws->threads[0],
                 ws->threads[1],
                 ws->threads[2],
                 ws->threads[3],
                 ws->threads[4],
                 ws->threads[5],
                 ws->threads[6],
                 ws->threads[7],
                 nbFindMaxDepthForDevice(di, st->workSizes, ctx->useQuad),

                 ctx->timestep,
                 ctx->eps2,
                 ctx->theta,
                 ctx->useQuad,

                 /* Set criterion */
                 ctx->criterion == NewCriterion,
                 ctx->criterion == SW93,
                 ctx->criterion == BH86,
                 ctx->criterion == Exact,


                 /* Set potential */
                 ctx->potentialType == EXTERNAL_POTENTIAL_DEFAULT,

                 p->disk.type == MiyamotoNagaiDisk,
                 p->disk.type == ExponentialDisk,
                 p->halo.type == LogarithmicHalo,
                 p->halo.type == NFWHalo,
                 p->halo.type == TriaxialHalo,

                 /* Set potential constants */
                 /* Spherical constants */
                 p->sphere[0].mass,
                 p->sphere[0].scale,

                 /* Disk constants */
                 p->disk.mass,
                 p->disk.scaleLength,
                 p->disk.scaleHeight,

                 /* Halo constants */
                 p->halo.vhalo,
                 p->halo.scaleLength,
                 p->halo.flattenZ,
                 p->halo.flattenY,
                 p->halo.flattenX,
                 p->halo.triaxAngle,
                 p->halo.c1,
                 p->halo.c2,
                 p->halo.c3,

                 /* Misc. other stuff */
                 mwHasNvidiaCompilerFlags(di) ? "-cl-nv-verbose" : "",
                 nbMaybeNvMaxRegCount(di, ctx),
                 mwNvidiaInlinePTXAvailable(st->ci->plat)
            ) < 1)
    {
        mw_printf("Error getting compile flags\n");
        return NULL;
    }

    return buf;
}

static cl_bool nbCreateKernels(cl_program program, NBodyKernels* kernels)
{
    kernels->boundingBox = mwCreateKernel(program, "boundingBox");
    kernels->buildTree = mwCreateKernel(program, "buildTree");
    kernels->summarization = mwCreateKernel(program, "summarization");
    kernels->quadMoments = mwCreateKernel(program, "quadMoments");
    kernels->sort = mwCreateKernel(program, "sort");
    kernels->forceCalculation = mwCreateKernel(program, "forceCalculation");
    kernels->integration = mwCreateKernel(program, "integration");
    kernels->forceCalculation_Exact = mwCreateKernel(program, "forceCalculation_Exact");

    return (   kernels->boundingBox
            && kernels->buildTree
            && kernels->summarization
            && kernels->quadMoments
            && kernels->sort
            && kernels->forceCalculation
            && kernels->integration
            && kernels->forceCalculation_Exact);
}

cl_bool nbLoadKernels(const NBodyCtx* ctx, NBodyState* st)
{
    CLInfo* ci = st->ci;
    cl_bool failed;
    char* compileFlags = NULL;
    cl_program program;
    const char* src = (const char*) nbody_kernels_cl;
    size_t srcLen = nbody_kernels_cl_len;

    compileFlags = nbGetCompileFlags(ctx, st, &ci->di);
    assert(compileFlags);

    program = mwCreateProgramFromSrc(ci, 1, &src, &srcLen, compileFlags);
    free(compileFlags);
    if (!program)
    {
        return CL_TRUE;
    }

    failed = nbCreateKernels(program, st->kernels);
    clReleaseProgram(program);

    return failed;
}

/* Return CL_FALSE if device isn't capable of running this */
cl_bool nbCheckDevCapabilities(const DevInfo* di, const NBodyCtx* ctx, cl_uint nbody)
{
    cl_ulong nNode = (cl_ulong) nbFindNNode(di, nbody) + 1;
    cl_ulong maxNodes = di->maxMemAlloc / (NSUB * sizeof(cl_int));

    (void) ctx;

    if (di->devType != CL_DEVICE_TYPE_GPU)
    {
        mw_printf("Device is not a GPU.\n");
        return CL_FALSE;
    }

    if (!mwIsNvidiaGPUDevice(di) && !mwIsAMDGPUDevice(di))
    {
        /* There is reliance on implementation details for Nvidia and
         * AMD GPUs. If some other kind of GPU decides to exist, it
         * would need to be tested.*/
        mw_printf("Only Nvidia and AMD GPUs are supported\n");
        return CL_FALSE;
    }

    if (DOUBLEPREC && !mwSupportsDoubles(di))
    {
        mw_printf("Device does not have usable double precision extension\n");
        return CL_FALSE;
    }

    if (   !strstr(di->exts, "cl_khr_global_int32_base_atomics")
        || !strstr(di->exts, "cl_khr_global_int32_extended_atomics")
        || !strstr(di->exts, "cl_khr_local_int32_base_atomics"))
    {
        mw_printf("Device lacks necessary atomics extensions\n");
        return CL_FALSE;
    }


    if (nNode > maxNodes)
    {
        mw_printf("Simulation of %u bodies requires "LLU" nodes, "
                  "however maximum allocation size only allows for "LLU"\n",
                  nbody,
                  nNode,
                  maxNodes
            );
        return CL_FALSE;
    }

    return CL_TRUE;
}

static cl_int nbReadTreeStatus(TreeStatus* tc, CLInfo* ci, NBodyBuffers* nbb)
{
    cl_int err;
    assert(tc);

    err = clEnqueueReadBuffer(ci->queue,
                              nbb->treeStatus,
                              CL_TRUE,
                              0, sizeof(*tc), tc,
                              0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error reading tree status");
    }

    return err;
}

static cl_int printBuffer(CLInfo* ci, cl_mem mem, size_t n, const char* name, int type)
{
    size_t i;
    void* p;

    p = mapBuffer(ci, mem, CL_MAP_READ, n * (type == 0 ? sizeof(real) : sizeof(int)));
    if (!p)
    {
        mw_printf("Fail to map buffer for printing\n");
        return MW_CL_ERROR;
    }

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
        const int* ip = (const int*) p;
        for (i = 0; i < n; ++i)
        {
            mw_printf("%s["ZU"] = %d\n", name, i, ip[i]);
        }
    }

    return clEnqueueUnmapMemObject(ci->queue, mem, p, 0, NULL, NULL);
}

static void stdDebugPrint(NBodyState* st, cl_bool children, cl_bool tree, cl_bool quads)
{
    cl_int err;
    CLInfo* ci = st->ci;
    NBodyBuffers* nbb = st->nbb;
    cl_uint nNode = nbFindNNode(&ci->di, st->effNBody);

    if (children)
    {
        mw_printf("--------------------------------------------------------------------------------\n");

        mw_printf("BEGIN CHILD\n");
        printBuffer(ci, nbb->child, NSUB * (nNode + 1), "child", 1);
        mw_printf("BEGIN START\n");
        printBuffer(ci, nbb->start, nNode, "start", 1);

        mw_printf("BEGIN MASS\n");
        printBuffer(ci, nbb->masses, nNode + 1, "mass", 0);

        mw_printf("BEGIN POSX\n");
        printBuffer(ci, nbb->pos[0], nNode + 1, "posX", 0);
        mw_printf("BEGIN POSY\n");
        printBuffer(ci, nbb->pos[1], nNode + 1, "posY", 0);
        mw_printf("BEGIN POSY\n");
        printBuffer(ci, nbb->pos[2], nNode + 1, "posZ", 0);

        mw_printf("BEGIN VELX\n");
        printBuffer(ci, nbb->vel[0], st->effNBody, "velX", 0);
        mw_printf("BEGIN VELY\n");
        printBuffer(ci, nbb->vel[1], st->effNBody, "velY", 0);
        mw_printf("BEGIN VELY\n");
        printBuffer(ci, nbb->vel[2], st->effNBody, "velZ", 0);

        mw_printf("BEGIN ACCX\n");
        printBuffer(ci, nbb->acc[0], st->effNBody, "accX", 0);
        mw_printf("BEGIN ACCY\n");
        printBuffer(ci, nbb->acc[1], st->effNBody, "accY", 0);
        mw_printf("BEGIN ACCY\n");
        printBuffer(ci, nbb->acc[2], st->effNBody, "accZ", 0);
    }

    if (quads)
    {
        mw_printf("BEGIN QUAD.XX\n");
        printBuffer(ci, nbb->quad.xx, nNode + 1, "quad.xx", 0);
        mw_printf("BEGIN QUAD.XY\n");
        printBuffer(ci, nbb->quad.xy, nNode + 1, "quad.xy", 0);
        mw_printf("BEGIN QUAD.XZ\n");
        printBuffer(ci, nbb->quad.xz, nNode + 1, "quad.xz", 0);

        mw_printf("BEGIN QUAD.YY\n");
        printBuffer(ci, nbb->quad.yy, nNode + 1, "quad.yy", 0);
        mw_printf("BEGIN QUAD.YZ\n");
        printBuffer(ci, nbb->quad.yz, nNode + 1, "quad.yz", 0);

        mw_printf("BEGIN QUAD.ZZ\n");
        printBuffer(ci, nbb->quad.zz, nNode + 1, "quad.zz", 0);
    }


    if (tree)
    {
        TreeStatus ts;
        memset(&ts, 0, sizeof(ts));
        err = nbReadTreeStatus(&ts, ci, nbb);
        if (err != CL_SUCCESS)
        {
            mwPerrorCL(err, "Reading tree status failed\n");
        }
        else
        {
            nbPrintTreeStatus(&ts);
            nbPrintDebug(&ts);
        }
    }

    mw_printf("--------------------------------------------------------------------------------\n");
}

static NBodyStatus nbKernelErrorToNBodyStatus(NBodyKernelError x)
{
    switch (x)
    {
        case NBODY_KERNEL_OK:
            return NBODY_SUCCESS;
        case NBODY_KERNEL_CELL_OVERFLOW:
            return NBODY_CELL_OVERFLOW_ERROR;
        case NBODY_KERNEL_TREE_INCEST:
            return NBODY_TREE_INCEST_FATAL; /* Somewhat inaccurate but shouldn't happen  */
        case NBODY_KERNEL_TREE_STRUCTURE_ERROR:
            return NBODY_TREE_STRUCTURE_ERROR;
        case NBODY_KERNEL_ERROR_OTHER:
            return NBODY_ERROR;
        default:
            return NBODY_ERROR;
    }
}


/* Check the error code */
static NBodyStatus nbCheckKernelErrorCode(const NBodyCtx* ctx, NBodyState* st)
{
    cl_int err;
    TreeStatus ts;
    CLInfo* ci = st->ci;
    NBodyBuffers* nbb = st->nbb;

    err = nbReadTreeStatus(&ts, ci, nbb);
    if (mw_unlikely(err != CL_SUCCESS))
    {
        return NBODY_CL_ERROR;
    }

    if (mw_unlikely(ts.assertionLine >= 0))
    {
        mw_printf("Kernel assertion failed: line %d\n", ts.assertionLine);
        return NBODY_ASSERTION_FAILURE;
    }

    if (mw_unlikely(ts.errorCode != 0))
    {
        /* Incest is special because we cagn choose to ignore it */
        if (ts.errorCode == NBODY_KERNEL_TREE_INCEST)
        {
            nbReportTreeIncest(ctx, st);
            return ctx->allowIncest ? NBODY_TREE_INCEST_NONFATAL : NBODY_TREE_INCEST_FATAL;
        }
        else
        {
            mw_printf("Kernel reported error: %d ", ts.errorCode);

            if (ts.errorCode > 0)
            {
                mw_printf("(%s (%u))\n",
                          showNBodyKernelError(ts.errorCode),
                          nbFindMaxDepthForDevice(&ci->di, st->workSizes, st->usesQuad));
                return NBODY_MAX_DEPTH_ERROR;
            }
            else
            {
                mw_printf("(%s)\n", showNBodyKernelError(ts.errorCode));
                return nbKernelErrorToNBodyStatus(ts.errorCode);
            }
        }
    }

    return NBODY_SUCCESS;
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

static void nbReportProgressWithTimings(const NBodyCtx* ctx, const NBodyState* st)
{
    NBodyWorkSizes* ws = st->workSizes;

    mw_mvprintw(0, 0,
                "Step %d (%f%%):\n"
                "  boundingBox:      %15f ms\n"
                "  buildTree:        %15f ms%15f ms\n"
                "  summarization:    %15f ms\n"
                "  sort:             %15f ms\n"
                "  quad moments:     %15f ms\n"
                "  forceCalculation: %15f ms%15f ms\n"
                "  integration:      %15f ms\n"
                "\n",
                st->step,
                100.0 * (double) st->step / (double) ctx->nStep,
                ws->timings[0],
                ws->timings[1], ws->chunkTimings[1],
                ws->timings[2],
                ws->timings[3],
                ws->timings[4],
                ws->timings[5], ws->chunkTimings[5],
                ws->timings[6]
        );
    mw_refresh();
}

/* Run kernels used only by tree versions. */
static cl_int nbExecuteTreeConstruction(NBodyState* st)
{
    cl_int err;
    size_t chunk;
    size_t nChunk;
    cl_int upperBound;
    size_t offset[1];
    cl_event boxEv, sumEv, sortEv, quadEv;
    CLInfo* ci = st->ci;
    NBodyWorkSizes* ws = st->workSizes;
    NBodyKernels* kernels = st->kernels;
    cl_int effNBody = st->effNBody;

    err = clEnqueueNDRangeKernel(ci->queue, kernels->boundingBox, 1,
                                 NULL, &ws->global[0], &ws->local[0],
                                 0, NULL, &boxEv);
    if (err != CL_SUCCESS)
        return err;

    nChunk     = st->ignoreResponsive ?        1 : mwDivRoundup((size_t) effNBody, ws->global[1]);
    upperBound = st->ignoreResponsive ? effNBody : (cl_int) ws->global[1];
    for (chunk = 0, offset[0] = 0; chunk < nChunk; ++chunk, offset[0] += ws->global[1])
    {
        cl_event ev;

        if (upperBound > effNBody)
            upperBound = effNBody;

        err = clSetKernelArg(kernels->buildTree, 28, sizeof(cl_int), &upperBound);
        if (err != CL_SUCCESS)
            return err;

        err = clEnqueueNDRangeKernel(ci->queue, kernels->buildTree, 1,
                                     offset, &ws->global[1], &ws->local[1],
                                     0, NULL, &ev);
        if (err != CL_SUCCESS)
            return err;

        upperBound += (cl_int) ws->global[1];
        ws->timings[1] += waitReleaseEventWithTime(ev);
    }


    err = clEnqueueNDRangeKernel(ci->queue, kernels->summarization, 1,
                                 NULL, &ws->global[2], &ws->local[2],
                                 0, NULL, &sumEv);
    if (err != CL_SUCCESS)
        return err;

    /* FIXME: This does not work unless ALL of the threads are
     * launched at once. This may be bad when we need
     * responsiveness. This also means it will always hang with
     * CPUs. It seems to be fast enough though in every case I've
     * tried. */
    err = clEnqueueNDRangeKernel(ci->queue, kernels->sort, 1,
                                 NULL, &ws->global[3], &ws->local[3],
                                 0, NULL, &sortEv);
    if (err != CL_SUCCESS)
        return err;

    if (st->usesQuad)
    {
        err = clEnqueueNDRangeKernel(ci->queue, kernels->quadMoments, 1,
                                     NULL, &ws->global[4], &ws->local[4],
                                     0, NULL, &quadEv);
        if (err != CL_SUCCESS)
            return err;
    }

    ws->timings[0] += mwReleaseEventWithTiming(boxEv);
    ws->chunkTimings[1] = ws->timings[1] / (double) nChunk;
    ws->timings[2] += waitReleaseEventWithTime(sumEv);
    ws->timings[3] += waitReleaseEventWithTime(sortEv);
    if (st->usesQuad)
    {
        ws->timings[4] += waitReleaseEventWithTime(quadEv);
    }

    return CL_SUCCESS;
}

/* Run force calculation and integration kernels */
static cl_int nbExecuteForceKernels(NBodyState* st, cl_bool updateState)
{
    cl_int err;
    size_t chunk;
    size_t nChunk;
    cl_int upperBound;
    size_t global[1];
    size_t local[1];
    size_t offset[1];
    cl_event integrateEv;
    cl_kernel forceKern;
    CLInfo* ci = st->ci;
    NBodyKernels* kernels = st->kernels;
    NBodyWorkSizes* ws = st->workSizes;
    cl_int effNBody = st->effNBody;


    if (st->usesExact)
    {
        forceKern = kernels->forceCalculation_Exact;
        global[0] = ws->global[7];
        local[0] = ws->local[7];
    }
    else
    {
        forceKern = kernels->forceCalculation;
        global[0] = ws->global[5];
        local[0] = ws->local[5];
    }

    nChunk = st->ignoreResponsive ? 1 : mwDivRoundup((size_t) effNBody, global[0]);
    upperBound = st->ignoreResponsive ? effNBody : (cl_int) global[0];
    for (chunk = 0, offset[0] = 0; chunk < nChunk; ++chunk, offset[0] += global[0])
    {
        cl_event ev;

        upperBound = (upperBound > effNBody) ? effNBody : upperBound;

        err = clSetKernelArg(forceKern, 28, sizeof(cl_int), &upperBound);
        if (err != CL_SUCCESS)
            return err;

        err = clEnqueueNDRangeKernel(ci->queue, forceKern, 1,
                                     offset, global, local,
                                     0, NULL, &ev);
        if (err != CL_SUCCESS)
            return err;

        upperBound += (cl_int) global[0];
        ws->timings[5] += waitReleaseEventWithTime(ev);
    }

    if (mw_likely(updateState))
    {
        err = clEnqueueNDRangeKernel(ci->queue, kernels->integration, 1,
                                     NULL, &ws->global[6], &ws->local[6],
                                     0, NULL, &integrateEv);
        if (err != CL_SUCCESS)
            return err;
    }


    ws->chunkTimings[5] = ws->timings[5] / (double) nChunk;
    if (mw_likely(updateState))
    {
        ws->timings[6] += waitReleaseEventWithTime(integrateEv);
    }

    return CL_SUCCESS;
}

NBodyStatus nbStepSystemCL(const NBodyCtx* ctx, NBodyState* st)
{
    cl_int err;
    cl_uint i;
    NBodyWorkSizes* ws = st->workSizes;

    memset(ws->timings, 0, sizeof(ws->timings));

    if (!st->usesExact)
    {
        err = nbExecuteTreeConstruction(st);
        if (err != CL_SUCCESS)
        {
            mwPerrorCL(err, "Error executing tree construction kernels");
            return NBODY_CL_ERROR;
        }
    }

    err = nbExecuteForceKernels(st, CL_TRUE);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error executing force kernels");
        return NBODY_CL_ERROR;
    }

    if (st->reportProgress)
    {
        nbReportProgressWithTimings(ctx, st);
    }

    for (i = 0; i < 7; ++i) /* Add timings to running totals */
    {
        ws->kernelTimings[i] += ws->timings[i];
    }

    return NBODY_SUCCESS;
}

/* We need to run a fake step to get the initial accelerations without
 * touching the positons/velocities */
static cl_int nbRunPreStep(NBodyState* st)
{
    static const cl_int trueVal = TRUE;    /* Need an lvalue */
    static const cl_int falseVal = FALSE;
    cl_kernel kernel = st->usesExact ? st->kernels->forceCalculation_Exact : st->kernels->forceCalculation;
    cl_int err;

    /* Only calculate accelerations*/
    err = clSetKernelArg(kernel, 29, sizeof(cl_int), &falseVal);
    if (err != CL_SUCCESS)
        return err;

    if (!st->usesExact)
    {
        err = nbExecuteTreeConstruction(st);
        if (err != CL_SUCCESS)
            return err;
    }

    err = nbExecuteForceKernels(st, CL_FALSE);
    if (err != CL_SUCCESS)
        return err;

    /* All later steps will be real timesteps */
    return clSetKernelArg(kernel, 29, sizeof(cl_int), &trueVal);
}

static NBodyStatus nbMainLoopCL(const NBodyCtx* ctx, NBodyState* st)
{
    NBodyStatus rc = NBODY_SUCCESS;
    cl_int err;

    err = nbRunPreStep(st);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error running pre step");
        return NBODY_CL_ERROR;
    }

    while (st->step < ctx->nStep)
    {
        st->dirty = TRUE;

        rc = nbCheckKernelErrorCode(ctx, st);
        if (nbStatusIsFatal(rc))
        {
            return rc;
        }

        rc = nbStepSystemCL(ctx, st);
        if (nbStatusIsFatal(rc))
        {
            return rc;
        }

        st->step++;
    }

    return rc;
}

/* This is dumb and errors if mem isn't set */
static cl_int clReleaseMemObject_quiet(cl_mem mem)
{
    return mem ? clReleaseMemObject(mem) : CL_SUCCESS;
}

static cl_int _nbReleaseBuffers(NBodyBuffers* nbb)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;

    if (!nbb)
    {
        return CL_SUCCESS;
    }

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


    err |= clReleaseMemObject_quiet(nbb->quad.xx);
    err |= clReleaseMemObject_quiet(nbb->quad.xy);
    err |= clReleaseMemObject_quiet(nbb->quad.xz);

    err |= clReleaseMemObject_quiet(nbb->quad.yy);
    err |= clReleaseMemObject_quiet(nbb->quad.yz);

    err |= clReleaseMemObject_quiet(nbb->quad.zz);

    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error releasing buffers");
    }

    return err;
}

cl_int nbReleaseBuffers(NBodyState* st)
{
    return _nbReleaseBuffers(st->nbb);
}

cl_int nbSetInitialTreeStatus(NBodyState* st)
{
    cl_int err;
    TreeStatus iniTreeStatus;
    CLInfo* ci = st->ci;
    NBodyBuffers* nbb = st->nbb;

    memset(&iniTreeStatus, 0, sizeof(iniTreeStatus));

    iniTreeStatus.radius = 0;
    iniTreeStatus.bottom = 0;
    iniTreeStatus.maxDepth = 1;
    iniTreeStatus.assertionLine = -1;
    iniTreeStatus.errorCode = 0;
    iniTreeStatus.blkCnt = 0;

    err = clEnqueueWriteBuffer(ci->queue,
                               nbb->treeStatus,
                               CL_TRUE,
                               0, sizeof(TreeStatus), &iniTreeStatus,
                               0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error writing initial tree status");
    }

    return err;
}

static cl_uint nbFindInc(cl_uint warpSize, cl_uint nbody)
{
    return (nbody + warpSize - 1) & (-warpSize);
}

/* In some cases to avoid conditionally barriering we want to round up to nearest workgroup size.
   Find the least common multiple necessary for kernels that need to avoid the issue
 */
cl_int nbFindEffectiveNBody(const NBodyWorkSizes* workSizes, cl_bool exact, cl_int nbody)
{
    if (exact)
    {
        /* Exact force kernel needs this */
        return mwNextMultiple((cl_int) workSizes->local[7], nbody);
    }
    else
    {
        /* Maybe tree construction will need this later */
        return nbody;
    }
}

cl_int nbCreateBuffers(const NBodyCtx* ctx, NBodyState* st)
{
    cl_uint i;
    CLInfo* ci = st->ci;
    NBodyBuffers* nbb = st->nbb;
    size_t massSize;
    cl_uint nNode = nbFindNNode(&ci->di, st->effNBody);

    for (i = 0; i < 3; ++i)
    {
        nbb->pos[i] = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
        nbb->vel[i] = mwCreateZeroReadWriteBuffer(ci, st->effNBody * sizeof(real));
        nbb->acc[i] = mwCreateZeroReadWriteBuffer(ci, st->effNBody * sizeof(real));

        if (!nbb->pos[i] || !nbb->vel[i] || !nbb->acc[i])
        {
            return MW_CL_ERROR;
        }

        if (ctx->criterion != Exact)
        {
            nbb->min[i] = mwCreateZeroReadWriteBuffer(ci, ci->di.maxCompUnits * sizeof(real));
            nbb->max[i] = mwCreateZeroReadWriteBuffer(ci, ci->di.maxCompUnits * sizeof(real));
            if (!nbb->min[i] || !nbb->max[i])
            {
                return MW_CL_ERROR;
            }
        }

        if (ctx->useQuad)
        {
            nbb->quad.xx = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
            nbb->quad.xy = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
            nbb->quad.xz = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));

            nbb->quad.yy = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
            nbb->quad.yz = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));

            nbb->quad.zz = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
            if (!nbb->quad.xx || !nbb->quad.xy || !nbb->quad.xz || !nbb->quad.yy || !nbb->quad.yz || !nbb->quad.zz)
            {
                return MW_CL_ERROR;
            }

        }
    }

    massSize = st->usesExact ? st->effNBody * sizeof(real) : (nNode + 1) * sizeof(real);
    nbb->masses = mwCreateZeroReadWriteBuffer(ci, massSize);
    if (!nbb->masses)
    {
        return MW_CL_ERROR;
    }

    nbb->treeStatus = mwCreateZeroReadWriteBuffer(ci, sizeof(TreeStatus));
    if (!nbb->treeStatus)
    {
        return MW_CL_ERROR;
    }


    /* If we are doing an exact Nbody, we don't need the rest */
    if (ctx->criterion != Exact)
    {
        nbb->start = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(cl_int));
        nbb->count = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(cl_int));
        nbb->sort = mwCreateZeroReadWriteBuffer(ci, st->effNBody * sizeof(cl_int));
        nbb->child = mwCreateZeroReadWriteBuffer(ci, NSUB * (nNode + 1) * sizeof(cl_int));

        if (!nbb->start || !nbb->count || !nbb->sort || !nbb->child)
        {
            return MW_CL_ERROR;
        }

        if (ctx->criterion == SW93 || ctx->criterion == NewCriterion)
        {
            /* This only is for cells, so we could subtract nbody if we wanted */
            nbb->critRadii = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
            if (!nbb->critRadii)
            {
                return MW_CL_ERROR;
            }
        }
    }

    return CL_SUCCESS;
}

static cl_int nbMapBodies(real* pos[3], real* vel[3], real** mass, NBodyBuffers* nbb, CLInfo* ci, cl_map_flags flags, NBodyState* st)
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

static cl_int nbUnmapBodies(real* pos[3], real* vel[3], real* mass, NBodyBuffers* nbb, CLInfo* ci)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;

    for (i = 0; i < 3; ++i)
    {
        if (pos[i])
        {
            err |= clEnqueueUnmapMemObject(ci->queue, nbb->pos[i], pos[i], 0, NULL, NULL);
        }

        if (vel[i])
        {
            err |= clEnqueueUnmapMemObject(ci->queue, nbb->vel[i], vel[i], 0, NULL, NULL);
        }
    }

    if (mass)
    {
        err |= clEnqueueUnmapMemObject(ci->queue, nbb->masses, mass, 0, NULL, NULL);
    }

    return err;
}

/* If last parameter is true, copy to the buffers. if false, copy from the buffers */
cl_int nbMarshalBodies(NBodyState* st, cl_bool marshalIn)
{
    cl_int i;
    cl_int err = CL_SUCCESS;
    const Body* b;
    real* pos[3] = { NULL, NULL, NULL };
    real* vel[3] = { NULL, NULL, NULL };
    real* mass = NULL;
    CLInfo* ci = st->ci;
    NBodyBuffers* nbb = st->nbb;
    cl_map_flags flags = marshalIn ? CL_MAP_WRITE : CL_MAP_READ;

    err = nbMapBodies(pos, vel, &mass, nbb, ci, flags, st);
    if (err != CL_SUCCESS)
    {
        nbUnmapBodies(pos, vel, mass, nbb, ci);
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

        st->dirty = FALSE;
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

    return nbUnmapBodies(pos, vel, mass, nbb, ci);
}

void nbPrintKernelTimings(const NBodyState* st)
{
    double totalTime = 0.0;
    cl_uint i;
    double nStep = (double) (st->step + 1);
    const double* kernelTimings = st->workSizes->kernelTimings;

    for (i = 0; i < 7; ++i)
    {
        totalTime += kernelTimings[i];
    }

    mw_printf("\n--------------------------------------------------------------------------------\n"
              "Total timing over %d (+1) steps:\n"
              "                         Average             Total            Fraction\n"
              "                    ----------------   ----------------   ----------------\n"
              "  boundingBox:      %16f   %16f   %15.4f%%\n"
              "  buildTree:        %16f   %16f   %15.4f%%\n"
              "  summarization:    %16f   %16f   %15.4f%%\n"
              "  sort:             %16f   %16f   %15.4f%%\n"
              "  quad moments:     %16f   %16f   %15.4f%%\n"
              "  forceCalculation: %16f   %16f   %15.4f%%\n"
              "  integration:      %16f   %16f   %15.4f%%\n"
              "  ==============================================================================\n"
              "  total             %16f   %16f   %15.4f%%\n"
              "\n--------------------------------------------------------------------------------\n"
              "\n",
              st->step,
              kernelTimings[0] / nStep, kernelTimings[0], 100.0 * kernelTimings[0] / totalTime,
              kernelTimings[1] / nStep, kernelTimings[1], 100.0 * kernelTimings[1] / totalTime,
              kernelTimings[2] / nStep, kernelTimings[2], 100.0 * kernelTimings[2] / totalTime,
              kernelTimings[3] / nStep, kernelTimings[3], 100.0 * kernelTimings[3] / totalTime,
              kernelTimings[4] / nStep, kernelTimings[4], 100.0 * kernelTimings[4] / totalTime,
              kernelTimings[5] / nStep, kernelTimings[5], 100.0 * kernelTimings[5] / totalTime,
              kernelTimings[6] / nStep, kernelTimings[6], 100.0 * kernelTimings[6] / totalTime,
              totalTime / nStep,        totalTime,        100.0 * totalTime / totalTime
        );
}

void nbPrintKernelLimits(NBodyState* st)
{
    WGInfo wgi;
    CLInfo* ci = st->ci;
    NBodyKernels* kernels = st->kernels;

    mw_printf("Bounding box:\n");
    mwGetWorkGroupInfo(kernels->boundingBox, ci, &wgi);
    mwPrintWorkGroupInfo(&wgi);

    mw_printf("Tree Build:\n");
    mwGetWorkGroupInfo(kernels->buildTree, ci, &wgi);
    mwPrintWorkGroupInfo(&wgi);

    mw_printf("Summarization:\n");
    mwGetWorkGroupInfo(kernels->summarization, ci, &wgi);
    mwPrintWorkGroupInfo(&wgi);

    mw_printf("Sort:\n");
    mwGetWorkGroupInfo(kernels->sort, ci, &wgi);
    mwPrintWorkGroupInfo(&wgi);

    mw_printf("Quad moments:\n");
    mwGetWorkGroupInfo(kernels->quadMoments, ci, &wgi);
    mwPrintWorkGroupInfo(&wgi);

    mw_printf("Force calculation:\n");
    mwGetWorkGroupInfo(kernels->forceCalculation, ci, &wgi);
    mwPrintWorkGroupInfo(&wgi);

    mw_printf("Integration:\n");
    mwGetWorkGroupInfo(kernels->integration, ci, &wgi);
    mwPrintWorkGroupInfo(&wgi);

    mw_printf("Force calculation (Exact):\n");
    mwGetWorkGroupInfo(kernels->forceCalculation_Exact, ci, &wgi);
    mwPrintWorkGroupInfo(&wgi);
}


NBodyStatus nbRunSystemCL(const NBodyCtx* ctx, NBodyState* st)
{
    NBodyStatus rc;
    cl_int err;

    rc = nbMainLoopCL(ctx, st);
    if (nbStatusIsFatal(rc))
    {
        return rc;
    }

    fflush(stdout); /* Try to prevent some of the GPU printfs from getting lost */
    fflush(stderr);

    err = nbMarshalBodies(st, CL_FALSE);
    if (err != CL_SUCCESS)
    {
        return NBODY_CL_ERROR;
    }

    nbPrintKernelTimings(st);

    return rc;
}

