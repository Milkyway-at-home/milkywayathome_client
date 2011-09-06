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

#include "milkyway_util.h"
#include "milkyway_extra.h"
#include "milkyway_cpp_util.h"
#include "milkyway_cl.h"
#include "mw_cl.h"
#include "setup_cl.h"
#include "separation_cl_buffers.h"
#include "separation_binaries.h"
#include "cl_compile_flags.h"


#include <assert.h>

#ifdef _WIN32
  #include <direct.h>
#endif /* _WIN32 */


#if SEPARATION_INLINE_KERNEL
extern char* inlinedIntegralKernelSrc;
#else
static char* inlinedIntegralKernelSrc = NULL;
#endif /* SEPARATION_INLINE_KERNEL */

cl_bool haveInlinedKernel()
{
    return (inlinedIntegralKernelSrc != NULL);
}

static void printRunSizes(const RunSizes* sizes, const IntegralArea* ia)
{
    warn("Range:          { nu_steps = %u, mu_steps = %u, r_steps = %u }\n"
         "Iteration area: "LLU"\n"
         "Chunk estimate: "ZU"\n"
         "Num chunks:     "ZU"\n"
         "Chunk size:     "ZU"\n"
         "Added area:     %u\n"
         "Effective area: "LLU"\n"
         ,
         ia->nu_steps, ia->mu_steps, ia->r_steps,
         sizes->area,
         sizes->nChunkEstimate,
         sizes->nChunk,
         sizes->chunkSize,
         sizes->extra,
         sizes->effectiveArea);
}

static cl_double estimateWUGFLOPsPerIter(const AstronomyParameters* ap, const IntegralArea* ia)
{
    cl_ulong perItem, perIter;
    cl_ulong tmp = 32 + ap->number_streams * 68;
    if (ap->aux_bg_profile)
        tmp += 8;

    perItem = tmp * ap->convolve + 1 + (ap->number_streams * 2);
    perIter = perItem * ia->mu_steps * ia->r_steps;

    return 1.0e-9 * (cl_double) perIter;
}

#define GPU_EFFICIENCY_ESTIMATE (0.95)

/* Based on the flops of the device and workunit, pick a target number of chunks */
static cl_uint findNChunk(const AstronomyParameters* ap,
                          const IntegralArea* ia,
                          const DevInfo* di,
                          const CLRequest* clr)
{
    cl_double gflops = deviceEstimateGFLOPs(di, DOUBLEPREC);
    cl_double effFlops = GPU_EFFICIENCY_ESTIMATE * (cl_double) gflops;
    cl_double iterFlops = estimateWUGFLOPsPerIter(ap, ia);

    cl_double estIterTime = 1000.0 * (cl_double) iterFlops / effFlops; /* milliseconds */

    cl_double timePerIter = 1000.0 / clr->targetFrequency;

    cl_uint nChunk = (cl_uint) (estIterTime / timePerIter);

    return nChunk == 0 ? 1 : nChunk;
}

/* Returns CL_TRUE on error */
cl_bool findRunSizes(RunSizes* sizes,
                     const CLInfo* ci,
                     const DevInfo* di,
                     const AstronomyParameters* ap,
                     const IntegralArea* ia,
                     const CLRequest* clr)
{
    WGInfo wgi;
    cl_int err;
    size_t nWavefrontPerCU;
    size_t blockSize; /* Size chunks should be multiples of */
    cl_bool forceOneChunk = clr->nonResponsive || di->nonOutput;

    /* I assume this is how this works for 1D limit */
    const cl_ulong maxWorkDim = (cl_ulong) di->maxWorkItemSizes[0] * di->maxWorkItemSizes[1] * di->maxWorkItemSizes[2];
    const cl_ulong r = (cl_ulong) ia->r_steps;
    const cl_ulong mu = (cl_ulong) ia->mu_steps;

    sizes->r = ia->r_steps;
    sizes->mu = ia->mu_steps;
    sizes->nu = ia->nu_steps;
    sizes->area = r * mu;

    if (r > CL_ULONG_MAX / mu)
    {
        warn("Integral area overflows cl_ulong\n");
        return CL_TRUE;
    }

    err = mwGetWorkGroupInfo(ci, &wgi);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to get work group info", err);
        return CL_TRUE;
    }

    if (clr->verbose)
    {
        mwPrintWorkGroupInfo(&wgi);
    }

    if (!mwDivisible(wgi.wgs, (size_t) di->warpSize))
    {
        warn("Kernel reported work group size ("ZU") not a multiple of warp size (%u)\n",
             wgi.wgs,
             di->warpSize);
        return CL_TRUE;
    }

    /* This should give a good occupancy. If the global size isn't a
     * multiple of this bad performance things happen. */
    nWavefrontPerCU = wgi.wgs / di->warpSize;

    /* Since we don't use any workgroup features, it makes sense to
     * use the wavefront size as the workgroup size */
    sizes->local[0] = di->warpSize;


    /* For maximum efficiency, we want global work sizes to be multiples of
     * (warp size) * (number compute units) * (number of warps for good occupancy)
     * Then we throw in another factor since we can realistically do more work at once
     */

    blockSize = nWavefrontPerCU * di->warpSize * di->maxCompUnits;
    {
        cl_uint magic = 1;
        sizes->nChunkEstimate = findNChunk(ap, ia, di, clr);

        /* If specified and acceptable, use a user specified factor for the
         * number of blocks to use. Otherwise, make a guess appropriate for the hardware. */

        if (clr->magicFactor < 0)
        {
            warn("Invalid magic factor %d. Magic factor must be >= 0\n", clr->magicFactor);
        }

        if (clr->magicFactor <= 0) /* Use default calculation */
        {
            /*   m * b ~= area / n   */
            magic = sizes->area / (sizes->nChunkEstimate * blockSize);
            if (magic == 0)
                magic = 1;
        }
        else   /* Use user setting */
        {
            magic = (cl_uint) clr->magicFactor;
        }

        sizes->chunkSize = magic * blockSize;
    }

    sizes->effectiveArea = sizes->chunkSize * mwDivRoundup(sizes->area, sizes->chunkSize);
    sizes->nChunk = forceOneChunk ? 1 : mwDivRoundup(sizes->effectiveArea, sizes->chunkSize);
    sizes->extra = (cl_uint) (sizes->effectiveArea - sizes->area);

    if (sizes->nChunk == 1) /* Magic factor probably too high or very small workunit, or nonresponsive */
    {
        /* Like using magic == 1 */
        sizes->effectiveArea = blockSize * mwDivRoundup(sizes->area, blockSize);
        sizes->chunkSize = sizes->effectiveArea;
        sizes->extra = sizes->effectiveArea - sizes->area;
    }

    warn("Using a block size of "ZU" with a magic factor of "ZU"\n",
         blockSize,
         sizes->chunkSize / blockSize);

    sizes->chunkSize = sizes->effectiveArea / sizes->nChunk;

    /* We should be hitting memory size limits before we ever get to this */
    if (sizes->chunkSize > maxWorkDim)
    {
        warn("Warning: Area too large for one chunk (max size = "LLU")\n", maxWorkDim);
        while (sizes->chunkSize > maxWorkDim)
        {
            sizes->nChunk *= 2;
            sizes->chunkSize = sizes->effectiveArea / sizes->nChunk;
        }

        if (!mwDivisible(sizes->chunkSize, sizes->local[0]))
        {
            warn("FIXME: I'm too lazy to handle very large workunits properly\n");
            return CL_TRUE;
        }
        else if (!mwDivisible(sizes->chunkSize, blockSize))
        {
            warn("FIXME: Very large workunit potentially slower than it should be\n");
        }
    }

    sizes->global[0] = sizes->chunkSize;

    printRunSizes(sizes, ia);

    if (sizes->effectiveArea < sizes->area)
    {
        warn("Effective area less than actual area!\n");
        return CL_TRUE;
    }

    return CL_FALSE;
}


/* Only sets the constant arguments, not the outputs which we double buffer */
cl_int separationSetKernelArgs(CLInfo* ci, SeparationCLMem* cm, const RunSizes* runSizes)
{
    cl_int err = CL_SUCCESS;

    /* Set output buffer arguments */
    err |= clSetKernelArg(ci->kern, 0, sizeof(cl_mem), &cm->outBg);
    err |= clSetKernelArg(ci->kern, 1, sizeof(cl_mem), &cm->outStreams);

    /* The constant, global arguments */
    err |= clSetKernelArg(ci->kern, 2, sizeof(cl_mem), &cm->rc);
    err |= clSetKernelArg(ci->kern, 3, sizeof(cl_mem), &cm->rPts);
    err |= clSetKernelArg(ci->kern, 4, sizeof(cl_mem), &cm->lTrig);
    err |= clSetKernelArg(ci->kern, 5, sizeof(cl_mem), &cm->bSin);

    /* The __constant arguments */
    err |= clSetKernelArg(ci->kern, 6, sizeof(cl_mem), &cm->ap);
    err |= clSetKernelArg(ci->kern, 7, sizeof(cl_mem), &cm->sc);
    err |= clSetKernelArg(ci->kern, 8, sizeof(cl_mem), &cm->sg_dx);

    err |= clSetKernelArg(ci->kern, 9, sizeof(cl_uint), &runSizes->extra);
    err |= clSetKernelArg(ci->kern, 10, sizeof(cl_uint), &runSizes->r);
    err |= clSetKernelArg(ci->kern, 11, sizeof(cl_uint), &runSizes->mu);
    err |= clSetKernelArg(ci->kern, 12, sizeof(cl_uint), &runSizes->nu);

    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error setting kernel arguments", err);
        return err;
    }

    return CL_SUCCESS;
}

/* Reading from the file is more convenient for actually working on
 * it. Inlining is more useful for releasing when we don't want to
 * deal with the hassle of distributing more files. */
char* findKernelSrc()
{
    char* kernelSrc = NULL;

    /* Try different alternatives */

    if (inlinedIntegralKernelSrc)
        return inlinedIntegralKernelSrc;

    kernelSrc = mwReadFileResolved("AllInOneFile.cl");
    if (kernelSrc)
        return kernelSrc;

    kernelSrc = mwReadFileResolved("integrals.cl");
    if (kernelSrc)
        return kernelSrc;

    kernelSrc = mwReadFile("../kernels/integrals.cl");
    if (kernelSrc)
        return kernelSrc;

    warn("Failed to read kernel file\n");

    return kernelSrc;
}

void freeKernelSrc(char* src)
{
    if (src != inlinedIntegralKernelSrc)
        free(src);
}

#define NUM_CONST_BUF_ARGS 5

/* Check that the device has the necessary resources */
static cl_bool separationCheckDevMemory(const DevInfo* di, const SeparationSizes* sizes)
{
    size_t totalConstBuf;
    size_t totalGlobalConst;
    size_t totalOut;
    size_t totalMem;

    totalOut = sizes->outBg + sizes->outStreams;
    totalConstBuf = sizes->ap + sizes->ia + sizes->sc + sizes->sg_dx;
    totalGlobalConst = sizes->lTrig + sizes->bSin + sizes->rPts + sizes->rc;

    totalMem = totalOut + totalConstBuf + totalGlobalConst;
    if (totalMem > di->memSize)
    {
        warn("Total required device memory ("ZU") > available ("LLU")\n", totalMem, di->memSize);
        return CL_FALSE;
    }

    /* Check individual allocations. Right now ATI has a fairly small
     * maximum allowed allocation compared to the actual memory
     * available. */
    if (totalOut > di->memSize)
    {
        warn("Device has insufficient global memory for output buffers\n");
        return CL_FALSE;
    }

    if (sizes->outBg > di->maxMemAlloc || sizes->outStreams > di->maxMemAlloc)
    {
        warn("An output buffer would exceed CL_DEVICE_MAX_MEM_ALLOC_SIZE\n");
        return CL_FALSE;
    }

    if (   sizes->lTrig > di->maxMemAlloc
        || sizes->bSin > di->maxMemAlloc
        || sizes->rPts > di->maxMemAlloc
        || sizes->rc > di->maxMemAlloc)
    {
        warn("A global constant buffer would exceed CL_DEVICE_MAX_MEM_ALLOC_SIZE\n");
        return CL_FALSE;
    }

    if (NUM_CONST_BUF_ARGS > di->maxConstArgs)
    {
        warn("Need more constant arguments than available\n");
        return CL_FALSE;
    }

    if (totalConstBuf > di-> maxConstBufSize)
    {
        warn("Device doesn't have enough constant buffer space\n");
        return CL_FALSE;
    }

    return CL_TRUE;
}

/* TODO: Should probably check for likelihood also */
cl_bool separationCheckDevCapabilities(const DevInfo* di, const AstronomyParameters* ap, const IntegralArea* ias)
{
    cl_int i;
    SeparationSizes sizes;

  #if DOUBLEPREC
    if (!mwSupportsDoubles(di))
    {
        warn("Device doesn't support double precision\n");
        return MW_CL_ERROR;
    }
  #endif /* DOUBLEPREC */

    for (i = 0; i < ap->number_integrals; ++i)
    {
        calculateSizes(&sizes, ap, &ias[i]);
        if (!separationCheckDevMemory(di, &sizes))
        {
            warn("Capability check failed for cut %u\n", i);
            return CL_FALSE;
        }
    }

    return CL_TRUE;
}

/* Estimate time for a nu step in milliseconds */
cl_double cudaEstimateIterTime(const DevInfo* di, cl_double flopsPerIter, cl_double flops)
{
    cl_double devFactor;

    /* Experimentally determined constants */
    devFactor = computeCapabilityIs(di, 1, 3) ? 1.87 : 1.53;

    /* Idea is this is a sort of efficiency factor for the
     * architecture vs. the theoretical FLOPs. We can then scale by
     * the theoretical flops compared to the reference devices. */

    return 1000.0 * devFactor * flopsPerIter / flops;
}

cl_int setupSeparationCL(CLInfo* ci,
                         const AstronomyParameters* ap,
                         const IntegralArea* ias,
                         const CLRequest* clr)
{
    cl_int err;
    char* compileFlags;
    char* kernelSrc;

    err = mwSetupCL(ci, clr);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error getting device and context", err);
        return err;
    }

    err = mwGetDevInfo(&ci->di, ci->dev);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get device info\n");
        return err;
    }

    if (clr->verbose)
    {
        mwPrintDevInfo(&ci->di);
    }
    else
    {
        mwPrintDevInfoShort(&ci->di);
    }

    if (!separationCheckDevCapabilities(&ci->di, ap, ias))
    {
        warn("Device failed capability check\n");
        return MW_CL_ERROR;
    }

    compileFlags = getCompilerFlags(ci, ap, clr);
    if (!compileFlags)
    {
        warn("Failed to get compiler flags\n");
        return MW_CL_ERROR;
    }

    kernelSrc = findKernelSrc();
    if (!kernelSrc)
    {
        warn("Failed to read CL kernel source\n");
        return MW_CL_ERROR;
    }

    warn("\nCompiler flags:\n%s\n\n", compileFlags);
    err = mwSetProgramFromSrc(ci, "probabilities", (const char**) &kernelSrc, 1, compileFlags);

    freeKernelSrc(kernelSrc);
    free(compileFlags);

    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating program from source", err);
        return err;
    }

    return CL_SUCCESS;
}

