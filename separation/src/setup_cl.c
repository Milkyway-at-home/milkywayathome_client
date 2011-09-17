/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
#include "replace_amd_il.h"
#include "il_kernels.h"

#include <assert.h>

#ifdef _WIN32
  #include <direct.h>
#endif /* _WIN32 */


#if SEPARATION_INLINE_KERNEL
extern char* inlinedIntegralKernelSrc;
#else
static char* inlinedIntegralKernelSrc = NULL;
#endif /* SEPARATION_INLINE_KERNEL */

cl_kernel _separationKernel = NULL;

static cl_int createSeparationKernel(const CLInfo* ci)
{
    cl_int err;

    _separationKernel = clCreateKernel(ci->prog, "mu_sum_kernel", &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating kernel '%s'", err, "mu_sum_kernel");
    }

    return err;
}

cl_int releaseSeparationKernel()
{
    cl_int err = CL_SUCCESS;

    if (_separationKernel)
        err = clReleaseKernel(_separationKernel);

    return err;
}

cl_bool haveInlinedKernel()
{
    return (inlinedIntegralKernelSrc != NULL);
}

static void printRunSizes(const RunSizes* sizes, const IntegralArea* ia, cl_bool verbose)
{
    mw_printf("Range:          { nu_steps = %u, mu_steps = %u, r_steps = %u }\n"
              "Iteration area: "LLU"\n"
              "Chunk estimate: "ZU"\n"
              "Num chunks:     "ZU"\n"
              "Chunk size:     "ZU"\n"
              "Added area:     %u\n"
              "Effective area: "LLU"\n",
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
        mw_printf("Integral area overflows cl_ulong\n");
        return CL_TRUE;
    }

    if (di->devType == CL_DEVICE_TYPE_CPU)
    {
        sizes->nChunk = sizes->nChunkEstimate = 1;
        sizes->chunkSize = sizes->effectiveArea = sizes->area;
        sizes->extra = 0;

        sizes->local[0] = 1;
        sizes->global[0] = sizes->area;

        return CL_FALSE;
    }

    err = mwGetWorkGroupInfo(_separationKernel, ci, &wgi);
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
        mw_printf("Kernel reported work group size ("ZU") not a multiple of warp size (%u)\n",
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
            mw_printf("Invalid magic factor %d. Magic factor must be >= 0\n", clr->magicFactor);
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

    mw_printf("Using a block size of "ZU" with a magic factor of "ZU"\n",
              blockSize,
              sizes->chunkSize / blockSize);

    sizes->chunkSize = sizes->effectiveArea / sizes->nChunk;

    /* We should be hitting memory size limits before we ever get to this */
    if (sizes->chunkSize > maxWorkDim)
    {
        mw_printf("Warning: Area too large for one chunk (max size = "LLU")\n", maxWorkDim);
        while (sizes->chunkSize > maxWorkDim)
        {
            sizes->nChunk *= 2;
            sizes->chunkSize = sizes->effectiveArea / sizes->nChunk;
        }

        if (!mwDivisible(sizes->chunkSize, sizes->local[0]))
        {
            mw_printf("FIXME: I'm too lazy to handle very large workunits properly\n");
            return CL_TRUE;
        }
        else if (!mwDivisible(sizes->chunkSize, blockSize))
        {
            mw_printf("FIXME: Very large workunit potentially slower than it should be\n");
        }
    }

    sizes->global[0] = sizes->chunkSize;

    printRunSizes(sizes, ia, clr->verbose);

    if (sizes->effectiveArea < sizes->area)
    {
        mw_printf("Effective area less than actual area!\n");
        return CL_TRUE;
    }

    return CL_FALSE;
}


/* Only sets the constant arguments, not the outputs which we double buffer */
cl_int separationSetKernelArgs(CLInfo* ci, SeparationCLMem* cm, const RunSizes* runSizes)
{
    cl_int err = CL_SUCCESS;

    /* Set output buffer arguments */
    err |= clSetKernelArg(_separationKernel, 0, sizeof(cl_mem), &cm->outBg);
    err |= clSetKernelArg(_separationKernel, 1, sizeof(cl_mem), &cm->outStreams);

    /* The constant, global arguments */
    err |= clSetKernelArg(_separationKernel, 2, sizeof(cl_mem), &cm->rc);
    err |= clSetKernelArg(_separationKernel, 3, sizeof(cl_mem), &cm->rPts);
    err |= clSetKernelArg(_separationKernel, 4, sizeof(cl_mem), &cm->lTrig);
    err |= clSetKernelArg(_separationKernel, 5, sizeof(cl_mem), &cm->bSin);

    /* The __constant arguments */
    err |= clSetKernelArg(_separationKernel, 6, sizeof(cl_mem), &cm->ap);
    err |= clSetKernelArg(_separationKernel, 7, sizeof(cl_mem), &cm->sc);
    err |= clSetKernelArg(_separationKernel, 8, sizeof(cl_mem), &cm->sg_dx);

    err |= clSetKernelArg(_separationKernel, 9, sizeof(cl_uint), &runSizes->extra);
    err |= clSetKernelArg(_separationKernel, 10, sizeof(cl_uint), &runSizes->r);
    err |= clSetKernelArg(_separationKernel, 11, sizeof(cl_uint), &runSizes->mu);
    err |= clSetKernelArg(_separationKernel, 12, sizeof(cl_uint), &runSizes->nu);

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

    mw_printf("Failed to read kernel file\n");

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
        mw_printf("Total required device memory ("ZU") > available ("LLU")\n", totalMem, di->memSize);
        return CL_FALSE;
    }

    /* Check individual allocations. Right now ATI has a fairly small
     * maximum allowed allocation compared to the actual memory
     * available. */
    if (totalOut > di->memSize)
    {
        mw_printf("Device has insufficient global memory for output buffers\n");
        return CL_FALSE;
    }

    if (sizes->outBg > di->maxMemAlloc || sizes->outStreams > di->maxMemAlloc)
    {
        mw_printf("An output buffer would exceed CL_DEVICE_MAX_MEM_ALLOC_SIZE\n");
        return CL_FALSE;
    }

    if (   sizes->lTrig > di->maxMemAlloc
        || sizes->bSin > di->maxMemAlloc
        || sizes->rPts > di->maxMemAlloc
        || sizes->rc > di->maxMemAlloc)
    {
        mw_printf("A global constant buffer would exceed CL_DEVICE_MAX_MEM_ALLOC_SIZE\n");
        return CL_FALSE;
    }

    if (NUM_CONST_BUF_ARGS > di->maxConstArgs)
    {
        mw_printf("Need more constant arguments than available\n");
        return CL_FALSE;
    }

    if (totalConstBuf > di-> maxConstBufSize)
    {
        mw_printf("Device doesn't have enough constant buffer space\n");
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
        mw_printf("Device doesn't support double precision\n");
        return MW_CL_ERROR;
    }
  #endif /* DOUBLEPREC */

    for (i = 0; i < ap->number_integrals; ++i)
    {
        calculateSizes(&sizes, ap, &ias[i]);
        if (!separationCheckDevMemory(di, &sizes))
        {
            mw_printf("Capability check failed for cut %u\n", i);
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

static char* replaceUAVIds(const char* ilSrc, size_t* lenOut, ...)
{
    char* buf;
    va_list argPtr;
    size_t len;
    int rc;

    len = strlen(ilSrc);
    buf = mwMalloc(len + 1);
    buf[len] = '\0';

    va_start(argPtr, lenOut);
    rc = vsprintf(buf, ilSrc, argPtr);
    va_end(argPtr);

    if ((size_t) rc != len)
    {
        free(buf);
        return NULL;
    }

    if (lenOut)
    {
        *lenOut = len;
    }

    return buf;
}

/* Different UAV IDs are used for different GPUs, and it must match */
static char* getILSrc(int nStream, MWCALtargetEnum target, size_t* len)
{
    char* ilSrc = NULL;
    cl_uint u = uavIdFromMWCALtargetEnum(target);

    /* Should we be checking which UAV is used from the binary? */
    if (u > 99)
    {
        /* We rely on the size of a 2 digit UAV id being the same as the '%d' format,
           but there are only a few UAVs anyway
         */
        mw_printf("UAV id %u is absurd\n", u);
        return NULL;
    }

    /* This is pretty special */
    switch (nStream)
    {
        case 1:   /* 9 */
            ilSrc = replaceUAVIds(ilKernelSrc1, len, u, u, u, u, u, u, u, u, u);
            break;

        case 2:  /* 11 */
            ilSrc = replaceUAVIds(ilKernelSrc2, len, u, u, u, u, u, u, u, u, u, u, u);
            break;

        case 3:  /* 13 */
            ilSrc = replaceUAVIds(ilKernelSrc3, len, u, u, u, u, u, u, u, u, u, u, u, u, u);
            break;

        case 4:  /* 15 */
            ilSrc = replaceUAVIds(ilKernelSrc4, len, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u);
            break;

        default:
            mw_unreachable();
    }

    if (!ilSrc)
    {
        mw_printf("Error getting processed IL kernel source\n");
    }

    return ilSrc;
}

static cl_int setProgramFromILKernel(CLInfo* ci, const AstronomyParameters* ap, const CLRequest* clr)
{
    unsigned char* bin;
    unsigned char* modBin;
    char* ilSrc;
    size_t ilLen = 0;
    size_t binSize = 0;
    size_t modBinSize = 0;
    cl_int err;

    bin = mwGetProgramBinary(ci, &binSize);
    if (!bin)
    {
        return MW_CL_ERROR;
    }

    err = clReleaseProgram(ci->prog);
    if (err != CL_SUCCESS)
    {
        free(bin);
        return err;
    }

    ilSrc = getILSrc(ap->number_streams, ci->di.calTarget, &ilLen);
    modBin = getModifiedAMDBinary(bin, binSize, ilSrc, ilLen, &modBinSize);
    free(bin);
    free(ilSrc);

    if (!modBin)
    {
        mw_printf("Error getting modified binary or IL source\n");
        return MW_CL_ERROR;
    }

    err = mwSetProgramFromBin(ci, modBin, modBinSize);
    free(modBin);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating program from binary", err);
        return err;
    }

    return CL_SUCCESS;
}

static cl_bool isILKernelTarget(const DevInfo* di)
{
    MWCALtargetEnum t = di->calTarget;

    return (t == MW_CAL_TARGET_770) || (t == MW_CAL_TARGET_CYPRESS) || (t == MW_CAL_TARGET_CAYMAN);
}

static cl_bool usingILKernelIsAcceptable(const CLInfo* ci, const AstronomyParameters* ap, const CLRequest* clr)
{
    const DevInfo* di = &ci->di;
    static const cl_int maxILKernelStreams = 4;

    if (!DOUBLEPREC)
        return CL_FALSE;

      if (clr->forceNoILKernel)
          return CL_FALSE;

    /* Supporting these unused options with the IL kernel is too much work */
    if (ap->number_streams > maxILKernelStreams || ap->aux_bg_profile)
        return CL_FALSE;

    /* Make sure an acceptable device */
    return (isAMDGPUDevice(di) && isILKernelTarget(di) && mwPlatformSupportsAMDOfflineDevices(ci));
}

cl_int setupSeparationCL(CLInfo* ci,
                         const AstronomyParameters* ap,
                         const IntegralArea* ias,
                         const CLRequest* clr)
{
    char* compileFlags;
    char* kernelSrc;
    cl_bool useILKernel;
    cl_int err;

    err = mwSetupCL(ci, clr);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error getting device and context", err);
        return err;
    }

    if (!separationCheckDevCapabilities(&ci->di, ap, ias))
    {
        mw_printf("Device failed capability check\n");
        return MW_CL_ERROR;
    }

    useILKernel = usingILKernelIsAcceptable(ci, ap, clr);
    kernelSrc = findKernelSrc();
    compileFlags = getCompilerFlags(ci, ap, useILKernel);

    if (!compileFlags || !kernelSrc)
    {
        mw_printf("Failed to get kernel source or compiler flags\n");
        err = MW_CL_ERROR;
        goto setup_exit;
    }

    mw_printf("\nCompiler flags:\n%s\n\n", compileFlags);
    err = mwSetProgramFromSrc(ci, (const char**) &kernelSrc, 1, compileFlags);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating program from source", err);
        goto setup_exit;
    }

    if (useILKernel)
    {
        mw_printf("Using AMD IL kernel\n");
        err = setProgramFromILKernel(ci, ap, clr);
        if (err != CL_SUCCESS)
        {
            /* Recompiles again but I don't really care. */
            mw_printf("Failed to create IL kernel. Falling back to source kernel\n");
            err = mwSetProgramFromSrc(ci, (const char**) &kernelSrc, 1, compileFlags);
            if (err != CL_SUCCESS)
            {
                mwCLWarn("Error creating program from source", err);
            }
        }
    }

    if (err == CL_SUCCESS)
    {
        err = mwCreateKernel(&_separationKernel, ci, "probabilities");
    }

setup_exit:
    freeKernelSrc(kernelSrc);
    free(compileFlags);

    return err;
}

