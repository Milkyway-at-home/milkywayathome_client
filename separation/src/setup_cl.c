/*
 *  Copyright (c) 2010-2012 Matthew Arsenault
 *  Copyright (c) 2010-2012 Rensselaer Polytechnic Institute
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
#include "milkyway_cl.h"
#include "setup_cl.h"
#include "separation_cl_buffers.h"
#include "separation_binaries.h"
#include "cl_compile_flags.h"
#include "replace_amd_il.h"

#include <assert.h>

#ifdef _WIN32
  #include <direct.h>
#endif /* _WIN32 */


static cl_program integrationProgram = NULL;
static cl_program summarizationProgram = NULL;


extern const unsigned char probabilities_kernel_cl[];
extern const size_t probabilities_kernel_cl_len;

extern const unsigned char summarization_kernel_cl[];
extern const size_t summarization_kernel_cl_len;


cl_kernel _separationKernel = NULL;
cl_kernel _summarizationKernel = NULL;

size_t _summarizationWorkgroupSize = 0;

cl_int releaseSeparationKernel(void)
{
    cl_int err = CL_SUCCESS;

    if (_separationKernel)
        err = clReleaseKernel(_separationKernel);

    if (_summarizationKernel)
        err |= clReleaseKernel(_summarizationKernel);

    if (integrationProgram)
        err |= clReleaseProgram(integrationProgram);

    if (integrationProgram)
        err |= clReleaseProgram(summarizationProgram);

    return err;
}

static void printRunSizes(const RunSizes* sizes, const IntegralArea* ia)
{
    mw_printf("Range:          { nu_steps = %u, mu_steps = %u, r_steps = %u }\n"
              "Iteration area: "LLU"\n"
              "Chunk estimate: "ZU"\n"
              "Num chunks:     "ZU"\n"
              "Chunk size:     "ZU"\n"
              "Added area:     %u\n"
              "Effective area: "LLU"\n"
              "Initial wait:   %u ms\n",
              ia->nu_steps, ia->mu_steps, ia->r_steps,
              sizes->area,
              sizes->nChunkEstimate,
              sizes->nChunk,
              sizes->chunkSize,
              sizes->extra,
              sizes->effectiveArea,
              sizes->initialWait
        );
}

static cl_double estimateWUGFLOPsPerIter(const AstronomyParameters* ap, const IntegralArea* ia)
{
    cl_ulong perItem, perIter;
    cl_ulong loop;

    loop = 58 * ap->number_streams + 56;
    if (ap->aux_bg_profile)
        loop += 8;

    perItem = ap->convolve * loop + (2 * ap->number_streams) + 4;
    perIter = perItem * ia->mu_steps * ia->r_steps;

    return 1.0e-9 * (cl_double) perIter;
}

/* Somewhat bullshit factor because of instructions not related to
 * actual work + other possible inefficiencies */
#define GPU_EFFICIENCY_ESTIMATE (0.80)

/* Based on the flops of the device and workunit, pick a target number of chunks */
static cl_uint findNChunk(const AstronomyParameters* ap,
                          const IntegralArea* ia,
                          const DevInfo* di,
                          const CLRequest* clr,
                          cl_uint* initialWaitTime)
{
    cl_double gflops = mwDeviceEstimateGFLOPs(di, DOUBLEPREC);
    cl_double effFlops = GPU_EFFICIENCY_ESTIMATE * (cl_double) gflops;
    cl_double iterFlops = estimateWUGFLOPsPerIter(ap, ia);

    cl_double estIterTime = 1000.0 * (cl_double) iterFlops / effFlops; /* milliseconds */

    cl_double timePerIter = 1000.0 / clr->targetFrequency;

    cl_uint nChunk = (cl_uint) (estIterTime / timePerIter);

    if (initialWaitTime)
    {
        /* Sleep for (wait factor) * estimated time before polling */
        *initialWaitTime = (cl_uint) (clr->gpuWaitFactor * estIterTime / nChunk);
    }

    return nChunk == 0 ? 1 : nChunk;
}

static void printPollMode(const CLInfo* ci, const RunSizes* sizes)
{
    cl_int mode = ci->pollingMode;

    if (mode == MW_POLL_CL_WAIT_FOR_EVENTS)
    {
        mw_printf("Using clWaitForEvents() for polling (mode %d)\n", mode);
    }
    else if (mode == MW_POLL_SLEEP_CL_WAIT_FOR_EVENTS)
    {
        mw_printf("Using clWaitForEvents() for polling with initial wait of %u ms (mode %d)\n",
                  sizes->initialWait,
                  mode
            );
    }
    else
    {
        mw_printf("Using manual polling with initial wait of %u ms (mode %d)\n",
                  sizes->initialWait,
                  mode);
    }
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
    cl_bool forceOneChunk = clr->nonResponsive || di->nonOutput || di->hasGraphicsQOS;

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
        mwPerrorCL(err, "Failed to get work group info");
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
        cl_uint nBlockPerChunk = 1;
        sizes->nChunkEstimate = findNChunk(ap, ia, di, clr, &sizes->initialWait);

        /* Make a guess appropriate for the hardware. */

        /* m * b ~= area / n   */
        nBlockPerChunk = sizes->area / (sizes->nChunkEstimate * blockSize);
        if (nBlockPerChunk == 0)
            nBlockPerChunk = 1;

        sizes->chunkSize = nBlockPerChunk * blockSize;
    }

    sizes->effectiveArea = sizes->chunkSize * mwDivRoundup(sizes->area, sizes->chunkSize);
    sizes->nChunk = forceOneChunk ? 1 : mwDivRoundup(sizes->effectiveArea, sizes->chunkSize);
    sizes->extra = (cl_uint) (sizes->effectiveArea - sizes->area);

    if (sizes->nChunk == 1) /* BlockPerChunk factor probably too high or very small workunit, or nonresponsive */
    {
        /* Like using nBlockPerChunk == 1 */
        sizes->effectiveArea = blockSize * mwDivRoundup(sizes->area, blockSize);
        sizes->chunkSize = sizes->effectiveArea;
        sizes->extra = sizes->effectiveArea - sizes->area;
    }

    mw_printf("Using a target frequency of %.1f\n"
              "Using a block size of "ZU" with "ZU" blocks/chunk\n",
              clr->targetFrequency,
              blockSize,
              sizes->chunkSize / blockSize
        );
    printPollMode(ci, sizes);

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

    printRunSizes(sizes, ia);

    if (sizes->effectiveArea < sizes->area)
    {
        mw_printf("Effective area less than actual area!\n");
        return CL_TRUE;
    }

    return CL_FALSE;
}


/* Only sets the constant arguments, not the outputs which we double buffer */
cl_int separationSetKernelArgs(SeparationCLMem* cm, const RunSizes* runSizes)
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

    err |= clSetKernelArg(_separationKernel, 9,  sizeof(cl_uint), &runSizes->extra);
    err |= clSetKernelArg(_separationKernel, 10, sizeof(cl_uint), &runSizes->r);
    err |= clSetKernelArg(_separationKernel, 11, sizeof(cl_uint), &runSizes->mu);
    err |= clSetKernelArg(_separationKernel, 12, sizeof(cl_uint), &runSizes->nu);

    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error setting kernel arguments");
        return err;
    }

    return CL_SUCCESS;
}

#define NUM_CONST_BUF_ARGS 5

/* Check that the device has the necessary resources */
static cl_bool separationCheckCutMemory(const DevInfo* di, const SeparationSizes* sizes)
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

static cl_bool separationCheckDevMemory(const DevInfo* di, const AstronomyParameters* ap, const IntegralArea* ias)
{
    cl_int i;
    SeparationSizes sizes;

    for (i = 0; i < ap->number_integrals; ++i)
    {
        calculateSizes(&sizes, ap, &ias[i]);
        if (!separationCheckCutMemory(di, &sizes))
        {
            mw_printf("Capability check failed for cut %u\n", i);
            return CL_FALSE;
        }
    }

    return CL_TRUE;
}

static cl_bool separationCheckDevCapabilities(const DevInfo* di)
{
    if (DOUBLEPREC && !mwSupportsDoubles(di))
    {
        mw_printf("Device doesn't support double precision\n");
        return CL_FALSE;
    }

    return CL_TRUE;
}

/* Estimate time for a nu step in milliseconds */
cl_double cudaEstimateIterTime(const DevInfo* di, cl_double flopsPerIter, cl_double flops)
{
    cl_double devFactor;

    /* Experimentally determined constants */
    devFactor = mwComputeCapabilityIs(di, 1, 3) ? 1.87 : 1.53;

    /* Idea is this is a sort of efficiency factor for the
     * architecture vs. the theoretical FLOPs. We can then scale by
     * the theoretical flops compared to the reference devices. */

    return 1000.0 * devFactor * flopsPerIter / flops;
}

static cl_int setProgramFromILKernel(CLInfo* ci, const AstronomyParameters* ap)
{
    unsigned char* bin;
    unsigned char* modBin;
    size_t binSize = 0;
    size_t modBinSize = 0;
    cl_program tmpProgram;

    bin = mwGetProgramBinary(integrationProgram, &binSize);
    if (!bin)
    {
        return MW_CL_ERROR;
    }

    modBin = getModifiedAMDBinary(bin, binSize, ap->number_streams, ci->di.calTarget, &modBinSize);
    free(bin);

    if (!modBin)
    {
        mw_printf("Error getting modified binary or IL source\n");
        return MW_CL_ERROR;
    }

    tmpProgram = mwCreateProgramFromBin(ci, modBin, modBinSize);
    free(modBin);
    if (!tmpProgram)
    {
        return MW_CL_ERROR;
    }

    /* Now that the IL program is OK replace the compiled one */
    clReleaseProgram(integrationProgram);
    integrationProgram = tmpProgram;

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

    if (!DOUBLEPREC || clr->forceNoILKernel)
        return CL_FALSE;

    /* Supporting these unused options with the IL kernel is too much work */
    if (ap->number_streams > maxILKernelStreams || ap->aux_bg_profile || ap->number_streams == 0)
        return CL_FALSE;

    /* Make sure an acceptable device */
    return (mwIsAMDGPUDevice(di) && isILKernelTarget(di) && mwPlatformSupportsAMDOfflineDevices(ci));
}

/* Return CL_TRUE on error */
static cl_bool setSummarizationWorkgroupSize(const CLInfo* ci)
{
    size_t maxGroupSize;
    size_t groupSize;
    size_t nextMultiple;
    cl_int err;

    err = clGetKernelWorkGroupInfo(_summarizationKernel,
                                   ci->dev,
                                   CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(maxGroupSize), &maxGroupSize,
                                   NULL);
    if (err != CL_SUCCESS)
    {
        return CL_TRUE;
    }

    if (maxGroupSize == 1)
    {
        /* Seems to be a problem on OS X CPU implementation */
        mw_printf("Workgroup size of 1 for summarization is not acceptable\n");
        return CL_TRUE;
    }

    /* OpenCL 1.1 has CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE
     * which we probably should be using if available  */

    nextMultiple = mwNextMultiple(64, ci->di.warpSize);
    if (nextMultiple <= maxGroupSize)
    {
        groupSize = nextMultiple;
    }
    else
    {
        groupSize = maxGroupSize;
    }

    if (groupSize > 128) /* Just in case */
    {
        groupSize = 128;
    }

    _summarizationWorkgroupSize = groupSize;

    return CL_FALSE;
}

cl_int setupSeparationCL(CLInfo* ci,
                         const AstronomyParameters* ap,
                         const IntegralArea* ias,
                         const CLRequest* clr)
{
    char* compileFlags;
    cl_bool useILKernel;
    cl_int err = MW_CL_ERROR;
    const char* kernSrc = (const char*) probabilities_kernel_cl;
    size_t kernSrcLen = probabilities_kernel_cl_len;

    const char* summarizationKernSrc = (const char*) summarization_kernel_cl;
    size_t summarizationKernSrcLen = summarization_kernel_cl_len;


    err = mwSetupCL(ci, clr);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error getting device and context");
        return err;
    }

    if (!separationCheckDevCapabilities(&ci->di))
    {
        return MW_CL_ERROR;
    }

    useILKernel = usingILKernelIsAcceptable(ci, ap, clr);
    compileFlags = getCompilerFlags(ci, ap, useILKernel);
    if (!compileFlags)
    {
        mw_printf("Failed to get CL compiler flags\n");
        return MW_CL_ERROR;
    }

    if (clr->verbose)
    {
        mw_printf("\nCompiler flags:\n%s\n\n", compileFlags);
    }

    integrationProgram = mwCreateProgramFromSrc(ci, 1, &kernSrc, &kernSrcLen, compileFlags);
    if (!integrationProgram)
    {
        mw_printf("Error creating integral program from source\n");
        err = MW_CL_ERROR;
        goto setup_exit;
    }

    summarizationProgram = mwCreateProgramFromSrc(ci, 1, &summarizationKernSrc, &summarizationKernSrcLen, compileFlags);
    if (!summarizationProgram)
    {
        mw_printf("Error creating summarization program from source\n");
        err = MW_CL_ERROR;
        goto setup_exit;
    }

    if (useILKernel)
    {
        mw_printf("Using AMD IL kernel\n");
        err = setProgramFromILKernel(ci, ap);
        if (err != CL_SUCCESS)
        {
            mw_printf("Failed to create IL kernel. Falling back to source kernel\n");
        }
    }

    if (err == CL_SUCCESS)
    {
        _separationKernel = mwCreateKernel(integrationProgram, "probabilities");
        _summarizationKernel = mwCreateKernel(summarizationProgram, "summarization");
        if (   !_separationKernel
            || !_summarizationKernel
            || setSummarizationWorkgroupSize(ci)
            || !separationCheckDevMemory(&ci->di, ap, ias))
        {
            err = MW_CL_ERROR;
        }
    }


setup_exit:
    free(compileFlags);

    return err;
}

