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
#include "separation_cl_defs.h"
#include "separation_binaries.h"


#include <assert.h>

#ifdef _WIN32
  #include <direct.h>
#endif /* _WIN32 */


#if SEPARATION_INLINE_KERNEL
extern char* inlinedIntegralKernelSrc;
#else
static char* inlinedIntegralKernelSrc = NULL;
#endif /* SEPARATION_INLINE_KERNEL */


/* Find an optimum block size to use and use that as the basis for the chunk estimate. Smaller number of chunks on a non-block size can be quite a bit slower. */
static cl_uint nvidiaNumChunks(const IntegralArea* ia, const DevInfo* di)
{
    cl_uint n;

    n = (ia->mu_steps * ia->r_steps) / mwBlockSize(di);

    /* The 480 seems to be fine with half, although this doesn't help that much. */
    if (minComputeCapabilityCheck(di, 2, 0))
        n /= 2;

    n += 1;  /* ceil, avoid 0 in tiny cases like the small tests */

    return (cl_uint) n;
}

static cl_uint chooseNumChunk(const IntegralArea* ia, const CLRequest* clr, const DevInfo* di)
{
    /* If not being used for output, 1 has the least overhead */
    if (di->nonOutput || clr->nonResponsive)
        return 1;

    if (di->vendorID == MW_NVIDIA)
        return nvidiaNumChunks(ia, di);

    return 1; /* FIXME: others */
}

void freeRunSizes(RunSizes* sizes)
{
    mwFreeA(sizes->chunkBorders);
}

static void printRunSizes(const RunSizes* sizes, const IntegralArea* ia, cl_bool verbose)
{
    size_t i;

    warn("Range:          { nu_steps = %u, mu_steps = %u, r_steps = %u }\n"
         "Iteration area: "ZU"\n"
         "Chunk estimate: "ZU"\n"
         "Num chunks:     "ZU"\n"
         "Added area:     "ZU"\n"
         "Effective area: "ZU"\n",
         ia->nu_steps, ia->mu_steps, ia->r_steps,
         sizes->area,
         sizes->nChunkEstimate,
         sizes->nChunk,
         sizes->extra,
         sizes->effectiveArea);

    if (!verbose)
        return;

    warn("Using "ZU" chunks with size(s): ", sizes->nChunk);
    for (i = 0; i < sizes->nChunk; ++i)
    {
        warn(" "ZU" ", sizes->chunkBorders[i + 1] - sizes->chunkBorders[i]);
    }
    warn("\n");

}

/* Returns CL_TRUE on error */
cl_bool findRunSizes(RunSizes* sizes,
                     const CLInfo* ci,
                     const DevInfo* di,
                     const IntegralArea* ia,
                     const CLRequest* clr)
{
    WGInfo wgi;
    cl_int err;
    cl_uint groupSize;
    size_t i, nMod;
    size_t sum = 0;

    err = mwGetWorkGroupInfo(ci, &wgi);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to get work group info", err);
        return CL_TRUE;
    }

    mwPrintWorkGroupInfo(&wgi);
    groupSize = mwFindGroupSize(di);

    sizes->local[0] = groupSize;
    sizes->local[1] = 1;

    sizes->nChunkEstimate = chooseNumChunk(ia, clr, di);
    sizes->nChunk = sizes->nChunkEstimate;

    /* Best for performance.
       Be a bit more flexible if chunking.
       Using the whole block size seems a bit better when attacking everything at once.
    */
    nMod = sizes->nChunk == 1 ? mwBlockSize(di) : groupSize * di->maxCompUnits;

    sizes->area = ia->r_steps * ia->mu_steps;
    sizes->effectiveArea = nMod * mwDivRoundup(sizes->area, nMod);
    sizes->extra = sizes->effectiveArea - sizes->area;

    warn("Keeping chunk boundaries as multiples of "ZU"\n", nMod);

    if (sizes->effectiveArea / sizes->nChunk < nMod)
    {
        sizes->nChunk = sizes->effectiveArea / nMod;
        warn("Warning: Estimated number of chunks ("ZU") too large. Using "ZU"\n", sizes->nChunkEstimate, sizes->nChunk);
    }

    if (nMod == 1)
    {
        /* When nMod = 1 we need to avoid losing pieces at the
         * beginning and end; the normal method doesn't quite work. */
        while (sizes->nChunk * (sizes->effectiveArea / sizes->nChunk) < sizes->effectiveArea)
        {
            sizes->nChunk++;
        }
        warn("Need to use "ZU" chunks to cover area using multiples of 1\n", sizes->nChunk);
    }

    sizes->chunkBorders = mwCallocA((sizes->nChunk + 1), sizeof(size_t));
    for (i = 0; i <= sizes->nChunk; ++i)
    {
        if (nMod == 1)
        {
            /* Avoid losing out the 0 border */
            sizes->chunkBorders[i] = i * (sizes->effectiveArea / sizes->nChunk);
        }
        else
        {
            sizes->chunkBorders[i] = (i * sizes->effectiveArea + sizes->nChunk) / (sizes->nChunk * nMod);
            sizes->chunkBorders[i] *= nMod;
        }

        if (sizes->chunkBorders[i] > sizes->effectiveArea)
            sizes->chunkBorders[i] = sizes->effectiveArea;

        if (i > 0)
        {
            sum += sizes->chunkBorders[i] - sizes->chunkBorders[i - 1];
            assert(sizes->chunkBorders[i] - sizes->chunkBorders[i - 1] > 0);
        }
    }

    printRunSizes(sizes, ia, clr->verbose);

    if (sum != sizes->effectiveArea)  /* Assert that the divisions aren't broken */
    {
        warn("Chunk total does not match: "ZU" != "ZU"\n", sum, sizes->effectiveArea);
        free(sizes->chunkBorders);
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

    /* The constant arguments */
    err |= clSetKernelArg(ci->kern, 2, sizeof(cl_mem), &cm->ap);
    err |= clSetKernelArg(ci->kern, 3, sizeof(cl_mem), &cm->ia);
    err |= clSetKernelArg(ci->kern, 4, sizeof(cl_mem), &cm->sc);
    err |= clSetKernelArg(ci->kern, 5, sizeof(cl_mem), &cm->rc);
    err |= clSetKernelArg(ci->kern, 6, sizeof(cl_mem), &cm->sg_dx);
    err |= clSetKernelArg(ci->kern, 7, sizeof(cl_mem), &cm->rPts);
    err |= clSetKernelArg(ci->kern, 8, sizeof(cl_mem), &cm->lbts);
    err |= clSetKernelArg(ci->kern, 9, sizeof(cl_uint), &runSizes->extra);

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
    totalConstBuf = sizes->ap + sizes->ia + sizes->sc + sizes->rc + sizes->sg_dx;
    totalGlobalConst = sizes->lbts + sizes->rPts;

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

    if (sizes->lbts > di->maxMemAlloc || sizes->rPts > di->maxMemAlloc)
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
    cl_uint i;
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

/* Return flag for Nvidia compiler for maximum registers to use. */
static const char* getNvidiaRegCount(const DevInfo* di)
{
    const char* regCount32 = "-cl-nv-maxrregcount=32 ";
    const char* regDefault = "";

    /* On the 285 GTX, max. 32 seems to help. Trying that on the 480
       makes it quite a bit slower. */

    if (computeCapabilityIs(di, 1, 3)) /* 1.3 == GT200 */
    {
        warn("Found a compute capability 1.3 device. Using %s\n", regCount32);
        return regCount32;
    }

    /* Higher or other is Fermi or unknown, */
    return regDefault;
}

/* Get string of options to pass to the CL compiler. */
static char* getCompilerFlags(const AstronomyParameters* ap, const DevInfo* di, cl_bool useImages)
{
    char* compileFlags = NULL;
    char cwd[1024] = "";
    char extraFlags[1024] = "";
    char includeFlags[1024] = "";
    char precBuf[1024] = "";
    char kernelDefBuf[4096] = "";

    /* Math options for CL compiler */
    const char mathFlags[] = "-cl-mad-enable "
                             "-cl-no-signed-zeros "
                             "-cl-strict-aliasing "
                             "-cl-finite-math-only ";

    /* Extra flags for different compilers */
    const char nvidiaOptFlags[] = "-cl-nv-verbose ";
    const char atiOptFlags[]    = "";

    const char kernelDefStr[] = "-DFAST_H_PROB=%d "
                                "-DAUX_BG_PROFILE=%d "

                                "-DNSTREAM=%u "
                                "-DCONVOLVE=%u "

                                "-DR0=%.15f "
                                "-DSUN_R0=%.15f "
                                "-DQ_INV_SQR=%.15f "
                                "-DUSE_IMAGES=%d ";

    const char includeStr[] = "-I%s "
                              "-I%s/../include ";

    if (snprintf(kernelDefBuf, sizeof(kernelDefBuf), kernelDefStr,
                 ap->fast_h_prob,
                 ap->aux_bg_profile,

                 ap->number_streams,
                 ap->convolve,

                 ap->r0,
                 ap->sun_r0,
                 ap->q_inv_sqr,

                 useImages) < 0)
    {
        warn("Error getting kernel constant definitions\n");
        return NULL;
    }

    snprintf(precBuf, sizeof(precBuf), "-DDOUBLEPREC=%d %s ",
             DOUBLEPREC, DOUBLEPREC ? "" : "-cl-single-precision-constant");

    /* FIXME: Device vendor not necessarily the platform vendor */
    if (di->vendorID == MW_NVIDIA)
    {
        if (snprintf(extraFlags, sizeof(extraFlags),
                     "%s %s ",
                     nvidiaOptFlags, getNvidiaRegCount(di)) < 0)
        {
            warn("Error getting extra Nvidia flags\n");
            return NULL;
        }
    }
    else if (di->vendorID == MW_AMD_ATI)
    {
        if (snprintf(extraFlags, sizeof(extraFlags), "%s", atiOptFlags) < 0)
        {
            warn("Error getting extra ATI flags\n");
            return NULL;
        }
    }
    else
        warn("Unknown vendor ID: 0x%x\n", di->vendorID);

    if (!inlinedIntegralKernelSrc)
    {
        if (!getcwd(cwd, sizeof(cwd)))
        {
            warn("Failed to get header directory\n");
            return NULL;
        }

        if (snprintf(includeFlags, sizeof(includeFlags), includeStr, cwd, cwd) < 0)
        {
            warn("Failed to get include flags\n");
            return NULL;
        }
    }

    if (asprintf(&compileFlags, "%s%s%s%s%s ",
                 includeFlags,
                 mathFlags,
                 extraFlags,
                 precBuf,
                 kernelDefBuf) < 0)
    {
        warn("Failed to get compile flags\n");
        return NULL;
    }

    return compileFlags;
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
                         const CLRequest* clr,
                         cl_int* useImages)
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

    *useImages = *useImages && ci->di.imgSupport;

    compileFlags = getCompilerFlags(ap, &ci->di, *useImages);
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
    err = mwSetProgramFromSrc(ci, "mu_sum_kernel", (const char**) &kernelSrc, 1, compileFlags);

    freeKernelSrc(kernelSrc);
    free(compileFlags);

    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating program from source", err);
        return err;
    }

    return CL_SUCCESS;
}

