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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef _WIN32
  #include <direct.h>
#endif /* _WIN32 */


#if SEPARATION_INLINE_KERNEL
extern char* inlinedIntegralKernelSrc;
#else
static char* inlinedIntegralKernelSrc = NULL;
#endif /* SEPARATION_INLINE_KERNEL */


typedef struct
{
    cl_uint x;  /* Extra area */
    cl_uint n;  /* Number chunks */
} SizeSolution;

static size_t findGroupSize(const DevInfo* di)
{
    return di->devType == CL_DEVICE_TYPE_CPU ? 1 : 64;
}

static cl_uint findGroupsPerCU(const DevInfo* di)
{
    if (di->devType == CL_DEVICE_TYPE_CPU)
        return 1;

    if (di->vendorID == MW_NVIDIA)
        return cudaCoresPerComputeUnit(di);

    return 1; /* TODO: ATI, etc. */
}

static cl_uint nvidiaNumChunks(const IntegralArea* ia, const DevInfo* di, cl_uint nthread)
{
    cl_uint n = (ia->mu_steps * ia->r_steps) / (nthread * di->maxCompUnits);

    /* The 480 seems to be fine with half, although this doesn't help that much. */
    if (minComputeCapabilityCheck(di, 2, 0))
        n /= 2;

    n += 1;  /* ceil, avoid 0 in tiny cases like the small tests */

    return n;
}

static cl_uint chooseNumChunk(const IntegralArea* ia, const CLRequest* clr, const DevInfo* di, cl_uint nthread)
{
    /* If not being used for output, 1 has the least overhead */
    if (di->nonOutput || clr->nonResponsive)
        return 1;

    if (clr->numChunk != 0)   /* Use manual override */
        return clr->numChunk;

    if (di->vendorID == MW_NVIDIA)
        return nvidiaNumChunks(ia, di, nthread);

    return 1; /* FIXME: others */
}

static cl_bool checkSolution(cl_uint area, cl_uint block, cl_uint n, cl_uint x)
{
    return (area + x) % (block * n) == 0;
}

/* Find the extra area added */
static cl_bool findMinWorkN(SizeSolution* sol, cl_uint tolerance, cl_uint area, cl_uint block, cl_uint n, cl_uint nThread)
{
    cl_uint x = 0;
    cl_bool isSolution = CL_FALSE;

    while (!isSolution && x < tolerance)
    {
      /* If we found a solution, it's the minimum amount of work already.
         Don't care about ones with the same n and more extra work. */
        isSolution = checkSolution(area, block, n, x);
        if (isSolution)
        {
            sol->x = x;
            sol->n = n;
        }

        x += nThread; /* x will be a multiple of the local size */
    }

    return isSolution;
}

/* Finds number of chunks that will divide the area if all else fails */
static cl_uint findFallbackSolution(const IntegralArea* ia, cl_uint area, cl_uint desiredChunkNum, cl_uint nthread)
{
    cl_bool minCloser;
    cl_uint nMin, nMax;

    nMin = nMax = desiredChunkNum;

    while (nMin > 1 && (area % nMin != 0))
        --nMin;

    while (nMax < ia->r_steps && (area % nMax != 0))
        ++nMax;

    /* Choose the closer one, unless the max is the max possible, which doesn't work very well */
    minCloser = desiredChunkNum - nMin <= nMax - desiredChunkNum;
    return (minCloser || nMax == ia->r_steps) ? nMin : nMax;
}

static SizeSolution chooseLowerOrHigher(SizeSolution max, SizeSolution min, cl_uint desiredChunkNum)
{
    cl_bool minCloser;

    /* Pick the direction to go. First try choosing the closest one */
    minCloser = max.n - desiredChunkNum >= desiredChunkNum - min.n;
    if (minCloser)
        return min;

    /* If the max solution adds many more chunks, prefer the lower one anyway. */
    if (max.n >= 3 * min.n / 2)
    {
        warn("Next highest solution has many more, prefering smaller solution\n");
        return min;
    }

    return max;
}

/* Brute force possible solutions */
static SizeSolution findSolution(const IntegralArea* ia,
                                 cl_uint desiredChunkNum,
                                 cl_uint nThread,
                                 cl_uint groupSize,
                                 cl_uint nCU)
{
    cl_uint n, maxN;
    cl_uint area, block, tolerance;

    SizeSolution min, max;
    SizeSolution sol = { 0, 1 }; /* Minimum possible */

    cl_bool haveSolution = CL_FALSE;
    cl_bool minSolution = CL_FALSE;
    cl_bool maxSolution = CL_FALSE;

    min = max = sol;

    area = ia->r_steps * ia->mu_steps;
    tolerance = area / 40; /* Solutions can only add max. ~2.5% more work */
    block = nThread * nCU;

    warn("Block size = %u\n", block);

    if (nCU == 0 || nCU == 1 || desiredChunkNum == 0 || desiredChunkNum == 1)
        return sol;

    warn("Desired = %u\n", desiredChunkNum);
    /* Find the closest one with a smaller n */
    n = desiredChunkNum;
    while (!minSolution && n > 1)
    {
        minSolution = findMinWorkN(&min, tolerance, area, block, n, groupSize);
        warn("Min sol: %u %u\n", min.n, min.x);

        --n;
    }

    /* Find the closest one with a larger n */
    n = desiredChunkNum;
    maxN = ia->r_steps / 2;
    while (!maxSolution && n < maxN)
    {
        maxSolution = findMinWorkN(&max, tolerance, area, block, n, groupSize);
        ++n;
    }

    if (minSolution)
        warn("Lower n solution: n = %u, x = %u\n", min.n, min.x);

    if (maxSolution)
        warn("Higher n solution: n = %u, x = %u\n", max.n, max.x);

    haveSolution = minSolution || maxSolution;
    if (!haveSolution)
    {
        /* This really shouldn't happen, but if it does we give up on
           adding extra work. Find the closest n which will divide the
           total area */
        sol.n = findFallbackSolution(ia, area, desiredChunkNum, nThread);
        sol.x = 0;
        warn("Didn't find a solution. Using fallback solution n = %u, x = %u\n", sol.n, sol.x);
        /* TODO: Try a bigger tolerance */
    }
    else
    {
        sol = chooseLowerOrHigher(max, min, desiredChunkNum);
    }

    return sol;
}

static SizeSolution chooseSolution(const SizeSolution* allSol,
                                   cl_uint numSolution,
                                   cl_uint desiredNumChunk)
{
    cl_uint i = 0;
    SizeSolution sol;

    sol.x = 0;
    sol.n = 1;

    if (!allSol || numSolution == 0) /* No solutions to look for */
        return sol;

    /* Use minimum number of chunks */
    if (desiredNumChunk == 1 || desiredNumChunk == 0)
        return sol;

    /* Find where we cross */
    while (allSol[i].n < desiredNumChunk && i < numSolution)
        ++i;

    if (i == numSolution)
    {
        warn("Desired number of chunks went beyond solution range. "
             "Choosing highest solution\n");
        sol = allSol[numSolution - 1];
    }
    else
    {
        sol = allSol[i];
    }

    return sol;
}

/* Returns CL_TRUE on error */
cl_bool findGoodRunSizes(RunSizes* sizes,
                         const CLInfo* ci,
                         const DevInfo* di,
                         const IntegralArea* ia,
                         const CLRequest* clr)
{
    /* To avoid lots of idle compute units near the end of each iteration, it helps to round the area up so that the total area is evenly divisible by the number of threads times the number of compute units when chunked.

       We can ensure this by finding solutions to this this equation:

       There exists some integers m and x such that:

       m == (mu_steps * r_steps + x) / (n_chunk * n_thread * n_cu)

       where n_chunk = number of chunks the iteration is divided into
             n_thread = number of threads used for local size, i.e 64
             n_cu = number of compute units of the GPU

        Good numbers of chunks have solutions for m and n, of which we only care what x is. Otherwise there will be more idle time.

        groupSize = nthreads
     */

    SizeSolution solution;
    cl_uint desiredNumChunk;
    WGInfo wgi;
    size_t localSize;
    cl_int err;
    cl_uint groupSize, groupsPerCU, threadsPerCU;

    err = mwGetWorkGroupInfo(ci, &wgi);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get work group info: %s\n", showCLInt(err));
        return CL_TRUE;
    }
    else
        mwPrintWorkGroupInfo(&wgi);

    groupSize   = findGroupSize(di);
    groupsPerCU = findGroupsPerCU(di);
    threadsPerCU = groupSize * groupsPerCU;

    warn("Group size = %u, per CU = %u, threads per CU = %u\n",
         groupSize, groupsPerCU, threadsPerCU);

    desiredNumChunk = chooseNumChunk(ia, clr, di, threadsPerCU);
    solution = findSolution(ia,
                            desiredNumChunk,
                            threadsPerCU,
                            groupSize,
                            di->maxCompUnits);
    warn("Using solution: n = %u, x = %u\n", solution.n, solution.x);

    sizes->groupSize     = groupSize;
    sizes->area          = ia->mu_steps * ia->r_steps;
    sizes->numChunks     = solution.n;
    sizes->extra         = solution.x;
    sizes->effectiveArea = sizes->area + sizes->extra;
    sizes->chunkSize     = sizes->effectiveArea / sizes->numChunks;

    sizes->global[0] = sizes->chunkSize;;
    sizes->global[1] = 1;

    sizes->local[0] = sizes->groupSize;
    sizes->local[1] = 1;

    warn("Range:          { nu_steps = %u, mu_steps = %u, r_steps = %u }\n"
         "Iteration area: "ZU"\n"
         "Chunk estimate: %u\n"
         "Num chunks:     "ZU"\n"
         "Added area:     %u\n"
         "Effective area: "ZU"\n",
         ia->nu_steps, ia->mu_steps, ia->r_steps,
         sizes->area,
         desiredNumChunk,
         sizes->numChunks,
         sizes->extra,
         sizes->effectiveArea);

    /* Check for error in solution just in case */
    if (sizes->effectiveArea % sizes->chunkSize != 0)
    {
        warn("Effective area ("ZU") not divisible by chunk size ("ZU")\n",
             sizes->effectiveArea, sizes->chunkSize);
        return CL_TRUE;
    }

    /* TODO: also check against CL_DEVICE_MAX_WORK_ITEM_SIZES */
    if (sizes->global[0] % sizes->local[0] || sizes->global[1] % sizes->local[1])
    {
        warn("Global dimensions not divisible by local\n");
        return CL_TRUE;
    }

    localSize = sizes->local[0] * sizes->local[1];
    if (localSize > wgi.wgs)
    {
        warn("Local size ("ZU") > maximum work group size ("ZU")\n", localSize, wgi.wgs);
        return CL_TRUE;
    }

    return CL_FALSE;
}


/* Only sets the constant arguments, not the outputs which we double buffer */
cl_int separationSetKernelArgs(CLInfo* ci, SeparationCLMem* cm, const RunSizes* runSizes)
{
    cl_int err = CL_SUCCESS;

    /* Set output buffer arguments */
    err |= clSetKernelArg(ci->kern, 0, sizeof(cl_mem), &cm->outMu);
    err |= clSetKernelArg(ci->kern, 1, sizeof(cl_mem), &cm->outProbs);

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
        warn("Error setting kernel arguments: %s\n", showCLInt(err));
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
cl_bool separationCheckDevCapabilities(const DevInfo* di, const SeparationSizes* sizes)
{
    size_t totalOut;
    size_t totalConstBuf;
    size_t totalGlobalConst;
    size_t totalMem;

    totalOut = 2 * sizes->outMu + 2 * sizes->outProbs; /* 2 buffers for double buffering */
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

    if (sizes->outMu > di->maxMemAlloc || sizes->outProbs > di->maxMemAlloc)
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
    char includeFlags[4096] = "";

    /* Math options for CL compiler */
    const char mathFlags[] = "-cl-mad-enable "
                             "-cl-no-signed-zeros "
                             "-cl-strict-aliasing "
                             "-cl-finite-math-only ";

    /* Build options used by milkyway_math stuff */
    const char mathOptions[] = "-DUSE_CL_MATH_TYPES=0 "
                               "-DUSE_MAD=1 "
                               "-DUSE_FMA=0 ";

    /* Extra flags for different compilers */
    const char nvidiaOptFlags[] = "-cl-nv-verbose ";
    const char atiOptFlags[]    = "";

  #if DOUBLEPREC
    const char precDefStr[]   = "-DDOUBLEPREC=1 ";
    const char atiPrecStr[]   = "";
    const char otherPrecStr[] = "";
  #else
    const char precDefStr[] = "-DDOUBLEPREC=0 ";
    const char atiPrecStr[] = "--single_precision_constant ";
    const char clPrecStr[]  = "-cl-single-precision-constant ";
  #endif

    /* Constants compiled into kernel. We need to define
     * MILKYWAY_MATH_COMPILATION since the header inclusion protection
     * doesn't really matter and doesn't work when everything is
     * dumped together. */
    const char kernelDefStr[] = "-DMILKYWAY_MATH_COMPILATION "
                                "-DNSTREAM=%u "
                                "-DFAST_H_PROB=%d "
                                "-DAUX_BG_PROFILE=%d "
                                "-DUSE_IMAGES=%d "
                                "-DI_DONT_KNOW_WHY_THIS_DOESNT_WORK_HERE=%d ";

    const char includeStr[] = "-I%s "
                              "-I%s/../include "
                              "-I%s/../../include "
                              "-I%s/../../milkyway/include ";

    /* Big enough. Also make sure to count for the extra characters of the format specifiers */
    char kernelDefBuf[sizeof(kernelDefStr) + 5 * 12 + 8];
    char precDefBuf[2 * sizeof(atiPrecStr) + sizeof(precDefStr) + 1];

    size_t totalSize = 4 * sizeof(cwd) + (sizeof(includeStr) + 8)
                     + sizeof(mathFlags)
                     + sizeof(precDefBuf)
                     + sizeof(kernelDefBuf)
                     + sizeof(mathOptions)
                     + sizeof(extraFlags);

    cl_bool isFermi = minComputeCapabilityCheck(di, 2, 0);

    if (snprintf(kernelDefBuf, sizeof(kernelDefBuf), kernelDefStr,
                 ap->number_streams,
                 ap->fast_h_prob,
                 ap->aux_bg_profile,
                 useImages,
                 isFermi) < 0)
    {
        warn("Error getting kernel constant definitions\n");
        return NULL;
    }

    /* Always use this flag */
    strncpy(precDefBuf, precDefStr, sizeof(precDefBuf));

  #if !DOUBLEPREC
    /* The ATI compiler rejects the one you're supposed to use, in
     * favor of a totally undocumented flag. */
    strncat(precDefBuf,
            di->vendorID != MW_AMD_ATI ? clPrecStr : atiPrecStr,
            2 * sizeof(atiPrecStr));
  #endif /* !DOUBLEPREC */

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
        strncat(extraFlags, atiOptFlags, sizeof(extraFlags));
    else
        warn("Unknown vendor ID: 0x%x\n", di->vendorID);

    if (!inlinedIntegralKernelSrc)
    {
        if (!getcwd(cwd, sizeof(cwd)))
        {
            warn("Failed to get header directory\n");
            return NULL;
        }

        if (snprintf(includeFlags, sizeof(includeFlags), includeStr, cwd, cwd, cwd, cwd) < 0)
        {
            warn("Failed to get include flags\n");
            return NULL;
        }
    }

    compileFlags = mwMalloc(totalSize);
    if (snprintf(compileFlags, totalSize, "%s%s%s%s%s%s ",
                 includeFlags,
                 mathFlags,
                 mathOptions,
                 extraFlags,
                 precDefBuf,
                 kernelDefBuf) < 0)
    {
        warn("Failed to get compile flags\n");
        free(compileFlags);
        return NULL;
    }

    return compileFlags;
}

/* Bad estimate */
cl_double estimateWUFLOPsPerIter(const AstronomyParameters* ap, const IntegralArea* ia)
{
    cl_ulong perItem, perIter; /* Needs 64 bit int */

    perItem = 4
            + 28 * ap->convolve
            + 4  * ap->number_streams
            + 51 * ap->convolve * ap->number_streams;

    perIter = perItem * ia->mu_steps * ia->r_steps;

    return (cl_double) perIter;
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
                         DevInfo* di,
                         const AstronomyParameters* ap,
                         const StreamConstants* sc,
                         const StreamGauss sg,
                         const CLRequest* clr,
                         cl_bool useImages)
{
    cl_int err;
    char* compileFlags;
    char* kernelSrc;

    err = mwSetupCL(ci, di, clr);
    if (err != CL_SUCCESS)
    {
        warn("Error getting device and context: %s\n", showCLInt(err));
        return err;
    }

    compileFlags = getCompilerFlags(ap, di, useImages);
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
        warn("Error creating program from source: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

