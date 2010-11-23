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

#include "milkyway_util.h"
#include "milkyway_cl.h"
#include "mw_cl.h"
#include "setup_cl.h"
#include "separation_cl_buffers.h"
#include "separation_cl_defs.h"
#include "separation_binaries.h"

#if SEPARATION_INLINE_KERNEL
  #include "integral_kernel.h"
#endif /* SEPARATION_INLINE_KERNEL */


/* Only sets the constant arguments, not the outputs which we double buffer */
cl_int separationSetKernelArgs(CLInfo* ci, SeparationCLMem* cm)
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

    if (err != CL_SUCCESS)
    {
        warn("Error setting kernel arguments: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

#if SEPARATION_INLINE_KERNEL

char* findKernelSrc()
{
    return integral_kernel_src;
}

void freeKernelSrc(char* src)
{
  #pragma unused(src)
}

#else

/* Reading from the file is more convenient for actually working on
 * it. Inlining is more useful for releasing when we don't want to
 * deal with the hassle of distributing more files. */
char* findKernelSrc()
{
    char* kernelSrc = NULL;
    kernelSrc = mwReadFile("../kernels/integrals.cl");
    if (!kernelSrc)
        warn("Failed to read kernel file\n");

    return kernelSrc;
}

void freeKernelSrc(char* src)
{
    free(src);
}

#endif


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
        warn("Total required device memory (%zu) > available (%zu)\n", totalMem, di->memSize);
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

  #if DOUBLEPREC
    if (!mwSupportsDoubles(di))
    {
        warn("Device doesn't support double precision\n");
        return CL_FALSE;
    }
  #endif

    return CL_TRUE;
}

/* Return flag for Nvidia compiler for maximum registers to use. */
static const char* getNvidiaRegCount(const DevInfo* di)
{
    cl_uint major, minor;
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
    const char* mathFlags = "-cl-mad-enable "
                            "-cl-no-signed-zeros "
                            "-cl-strict-aliasing "
                            "-cl-finite-math-only ";

    /* Build options used by milkyway_math stuff */
    const char* mathOptions = "-DUSE_CL_MATH_TYPES=0 "
                              "-DUSE_MAD=0 "
                              "-DUSE_FMA=0 ";

    /* Extra flags for different compilers */
    const char* nvidiaOptFlags = "-cl-nv-verbose ";
    const char* atiOptFlags    = "";

  #if DOUBLEPREC
    const char* precDefStr   = "-DDOUBLEPREC=1 ";
    const char* atiPrecStr   = "";
    const char* otherPrecStr = "";
  #else
    const char* precDefStr = "-DDOUBLEPREC=0 ";
    const char* atiPrecStr = "--single_precision_constant ";
    const char* clPrecStr  = "-cl-single-precision-constant ";
  #endif

    /* Constants compiled into kernel */
    const char* kernelDefStr = "-DNSTREAM=%u "
                               "-DFAST_H_PROB=%d "
                               "-DAUX_BG_PROFILE=%d "
                               "-DUSE_IMAGES=%d ";

    const char* includeStr = "-I%s/../include "
                             "-I%s/../../include "
                             "-I%s/../../milkyway/include ";

    /* Big enough. Also make sure to count for the extra characters of the format specifiers */
    char kernelDefBuf[sizeof(kernelDefStr) + 4 * 12 + 8];
    char precDefBuf[2 * sizeof(atiPrecStr) + sizeof(precDefStr)];

    size_t totalSize = 3 * sizeof(cwd) + (sizeof(includeStr) + 6)
                     + sizeof(mathFlags)
                     + sizeof(precDefBuf)
                     + sizeof(kernelDefBuf)
                     + sizeof(mathOptions)

                     + sizeof(extraFlags);

    if (snprintf(kernelDefBuf, sizeof(kernelDefBuf), kernelDefStr,
                 ap->number_streams,
                 ap->fast_h_prob,
                 ap->aux_bg_profile,
                 useImages) < 0)
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

    if (!getcwd(cwd, sizeof(cwd)))
    {
        perror("getcwd");
        return NULL;
    }

    if (snprintf(includeFlags, sizeof(includeFlags), includeStr, cwd, cwd, cwd) < 0)
    {
        warn("Failed to get include flags\n");
        return NULL;
    }

    compileFlags = mallocSafe(totalSize);
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
        return -1;
    }

    kernelSrc = findKernelSrc();
    if (!kernelSrc)
    {
        warn("Failed to read CL kernel source\n");
        return -1;
    }

    warn("\nCompiler flags:\n%s\n\n", compileFlags);
    err = mwSetProgramFromSrc(ci, "mu_sum_kernel", &kernelSrc, 1, compileFlags);

    freeKernelSrc(kernelSrc);
    free(compileFlags);

    if (err != CL_SUCCESS)
    {
        warn("Error creating program from source: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

