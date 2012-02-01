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

#ifndef _SETUP_CL_H_
#define _SETUP_CL_H_

#include "milkyway_cl.h"
#include "separation_types.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef struct
{
    size_t summarizationBufs[2];
    size_t outBg;
    size_t outStreams;

    size_t ap;        /* Constants */
    size_t sc;
    size_t ia;
    size_t rc;
    size_t rPts;
    size_t sg_dx;
    size_t lTrig;
    size_t bSin;

    cl_int nStream;
} SeparationSizes;

typedef struct
{
    size_t global[1];
    size_t local[1];
    size_t groupSize;
    size_t nChunkEstimate;  /* Target number of chunks to use */
    size_t nChunk;          /* Number of chunks to divide each iteration into */
    cl_uint extra;          /* Extra area added */
    cl_uint initialWait;    /* If manually polling for kernel completion how long to initially wait */

    cl_uint r, mu, nu;
    cl_ulong area;
    cl_ulong effectiveArea;
    size_t chunkSize;        /* effectiveArea / numChunks */
} RunSizes;

/* The various buffers needed by the integrate function. */
typedef struct
{
    cl_mem summarizationBufs[2];
    cl_mem outBg;
    cl_mem outStreams;  /* stream_probs */

    cl_mem rc;          /* r constants */
    cl_mem rPts;
    cl_mem lTrig;
    cl_mem bSin;

    /* constant, read only buffers */
    cl_mem ap;
    cl_mem sc;        /* Stream Constants */
    cl_mem sg_dx;
} SeparationCLMem;

#define EMPTY_SEPARATION_CL_MEM { { NULL, NULL }, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }

cl_int setupSeparationCL(CLInfo* ci,
                         const AstronomyParameters* ap,
                         const IntegralArea* ias,
                         const CLRequest* clr);


cl_int separationSetKernelArgs(SeparationCLMem* cm, const RunSizes* runSizes);

cl_bool findRunSizes(RunSizes* sizes,
                     const CLInfo* ci,
                     const DevInfo* di,
                     const AstronomyParameters* ap,
                     const IntegralArea* ia,
                     const CLRequest* clr);

cl_double cudaEstimateIterTime(const DevInfo* di, cl_double flopsPerIter, cl_double flops);

extern cl_kernel _separationKernel;
extern cl_kernel _summarizationKernel;
extern size_t _summarizationWorkgroupSize;

cl_int releaseSeparationKernel(void);

#ifdef __cplusplus
}
#endif

#endif /* _SETUP_CL_H_ */

