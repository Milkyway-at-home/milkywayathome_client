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

#ifndef _SETUP_CL_H_
#define _SETUP_CL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_cl.h"
#include "mw_cl.h"
#include "separation_types.h"


typedef struct
{
    size_t outBg;
    size_t outStreams;

    size_t ap;        /* Constants */
    size_t sc;
    size_t ia;
    size_t rc;
    size_t rPts;
    size_t sg_dx;
    size_t lbts;
} SeparationSizes;

typedef struct
{

    size_t* chunkBorders;

    size_t global[2];

    size_t local[2];
    size_t groupSize;
    size_t nChunkEstimate;  /* Target number of chunks to use */
    size_t nChunk;     /* Number of chunks to divide each iteration into */
    size_t extra;      /* Extra area added */
    size_t area;
    size_t effectiveArea;
    size_t chunkSize;   /* effectiveArea / numChunks */
} RunSizes;

/* The various buffers needed by the integrate function. */
typedef struct
{
    cl_mem outBg;
    cl_mem outStreams;  /* stream_probs * V * reff_xr_rp3 */

    /* constant, read only buffers */
    cl_mem ap;
    cl_mem ia;
    cl_mem sc;        /* Stream Constants */
    cl_mem rc;        /* r constants */
    cl_mem rPts;
    cl_mem sg_dx;
    cl_mem lbts;
} SeparationCLMem;

#define EMPTY_SEPARATION_CL_MEM { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }


cl_int setupSeparationCL(CLInfo* ci,
                         const AstronomyParameters* ap,
                         const IntegralArea* ias,
                         const CLRequest* clr,
                         cl_int* useImages);

cl_bool separationCheckDevCapabilities(const DevInfo* di, const AstronomyParameters* ap, const IntegralArea* ias);

cl_int separationSetKernelArgs(CLInfo* ci, SeparationCLMem* cm, const RunSizes* runSizes);

void freeRunSizes(RunSizes* sizes);

cl_bool findRunSizes(RunSizes* sizes,
                     const CLInfo* ci,
                     const DevInfo* di,
                     const IntegralArea* ia,
                     const CLRequest* clr);

cl_double cudaEstimateIterTime(const DevInfo* di, cl_double flopsPerIter, cl_double flops);

#ifdef __cplusplus
}
#endif

#endif /* _SETUP_CL_H_ */

