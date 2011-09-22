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

#ifndef _NBODY_CL_H_
#define _NBODY_CL_H_

#include "nbody.h"
#include "mw_cl.h"


#ifdef __cplusplus
extern "C" {
#endif



#define NBODY_MAXDEPTH 26

typedef struct
{
    cl_mem pos[3];
    cl_mem vel[3];
    cl_mem acc[3];
    cl_mem max[3];
    cl_mem min[3];
    cl_mem masses;
    cl_mem treeStatus;

    cl_mem start; /* TODO: We can reuse other buffers with this later to save memory */
    cl_mem count;
    cl_mem child;
    cl_mem sort;

    cl_mem critRadii; /* Used by the alternative cell opening criterion.
                         Unnecessary for BH86.
                         BH86 will be the fastest option since it won't need to load from this
                       */

    cl_mem debug;
} NBodyBuffers;


typedef struct
{
    size_t factors[6];
    size_t threads[6];
    double timings[6];        /* In a single iteration */
    double chunkTimings[6];   /* Average time per chunk */
    double kernelTimings[6];  /* Running totals */

    size_t global[6];
    size_t local[6];
} NBodyWorkSizes;



NBodyStatus runSystemCL(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf);


#ifdef __cplusplus
}
#endif

#endif /* _NBODY_CL_H_ */

