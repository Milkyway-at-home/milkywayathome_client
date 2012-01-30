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


cl_bool nbSetWorkSizes(NBodyWorkSizes* ws, const DevInfo* di);
cl_bool nbSetThreadCounts(NBodyWorkSizes* ws, const DevInfo* di, const NBodyCtx* ctx);
cl_int nbFindEffectiveNBody(const NBodyWorkSizes* ws, cl_bool exact, cl_int nbody);

cl_bool nbLoadKernels(const NBodyCtx* ctx, NBodyState* st);
cl_bool nbCheckDevCapabilities(const DevInfo* di, const NBodyCtx* ctx, cl_uint nbody);

cl_int nbSetInitialTreeStatus(NBodyState* st);
cl_int nbCreateBuffers(const NBodyCtx* ctx, NBodyState* st);
cl_int nbReleaseBuffers(NBodyState* st);
cl_int nbReleaseKernels(NBodyState* st);

cl_int nbSetAllKernelArguments(NBodyState* st);

cl_int nbMarshalBodies(NBodyState* st, cl_bool marshalIn);
void nbPrintKernelTimings(const NBodyState* st);


NBodyStatus nbStepSystemCL(const NBodyCtx* ctx, NBodyState* st);
NBodyStatus nbRunSystemCL(const NBodyCtx* ctx, NBodyState* st);


#ifdef __cplusplus
}
#endif

#endif /* _NBODY_CL_H_ */

