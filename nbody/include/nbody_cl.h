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


cl_bool setWorkSizes(NBodyWorkSizes* ws, const DevInfo* di);
cl_bool setThreadCounts(NBodyWorkSizes* ws, const DevInfo* di);
cl_int nbodyLoadKernels(const NBodyCtx* ctx, NBodyState* st);
cl_int nbodyCreateKernels(NBodyState* st);
cl_bool nbodyCheckDevCapabilities(const DevInfo* di, const NBodyCtx* ctx, NBodyState* st);

cl_int nbodySetInitialTreeStatus(NBodyState* st);
cl_int nbodyCreateBuffers(const NBodyCtx* ctx, NBodyState* st);
cl_int nbodyReleaseBuffers(NBodyState* st);
cl_int nbodyReleaseKernels(NBodyState* st);

cl_int nbodySetAllKernelArguments(NBodyState* st);

NBodyStatus runSystemCL(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf);

cl_int nbodyMarshalBodies(NBodyState* st, cl_bool marshalIn);


#ifdef __cplusplus
}
#endif

#endif /* _NBODY_CL_H_ */

