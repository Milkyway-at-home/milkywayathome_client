/*
 *  Copyright (c) 2011 Matthew Arsenault
 *  Copyright (c) 2011 Rensselaer Polytechnic Institute.
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

#ifndef _NBODY_UTIL_H_
#define _NBODY_UTIL_H_

#include "nbody_types.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define _nbValidPositionItem(x) (!isinf(x) && !isnan(x))
#define nbPositionValid(r) (_nbValidPositionItem(r.x) && _nbValidPositionItem(r.y) && _nbValidPositionItem(r.z))

real nbCorrectTimestep(real timeEvolve, real dt);
mwvector nbCenterOfMass(const NBodyState* st);
mwvector nbCenterOfMom(const NBodyState* st);

real nbEstimateNumberFlops(const NBodyCtx* ctx, int nbody);
real nbEstimateTime(const NBodyCtx* ctx, int nbody, real flops);

void nbReportTreeIncest(const NBodyCtx* ctx, NBodyState* st);

#ifdef _OPENMP
#define nbGetMaxThreads() omp_get_max_threads()
#else
#define nbGetMaxThreads() (1)
#endif

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_UTIL_H_ */

