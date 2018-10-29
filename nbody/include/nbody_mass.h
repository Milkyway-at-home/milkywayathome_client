/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
 * Copyright (c) 2016-2018 Siddhartha Shelton
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

#ifndef _NBODY_MASS_H_
#define _NBODY_MASS_H_

#include "nbody_types.h"
#include "milkyway_math.h"

#ifdef __cplusplus
extern "C" {
#endif

real probability_match(int n, real k, real pobs);

real GammaFunc(const real z);

real IncompleteGammaFunc(real a, real x);

real nbCostComponent(const NBodyHistogram* data, const NBodyHistogram* histogram);

real calc_vLOS(const mwvector v, const mwvector p, real sunGCdist);

real calc_distance(const mwvector p, real sunGCdist);	

void nbCalcVelDisp(NBodyHistogram* histogram, mwbool initial, real correction_factor);
void nbCalcBetaDisp(NBodyHistogram* histogram, mwbool initial, real correction_factor);

void nbRemoveVelOutliers(const NBodyState* st, NBodyHistogram* histogram, real * use_velbody, real * vlos, real sigma_cutoff, real sunGCdist);
void nbRemoveBetaOutliers(const NBodyState* st, NBodyHistogram* histogram, real * use_betabody, real * betas, real sigma_cutoff, real sunGCdist);

real nbVelocityDispersion(const NBodyHistogram* data, const NBodyHistogram* histogram);
real nbBetaDispersion(const NBodyHistogram* data, const NBodyHistogram* histogram);
real nbBetaComponent(const NBodyHistogram* data, const NBodyHistogram* histogram);
real nbLOSVelocityComponent(const NBodyHistogram* data, const NBodyHistogram* histogram);
real nbDistanceComponent(const NBodyHistogram* data, const NBodyHistogram* histogram);

#ifdef __cplusplus
}
#endif

#endif
