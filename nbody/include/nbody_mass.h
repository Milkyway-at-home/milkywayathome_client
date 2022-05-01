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

real_0 probability_match(int n, real_0 ktmp, real_0 pobs);

real_0 GammaFunc(const real_0 z);

real_0 IncompleteGammaFunc(real_0 a, real_0 x);

real nbCostComponent(const NBodyHistogram* data, const NBodyHistogram* histogram);

real calc_vLOS(const mwvector* v, const mwvector* p, real_0 sunGCdist);

real calc_distance(const mwvector* p, real_0 sunGCdist);

mwvector LBtoY(real* lambda, real* beta);

real ln_FB5_dist(mwvector* x, mwvector* gamma1, mwvector* gamma2, mwvector* gamma3, real* kappa, real* beta);

real add_logspace(real* r1, real* r2);

real getBodyBinFrac(const NBodyCtx* ctx, const HistogramParams* hp, const Body* p, int HistIndex);

void nbCalcDisp(NBodyHistogram* histogram, mwbool initial, real_0 correction_factor);

void nbRemoveOutliers(const NBodyState* st, NBodyHistogram* histogram, real_0 * use_body, real * var, real_0 sigma_cutoff, int histBins, mwbool contBins, real * bodyFraction);

real nbLikelihood(const NBodyHistogram* data, const NBodyHistogram* histogram);

#ifdef __cplusplus
}
#endif

#endif
