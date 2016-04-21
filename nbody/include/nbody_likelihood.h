/*
 * Copyright (c) 2010, 2011 Ben Willett
 * Copyright (c) 2010, 2011 Matthew Arsenault
 * Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute.
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

#ifndef _NBODY_LIKELIHOOD_H_
#define _NBODY_LIKELIHOOD_H_

#include "nbody_types.h"
#include "nbody.h"

#ifdef __cplusplus
extern "C" {
#endif

real nbSystemLikelihood(const NBodyState* st,
                     const NBodyHistogram* data,
                     const NBodyHistogram* histogram,
                     NBodyLikelihoodMethod method);

int nbGetLikelihoodInfo(const NBodyFlags* nbf, HistogramParams* hp, NBodyLikelihoodMethod* method);

real nbMatchHistogramFiles(const char* datHist, const char* matchHist);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_LIKELIHOOD_H_ */
