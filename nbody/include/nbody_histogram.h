/*
 * Copyright (c) 2010, 2011 Ben Willett
 * Copyright (c) 2010, 2011 Matthew Arsenault
 * Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute.
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

#ifndef _NBODY_HISTOGRAM_H_
#define _NBODY_HISTOGRAM_H_

#include "nbody_types.h"
#include "nbody.h"


#ifdef __cplusplus
extern "C" {
#endif

AllHistograms* nbReadHistogram(const char* histogramFile);

AllHistograms* nbCreateHistogram(const NBodyCtx* ctx, const NBodyState* st, const HistogramParams* hp);

void nbPrintHistogram(FILE* f, const AllHistograms* histogram);

void nbWriteHistogram(const char* histoutFileName,
                      const NBodyCtx* ctx,
                      const NBodyState* st,
                      const AllHistograms* histogram);

real nbCorrectRenormalizedInHistogram(const NBodyHistogram* histogram, const NBodyHistogram* data);

real nbNormalizedHistogramError(unsigned int n, real total);

unsigned int nbCorrectTotalNumberInHistogram(const NBodyHistogram* histogram, /* Generated histogram */
					     const NBodyHistogram* data);      /* Data histogram */


#ifdef __cplusplus
}
#endif

#endif /* _NBODY_HISTOGRAM_H_ */

