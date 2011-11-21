/* Copyright 2010 Ben Willett, Matthew Arsenault, Travis Desell,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#ifndef _CHISQ_H_
#define _CHISQ_H_

#include "nbody_types.h"
#include "nbody.h"

typedef enum
{
    NBODY_LIKELIHOOD,
    NBODY_ALT_LIKELIHOOD,
    NBODY_CHISQ_ALT_LIKELIHOOD,
    NBODY_POISSON_LIKELIHOOD,
    NBODY_KOLMOGOROV_DISTANCE,
    NBODY_KULLBACK_LEIBLER_DISTANCE,
    NBODY_W_LIKELIHOOD,
    NBODY_W_ALT_LIKELIHOOD
} NBodyLikelihoodMethod;



#ifdef __cplusplus
extern "C" {
#endif

double nbCalcChisq(const NBodyHistogram* data,        /* Data histogram */
                   const NBodyHistogram* histogram,   /* Generated histogram */
                   NBodyLikelihoodMethod method);
double nbMatchEMD(const NBodyHistogram* data, const NBodyHistogram* histogram);

double nbSystemChisq(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf);

NBodyHistogram* nbReadHistogram(const char* histogramFile);

#ifdef __cplusplus
}
#endif

#endif /* _CHISQ_H_ */

