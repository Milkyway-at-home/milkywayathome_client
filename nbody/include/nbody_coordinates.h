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

#ifndef _NBODY_COORDINATES_H_
#define _NBODY_COORDINATES_H_

#include "milkyway_math.h"
#include "nbody_types.h"

typedef struct
{
    real cosphi;
    real sinphi;
    real sinpsi;
    real cospsi;
    real costh;
    real sinth;
} NBHistTrig;

#ifdef __cplusplus
extern "C" {
#endif

/* Coordinate conversion */
mwvector cartesianToLbr(mwvector r, real sunGCDist);
mwvector cartesianToLbr_rad(mwvector r, real sunGCDist);
mwvector lbrToCartesian(mwvector lbr, real sunGCDist);
mwvector lbrToCartesian_rad(mwvector lbr, real sunGCDist);

void nbGetHistTrig(NBHistTrig* ht, const HistogramParams* hp);
real nbXYZToLambda(const NBHistTrig* ht, mwvector xyz, real runGCDist);

#ifdef __cplusplus
}
#endif


#endif /* _NBODY_COORDINATES_H_ */

