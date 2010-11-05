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

#ifndef _R_POINTS_H_
#define _R_POINTS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_types.h"
#include "separation_constants.h"
#include "coordinates.h"
#include "milkyway_cl.h"
#include "milkyway_extra.h"

real calcG(const real coords);
real calcReffXrRp3(const real coords, const real gPrime);

void set_r_points(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                  const STREAM_GAUSS sg,
                  const unsigned int n_convolve,
                  const real coords,
                  __MW_LOCAL R_POINTS* r_pts);

R_POINTS* precalculate_r_pts(const ASTRONOMY_PARAMETERS* ap,
                             const INTEGRAL_AREA* ia,
                             const STREAM_GAUSS sg,
                             R_CONSTS** rc_out);

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRAL_CONSTANTS_H_ */

