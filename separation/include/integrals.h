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

#ifndef _INTEGRALS_H_
#define _INTEGRALS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_types.h"
#include "evaluation_state.h"

real bg_probability_fast_hprob(const ASTRONOMY_PARAMETERS* ap,
                               const STREAM_CONSTANTS* sc,
                               const R_POINTS* r_pts,
                               const real* sg_dx,
                               const LB_TRIG lbt,
                               const real gPrime,
                               const int aux_bg_profile,
                               const unsigned int convolve,

                               real* st_probs);

real bg_probability_slow_hprob(const ASTRONOMY_PARAMETERS* ap,
                               const STREAM_CONSTANTS* sc,
                               const R_POINTS* r_pts,
                               const real* sg_dx,
                               const LB_TRIG lbt,
                               const real gPrime,
                               const int aux_bg_profile,
                               const unsigned int convolve,
                               real* st_probs);

NU_ID calc_nu_step(const INTEGRAL_AREA* ia, const unsigned int nu_step);

real integrate(const ASTRONOMY_PARAMETERS* ap,
               const INTEGRAL_AREA* ia,
               const STREAM_CONSTANTS* sc,
               const STREAM_GAUSS sg,
               KAHAN* probs,
               EVALUATION_STATE* es);

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRALS_H_ */

