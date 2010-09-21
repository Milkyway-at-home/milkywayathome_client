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
#include "milkyway_cl.h"
#include "milkyway_extra.h"


/* Literals are assumed to be doubles by default, and the
 * -cl-single-precision-constant flag seems to not be working when
 * trying to use float */
#define RL1 ((real) 1.0)
#define RL2 ((real) 2.0)
#define RL3 ((real) 3.0)
#define RL5 ((real) 5.0)
#define RL10 ((real) 10.0)
#define RL1000 ((real) 1000.0)

ALWAYS_INLINE CONST_F
inline real distance_magnitude(const real m)
{
    return mw_powr(RL10, (m - (real) 14.2) / RL5);
}

inline R_PRIME calcRPrime(__MW_CONSTANT INTEGRAL_AREA* ia, const unsigned int r_step)
{
    real r, next_r, log_r;
    R_PRIME ret;

    log_r = ia->r_min + (r_step * ia->r_step_size);
    r = distance_magnitude(log_r);
    next_r = distance_magnitude(log_r + ia->r_step_size);

    ret.irv = d2r(((cube(next_r) - cube(r)) / RL3) * ia->mu_step_size);
    ret.rPrime = (next_r + r) / RL2;

    return ret;
}

inline real calcGPrime(const real coords)
{
    return RL5 * (mw_log10(coords * RL1000) - RL1) + absm;
}

inline real calcReffXrRp3(const real coords)
{
    _MW_STATIC const real sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };
    const real gPrime = calcGPrime(coords);

    /* REFF */
    const real exp_result = mw_exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    const real reff_value = sigmoid_curve_params[0] / (exp_result + RL1);
    const real rPrime3 = cube(coords);
    const real reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

void set_r_points(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                  __MW_CONSTANT STREAM_GAUSS* sg,
                  const unsigned int n_convolve,
                  const real coords,
                  __MW_LOCAL R_POINTS* r_pts);

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRAL_CONSTANTS_H_ */

