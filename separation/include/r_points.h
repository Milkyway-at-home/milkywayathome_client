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

ALWAYS_INLINE CONST_F OLD_GCC_EXTERNINLINE
inline real distance_magnitude(const real m)
{
    return mw_exp10((m - (real) 14.2) * 0.2);
}

OLD_GCC_EXTERNINLINE
inline R_PRIME calcRPrime(__MW_CONSTANT INTEGRAL_AREA* ia, const unsigned int r_step)
{
    real r, next_r, log_r;
    R_PRIME ret;

    log_r = ia->r_min + (r_step * ia->r_step_size);
    r = distance_magnitude(log_r);
    next_r = distance_magnitude(log_r + ia->r_step_size);

    ret.irv = d2r(((cube(next_r) - cube(r)) * (1.0 / 3.0)) * ia->mu_step_size);
    ret.rPrime = 0.5 * (next_r + r);

    return ret;
}

OLD_GCC_EXTERNINLINE
inline real calcReffXrRp3(const real coords, const real gPrime)
{
    static const real sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

    /* REFF */
    const real exp_result = mw_exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    const real reff_value = sigmoid_curve_params[0] / (exp_result + 1.0);
    const real rPrime3 = cube(coords);
    const real reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

ALWAYS_INLINE OLD_GCC_EXTERNINLINE
inline R_CONSTS calcRConsts(R_PRIME rp)
{
    R_CONSTS rc;
    rc.gPrime = calcG(rp.rPrime);
    rc.reff_xr_rp3 = calcReffXrRp3(rp.rPrime, rc.gPrime);
    rc.irv = rp.irv;

    return rc;
}

ALWAYS_INLINE HOT OLD_GCC_EXTERNINLINE
inline R_POINTS calc_r_point(const real dx, const real qgaus_W, const real gPrime, const real coeff)
{
    R_POINTS r_pt;
    real g, exponent, r3, N;

    g = gPrime + dx;

    /* MAG2R */
    r_pt.r_point = 0.001 * mw_exp10(0.2 * (g - absm) + 1.0);

    r3 = cube(r_pt.r_point);
    exponent = sqr(g - gPrime) * inv(2.0 * sqr(stdev));
    N = coeff * mw_exp(-exponent);
    r_pt.qw_r3_N = qgaus_W * r3 * N;

    return r_pt;
}

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

