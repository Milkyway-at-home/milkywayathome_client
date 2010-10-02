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
#define RL_1_2 ((real) 0.5)
#define RL3 ((real) 3.0)
#define RL_1_3 ((real) 0.33333333333333333333333333)
#define RL5 ((real) 5.0)
#define RL_1_5 ((real) 0.2)
#define RL1000 ((real) 1000.0)
#define RL_1_1000 ((real) 0.001)

ALWAYS_INLINE CONST_F
inline real distance_magnitude(const real m)
{
    return mw_exp10((m - (real) 14.2) * RL_1_5);
}

inline R_PRIME calcRPrime(__MW_CONSTANT INTEGRAL_AREA* ia, const unsigned int r_step)
{
    real r, next_r, log_r;
    R_PRIME ret;

    log_r = ia->r_min + (r_step * ia->r_step_size);
    r = distance_magnitude(log_r);
    next_r = distance_magnitude(log_r + ia->r_step_size);

    ret.irv = d2r(((cube(next_r) - cube(r)) * RL_1_3) * ia->mu_step_size);
    ret.rPrime = RL_1_2 * (next_r + r);

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

ALWAYS_INLINE
inline R_CONSTS calcRConsts(R_PRIME rp)
{
    R_CONSTS rc;

    rc.reff_xr_rp3 = calcReffXrRp3(rp.rPrime);
    rc.irv = rp.irv;

    return rc;
}

ALWAYS_INLINE HOT
inline R_POINTS calc_r_point(__MW_CONSTANT STREAM_GAUSS* sg, const real gPrime, const real coeff)
{
    R_POINTS r_pt;
    real g, exponent, r3, N;

    g = gPrime + sg->dx;

    /* MAG2R */
    r_pt.r_in_mag = g;
    r_pt.r_in_mag2 = sqr(g);
    r_pt.r_point = RL_1_1000 * mw_exp10(RL_1_5 * (g - absm) + RL1);

    r3 = cube(r_pt.r_point);
    exponent = sqr(g - gPrime) * inv(RL2 * sqr(stdev));
    N = coeff * mw_exp(-exponent);
    r_pt.qw_r3_N = sg->qgaus_W * r3 * N;

    return r_pt;
}

void set_r_points(__MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                  __MW_CONSTANT STREAM_GAUSS* sg,
                  const unsigned int n_convolve,
                  const real coords,
                  __MW_LOCAL R_POINTS* r_pts);

R_POINTS* precalculate_r_pts(const ASTRONOMY_PARAMETERS* ap,
                             const INTEGRAL_AREA* ia,
                             const STREAM_GAUSS* sg,
                             R_CONSTS** rc_out);

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRAL_CONSTANTS_H_ */

