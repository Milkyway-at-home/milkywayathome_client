/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole, Matthew
Arsenault, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#include "separation_types.h"
#include "milkyway_cl.h"
#include "milkyway_extra.h"
#include "r_points.h"

#include "milkyway_util.h"

static inline real distance_magnitude(const real m)
{
    return mw_exp10((m - (real) 14.2) * 0.2);
}

static R_PRIME calcRPrime(const INTEGRAL_AREA* ia, const unsigned int r_step)
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

real calcReffXrRp3(const real coords, const real gPrime)
{
    static const real sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

    /* REFF */
    const real exp_result = mw_exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    const real reff_value = sigmoid_curve_params[0] / (exp_result + 1.0);
    const real rPrime3 = cube(coords);
    const real reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

real calcG(const real coords)
{
    return 5.0 * (mw_log10(1000.0 * coords) - 1.0) + absm;
}

static inline R_POINTS calc_r_point(const real dx, const real qgaus_W, const real gPrime, const real coeff)
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

static inline R_CONSTS calcRConsts(R_PRIME rp)
{
    R_CONSTS rc;

    rc.gPrime = calcG(rp.rPrime);
    rc.irv_reff_xr_rp3 = rp.irv * calcReffXrRp3(rp.rPrime, rc.gPrime);

    return rc;
}

void set_r_points(const ASTRONOMY_PARAMETERS* ap,
                  const STREAM_GAUSS sg,
                  const unsigned int n_convolve,
                  const real gPrime,
                  R_POINTS* r_pts)
{
    unsigned int i;

    for (i = 0; i < n_convolve; ++i)
        r_pts[i] = calc_r_point(sg.dx[i], sg.qgaus_W[i], gPrime, ap->coeff);
}

R_POINTS* precalculate_r_pts(const ASTRONOMY_PARAMETERS* ap,
                             const INTEGRAL_AREA* ia,
                             const STREAM_GAUSS sg,
                             R_CONSTS** rc_out)
{
    unsigned int i;
    R_POINTS* r_pts;
    R_PRIME rp;
    R_CONSTS* rc;

    size_t rPtsSize = sizeof(R_POINTS) * ap->convolve * ia->r_steps;
    size_t rConstsSize = sizeof(R_CONSTS) * ia->r_steps;

    r_pts = (R_POINTS*) mwMallocAligned(rPtsSize, sizeof(R_POINTS));
    rc = (R_CONSTS*) mwMallocAligned(rConstsSize, sizeof(R_CONSTS));

    for (i = 0; i < ia->r_steps; ++i)
    {
        rp = calcRPrime(ia, i);
        rc[i] = calcRConsts(rp);
        set_r_points(ap, sg, ap->convolve, rc[i].gPrime, &r_pts[i * ap->convolve]);
    }

    *rc_out = rc;
    return r_pts;
}

