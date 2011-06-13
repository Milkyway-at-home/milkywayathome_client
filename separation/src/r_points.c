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
#include "milkyway_extra.h"
#include "r_points.h"

#include "milkyway_util.h"

static inline real distance_magnitude(const real m)
{
    return mw_exp10((m - (real) 14.2) * 0.2);
}

static RPrime calcRPrime(const IntegralArea* ia, unsigned int r_step)
{
    real r, next_r, log_r;
    RPrime ret;

    log_r = ia->r_min + (r_step * ia->r_step_size);
    r = distance_magnitude(log_r);
    next_r = distance_magnitude(log_r + ia->r_step_size);

    ret.irv = d2r(((cube(next_r) - cube(r)) * (1.0 / 3.0)) * ia->mu_step_size);
    ret.rPrime = 0.5 * (next_r + r);

    return ret;
}

/* This applies a sigmoid to account for the falloff in the SDSS data. Heidi
   Described it in 2002. Check Nates MW Bible */
real calcReffXrRp3(real coords, real gPrime)
{
    static const real sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

    /* REFF */
    const real exp_result = mw_exp(sigmoid_curve_params[1] * (gPrime - sigmoid_curve_params[2]));
    const real reff_value = sigmoid_curve_params[0] / (exp_result + 1.0);
    const real rPrime3 = cube(coords);
    const real reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

real calcG(real coords)
{
    return 5.0 * (mw_log10(1000.0 * coords) - 1.0) + absm;
}

static inline RPoints calc_r_point(real dx, real qgaus_W, real gPrime, real coeff)
{
    RPoints r_pt;

    real g, exponent, r3, N, stddev_l, stddev_r, stddev_i, A;

    g = gPrime + dx;

    /* MAG2R */
    r_pt.r_point = 0.001 * mw_exp10(0.2 * (g - absm) + 1.0);
    r3 = cube(r_pt.r_point);

    exponent = sqr(g - gPrime) * inv(2.0 * sqr(stdev));
    N = coeff * mw_exp(-exponent);
    r_pt.qw_r3_N = qgaus_W * r3 * N;

#if 0
    /* Reimplemented to account for matt newbys f_turnoff distribution insights */
    stddev_l = 0.315;

    /* Function from Matt
       \alpha = .52, \beta=12.0 \gamma=0.76
       Get d_eff, I assumed it was r_pt.r_point given the simularities
    */
    stddev_r = 0.52 * inv(1.0 + mw_exp(12.0 - r_pt.r_point)) + 0.76;

    /* if g <= \mu = 4.2 we use a constant
       however is gPrime equal to \mu? It seems to be used in the same way
    */

    stddev_i = (g <= absm) ? stddev_l : stddev_r;

    /* Note see previous uncertainty about gPrime versus \mu */
    exponent = sqr(g - gPrime) * inv(2.0 * sqr(stddev_i));

    A = inv(2.0 * M_PI * (stddev_l + stddev_r) * inv(2.0));
    N = A * mw_exp(-exponent);
#endif /* 0 */

    return r_pt;
}

static RConsts calcRConsts(RPrime rp)
{
    RConsts rc;

    rc.gPrime = calcG(rp.rPrime);
    rc.irv_reff_xr_rp3 = rp.irv * calcReffXrRp3(rp.rPrime, rc.gPrime);

    return rc;
}


void setRPoints(const AstronomyParameters* ap,
                const StreamGauss sg,
                unsigned int n_convolve,
                real gPrime,
                RPoints* r_pts)
{
    unsigned int i;

    for (i = 0; i < n_convolve; ++i)
        r_pts[i] = calc_r_point(sg.dx[i], sg.qgaus_W[i], gPrime, ap->coeff);
}

void setSplitRPoints(const AstronomyParameters* ap,
                     const StreamGauss sg,
                     unsigned int n_convolve,
                     real gPrime,
                     real* RESTRICT r_points,
                     real* RESTRICT qw_r3_N)
{
    unsigned int i;
    RPoints rPt;

    for (i = 0; i < n_convolve; ++i)
    {
        rPt = calc_r_point(sg.dx[i], sg.qgaus_W[i], gPrime, ap->coeff);
        r_points[i] = rPt.r_point;
        qw_r3_N[i] = rPt.qw_r3_N;
    }
}

RPoints* precalculateRPts(const AstronomyParameters* ap,
                          const IntegralArea* ia,
                          const StreamGauss sg,
                          RConsts** rc_out,
                          int transpose)
{
    unsigned int i, j, idx;
    RPoints* r_pts;
    RPrime rp;
    RConsts* rc;
    RPoints r_pt;

    size_t rPtsSize = sizeof(RPoints) * ap->convolve * ia->r_steps;
    size_t rConstsSize = sizeof(RConsts) * ia->r_steps;

    r_pts = (RPoints*) mwMallocA(rPtsSize);
    rc = (RConsts*) mwMallocA(rConstsSize);

    for (i = 0; i < ia->r_steps; ++i)
    {
        rp = calcRPrime(ia, i);
        rc[i] = calcRConsts(rp);

        for (j = 0; j < ap->convolve; ++j)
        {
            r_pt = calc_r_point(sg.dx[j], sg.qgaus_W[j], rc[i].gPrime, ap->coeff);
            idx = transpose ? j * ia->r_steps + i : i * ap->convolve + j;
            r_pts[idx] = r_pt;
        }
    }

    *rc_out = rc;
    return r_pts;
}

