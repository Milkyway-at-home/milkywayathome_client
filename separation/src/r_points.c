/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole, Dave Przybylo
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

    /*See the discussion of xr where it is declared in separation_constants.h for more details.  We are leaving the
    function and variable names as reff_xr_rp3 so that reverting back to include this value is trivial.  For future
    code clean-up it should be noted that this can instead be changed to reff_rp3 the way it is now to more clearly
    reflect the quantity.*/
    //const real reff_xr_rp3 = reff_value * xr / rPrime3;
    const real reff_xr_rp3 = reff_value / rPrime3;

    return reff_xr_rp3;
}

real calcG(real coords)
{
    return 5.0 * (mw_log10(1000.0 * coords) - 1.0) + absm;
}

static inline RPoints calc_r_point(real dx, real qgaus_W, RConsts* rc)
{
    RPoints r_pt;
    real g, exponent, r3, N;
    real absm_u, stdev_l, stdev_i;

    absm_u = absm;
    stdev_i = stdev;

    g = rc->gPrime + dx;

//#if sf.modfit
#if 1
    /* Implement modified f_turnoff distribution described in Newby 2011*/
    absm_u = 4.18;
    stdev_l = 0.36;
    stdev_i = (g <= rc->gPrime) ? stdev_l : rc->stdev_r;
#endif 

    /* MAG2R */
    r_pt.r_point = 0.001 * mw_exp10(0.2 * (g - absm_u) + 1.0);
    r3 = cube(r_pt.r_point);

    exponent = sqr(g - rc->gPrime) * inv(2.0 * sqr(stdev_i));
    N = rc->coeff * mw_exp(-exponent);
    r_pt.qw_r3_N = qgaus_W * r3 * N;

    return r_pt;
}

//Ignore calculation of irv_reff_xr_rp3.  Use for efficiency when this is not needed.
RConsts calcRConsts(real coords)
{
    RConsts rc;
    real stdev_o;

    rc.gPrime = calcG(coords);
    rc.irv_reff_xr_rp3 = 0.0;
    rc.stdev_r = stdev;
    stdev_o = rc.stdev_r;

    //#if sf.modfit
    #if 1
        /* Implement modified f_turnoff distribution described in Newby 2011*/
        const real stdev_l = 0.36;
        const real alpha = 0.52;
        const real beta = 12.0;
        const real gam = 0.76;

        rc.stdev_r = alpha * inv(1.0 + mw_exp(beta - rc.gPrime)) + gam;
        stdev_o = 0.5 * (stdev_l + rc.stdev_r);
    #endif

    rc.coeff = 1.0 / (stdev_o * SQRT_2PI);
    
    return rc;
}

//Include calculation of irv_reff_xr_rp3
static RConsts calcRConstsplus(RPrime rp)
{
    RConsts rc;
    real stdev_o;

    rc.gPrime = calcG(rp.rPrime);
    rc.irv_reff_xr_rp3 = rp.irv * calcReffXrRp3(rp.rPrime, rc.gPrime);
    rc.stdev_r = stdev;
    stdev_o = rc.stdev_r;

    //#if sf.modfit
    #if 1
        /* Implement modified f_turnoff distribution described in Newby 2011*/
        const real stdev_l = 0.36;
        const real alpha = 0.52;
        const real beta = 12.0;
        const real gam = 0.76;

        rc.stdev_r = alpha * inv(1.0 + mw_exp(beta - rc.gPrime)) + gam;
        stdev_o = 0.5 * (stdev_l + rc.stdev_r);
    #endif

    rc.coeff = 1.0 / (stdev_o * SQRT_2PI);
    
    return rc;
}

void setRPoints(const StreamGauss sg,
                unsigned int n_convolve,
                RConsts* rc,
                RPoints* r_pts)
{
    unsigned int i;

    for (i = 0; i < n_convolve; ++i)
        r_pts[i] = calc_r_point(sg.dx[i], sg.qgaus_W[i], rc);
}

void setSplitRPoints(const StreamGauss sg,
                     unsigned int n_convolve,
                     RConsts* rc,
                     real* RESTRICT r_points,
                     real* RESTRICT qw_r3_N)
{
    unsigned int i;
    RPoints rPt;

    for (i = 0; i < n_convolve; ++i)
    {
        rPt = calc_r_point(sg.dx[i], sg.qgaus_W[i], rc);
        r_points[i] = rPt.r_point;
        qw_r3_N[i] = rPt.qw_r3_N;
    }
}

RPoints* precalculateRPts(unsigned int n_convolve,
                          const IntegralArea* ia,
                          const StreamGauss sg,
                          RConsts** rc_out,
                          int transpose)
{
    unsigned int i, idx;
    int j;
    RPoints* r_pts;
    RPrime rp;
    RConsts* rc;
    RPoints r_pt;

    size_t rPtsSize = sizeof(RPoints) * n_convolve * ia->r_steps;
    size_t rConstsSize = sizeof(RConsts) * ia->r_steps;

    r_pts = (RPoints*) mwMallocA(rPtsSize);
    rc = (RConsts*) mwMallocA(rConstsSize);

    for (i = 0; i < ia->r_steps; ++i)
    {
        rp = calcRPrime(ia, i);
        rc[i] = calcRConstsplus(rp);

        for (j = 0; j < n_convolve; ++j)
        {
            r_pt = calc_r_point(sg.dx[j], sg.qgaus_W[j], &rc[i]);
            idx = transpose ? j * ia->r_steps + i : i * n_convolve + j;
            r_pts[idx] = r_pt;
        }
    }

    *rc_out = rc;
    return r_pts;
}

