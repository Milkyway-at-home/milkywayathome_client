/*
Copyright 2008-2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
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

#if DOUBLEPREC
  #ifdef __ATI_CL__
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
  #else
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
  #endif /* __ATI_CL__ */
#endif /* DOUBLEPREC */

#ifdef __FAST_RELAXED_MATH__
  #error "Bad bad bad bad bad"
#endif /* __FAST_RELAXED_MATH__ */


#include "separation_types.h"
//#include "r_points.h"
#include "r_points.c"

#include "milkyway_cl.h"
#include "milkyway_math.h"
#include "milkyway_extra.h"

#include "integrals_likelihood.h"
#include "coordinates.h"
#include "integrals_common.h"

__attribute__ ((always_inline))
inline _MW_STATIC BG_PROB nu_sum(__constant ASTRONOMY_PARAMETERS* ap,
                                 __constant STREAM_CONSTANTS* sc,
                                 __constant INTEGRAL_AREA* ia,
                                 const real irv,
                                 const real reff_xr_rp3,
                                 __local const R_POINTS* r_pts,
                                 __constant NU_CONSTANTS* nu_consts,
                                 __local real* st_probs,
                                 __local ST_PROBS* probs)
{
    unsigned int nu_step;
    BG_PROB mu_result;
    BG_PROB nu_acc = ZERO_BG_PROB;

    const unsigned int nu_steps = ia->nu_steps;
    const unsigned int mu_steps = ia->mu_steps;
    const real mu_min = ia->mu_min;
    const real mu_step_size = ia->mu_step_size;

    for (nu_step = 0; nu_step < nu_steps; ++nu_step)
    {
        mu_result = mu_sum(ap,
                           sc,
                           r_pts,
                           irv,
                           reff_xr_rp3,
                           nu_consts[nu_step].id,
                           nu_consts[nu_step].nu,
                           mu_steps,
                           mu_step_size,
                           mu_min,
                           st_probs,
                           probs);

        INCADD_BG_PROB(nu_acc, mu_result);
    }

    return nu_acc;
}

__attribute__ ((always_inline))
inline _MW_STATIC BG_PROB r_sum(__constant ASTRONOMY_PARAMETERS* ap,
                                __constant INTEGRAL_AREA* ia,
                                __constant STREAM_CONSTANTS* sc,
                                __constant STREAM_GAUSS* sg,
                                __constant NU_CONSTANTS* nu_consts,
                                __local R_POINTS* r_pts,
                                __local real* st_probs,
                                __local ST_PROBS* probs,
                                const unsigned int r_step)
{
    R_PRIME rp;
    real reff_xr_rp3;

    rp = calcRPrime(ia, r_step);
    reff_xr_rp3 = calcReffXrRp3(rp.rPrime);

    return nu_sum(ap, sc, ia, rp.irv, reff_xr_rp3, r_pts, nu_consts, st_probs, probs);
}

__kernel void r_sum_kernel(__global BG_PROB* nu_out,
                           __global ST_PROBS* probs_out,

                           __constant ASTRONOMY_PARAMETERS* ap,
                           __constant INTEGRAL_AREA* ia,
                           __constant STREAM_CONSTANTS* sc,
                           __constant STREAM_GAUSS* sg,
                           __constant NU_CONSTANTS* nu_consts,
                           __constant R_POINTS* r_pts_all,

                           __local real* st_probs,
                           __local ST_PROBS* probs,
                           __local R_POINTS* r_pts)
{
    unsigned int i;
    BG_PROB nu_result;
    __constant R_POINTS* r_pts_this;
    size_t r_step = get_global_id(0);

    if (r_step > ia->r_steps)
        return;

    for (i = 0; i < ap->number_streams; ++i)
    {
        probs[i].st_prob_int = 0.0;
        probs[i].st_prob_int_c = 0.0;
    }

    /* Load r_pts into local memory */
    r_pts_this = &r_pts_all[r_step * ap->convolve];
    for (i = 0; i < ap->convolve; ++i)
        r_pts[i] = r_pts_this[i];

    nu_result = r_sum(ap, ia, sc, sg, nu_consts, r_pts, st_probs, probs, r_step);

    /* Write results back */
    nu_out[r_step] = nu_result;
    for (i = 0; i < ap->number_streams; ++i)
        probs_out[r_step * ap->number_streams + i] = probs[i];
}

