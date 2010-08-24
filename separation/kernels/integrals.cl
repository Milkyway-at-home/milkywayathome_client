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

#ifdef __FAST_RELAXED_MATH__
  #error "Bad bad bad"
#endif /* __FAST_RELAXED_MATH__ */

#if DOUBLEPREC
  #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif /* DOUBLEPREC */

#include "separation_types.h"
//#include "r_points.h"

#include "milkyway_cl.h"
#include "milkyway_math.h"
#include "milkyway_extra.h"

#include "integrals_likelihood.h"
#include "coordinates.h"
#include "integrals_common.h"


__attribute__ ((always_inline))
inline _MW_STATIC BG_PROB nu_sum(__private const ASTRONOMY_PARAMETERS* ap,
                                 __constant const STREAM_CONSTANTS* sc,
                                 __private const INTEGRAL_AREA* ia,
                                 const real irv,
                                 const real reff_xr_rp3,
                                 __local const R_POINTS* r_pts,
                                 __constant const NU_CONSTANTS* nu_consts,
                                 __local ST_PROBS* probs,
                                 __local vector* xyz)
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
                           probs,
                           xyz);
        INCADD_BG_PROB(nu_acc, mu_result);
    }

    return nu_acc;
}

__attribute__ ((always_inline))
inline _MW_STATIC BG_PROB r_sum(__private const ASTRONOMY_PARAMETERS* ap,
                                __private const INTEGRAL_AREA* ia,
                                __constant const STREAM_CONSTANTS* sc,
                                __constant const STREAM_GAUSS* sg,
                                __constant const NU_CONSTANTS* nu_consts,
                                __local R_POINTS* r_pts,
                                __local ST_PROBS* probs,
                                __local vector* xyz,
                                const unsigned int r_step)
{
    BG_PROB nu_result;
    real r, next_r, rPrime;
    real irv, reff_xr_rp3;

  #ifdef USE_KPC
    const real r_max           = ia->r_min + ia->r_step_size * r_steps;
    const real r_min_kpc       = distance_magnitude(ia->r_min);
    const real r_max_kpc       = distance_magnitude(ia->r_max);
    const real r_step_size_kpc = (r_max_kpc - r_min_kpc) / r_steps;
    r = r_min_kpc + (r_step * r_step_size_kpc);
    next_r = r + r_step_size_kpc;
  #else
    real log_r = ia->r_min + (r_step * ia->r_step_size);
    r = distance_magnitude(log_r);
    next_r = distance_magnitude(log_r + ia->r_step_size);
  #endif /* USE_KPC */

    irv = d2r(((cube(next_r) - cube(r)) / (real) 3.0) * ia->mu_step_size);
    rPrime = (next_r + r) / (real) 2.0;

    reff_xr_rp3 = set_r_points(ap, sg, ap->convolve, rPrime, r_pts);

    nu_result = nu_sum(ap, sc, ia, irv, reff_xr_rp3, r_pts, nu_consts, probs, xyz);

    return nu_result;
}

__kernel void r_sum_kernel(__global BG_PROB* nu_out,
                           __global ST_PROBS* probs_out,

                           __private const ASTRONOMY_PARAMETERS ap,
                           __private const INTEGRAL_AREA ia,
                           __constant STREAM_CONSTANTS* sc,
                           __constant STREAM_GAUSS* sg,
                           __constant NU_CONSTANTS* nu_consts,

                           __local ST_PROBS* probs,
                           __local vector* xyz,
                           __local R_POINTS* r_pts)
{
    unsigned int i;
    BG_PROB nu_result;
    size_t r_step = get_global_id(0);

    if (r_step > ia.r_steps)
        return;

    for (i = 0; i < ap.number_streams; ++i)
    {
        probs[i].st_prob_int = 0.0;
        probs[i].st_prob_int_c = 0.0;
    }

    nu_result = r_sum(&ap, &ia, sc, sg, nu_consts, r_pts, probs, xyz, r_step);

    /* Write results back */
    nu_out[r_step] = nu_result;
    for (i = 0; i < ap.number_streams; ++i)
        probs_out[r_step * ap.number_streams + i] = probs[i];
}

