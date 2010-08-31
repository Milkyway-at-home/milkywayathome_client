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
inline _MW_STATIC BG_PROB nu_sum(__global BG_PROB* mu_out,
                                 __global ST_PROBS* probs_out,
                                 __constant ASTRONOMY_PARAMETERS* ap,
                                 __constant STREAM_CONSTANTS* sc,
                                 __constant INTEGRAL_AREA* ia,
                                 const real irv,
                                 const real reff_xr_rp3,
                                 __local const R_POINTS* r_pts,
                                 __constant NU_CONSTANTS* nu_consts,
                                 __local real* st_probs,
                                 __local ST_PROBS* probs,
                                 const size_t r_step)
{
    unsigned int i;
    BG_PROB mu_result;
    size_t nu_step = get_global_id(1);

    mu_result = mu_sum(ap,
                       sc,
                       r_pts,
                       irv,
                       reff_xr_rp3,
                       nu_consts[nu_step].id,
                       nu_consts[nu_step].nu,
                       ia->mu_steps,
                       ia->mu_step_size,
                       ia->mu_min,
                       st_probs,
                       probs);

    mu_out[r_step * ia->nu_steps + nu_step] = mu_result;

    /* Write results back */
    for (i = 0; i < ap->number_streams; ++i)
        probs_out[nu_step * ap->number_streams + i] = probs[i];
}

__attribute__ ((always_inline))
inline _MW_STATIC void r_sum(__global BG_PROB* mu_out,
                             __global ST_PROBS* probs_out,
                             __constant ASTRONOMY_PARAMETERS* ap,
                             __constant INTEGRAL_AREA* ia,
                             __constant STREAM_CONSTANTS* sc,
                             __constant NU_CONSTANTS* nu_consts,
                             __local R_POINTS* r_pts,
                             __local real* st_probs,
                             __local ST_PROBS* probs,
                             const size_t r_step)
{
    R_PRIME rp;
    real reff_xr_rp3;

    rp = calcRPrime(ia, r_step);
    reff_xr_rp3 = calcReffXrRp3(rp.rPrime);

    nu_sum(mu_out, probs_out, ap, sc, ia, rp.irv, reff_xr_rp3, r_pts, nu_consts, st_probs, probs, r_step);
}

__kernel void r_sum_kernel(__global BG_PROB* mu_out,
                           __global ST_PROBS* probs_out,

                           __constant ASTRONOMY_PARAMETERS* ap,
                           __constant INTEGRAL_AREA* ia,
                           __constant STREAM_CONSTANTS* sc,
                           __constant NU_CONSTANTS* nu_consts,
                           __constant R_POINTS* r_pts_all,

                           __local real* st_probs,
                           __local ST_PROBS* probs,
                           __local R_POINTS* r_pts)
{
    unsigned int i;
    BG_PROB nu_result;
    __constant R_POINTS* r_pts_this;
    __global ST_PROBS* nu_probs;
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

    barrier(CLK_LOCAL_MEM_FENCE);

    nu_probs = &probs_out[r_step * ia->nu_steps * ap->number_streams];
    r_sum(mu_out, nu_probs, ap, ia, sc, nu_consts, r_pts, st_probs, probs, r_step);
}

