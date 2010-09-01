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
#include "r_points.h"

#include "milkyway_cl.h"
#include "milkyway_math.h"
#include "milkyway_extra.h"

#include "integrals_likelihood.h"
#include "coordinates.h"
#include "integrals_common.h"

__attribute__ ((always_inline))
inline void write_st_probs(__global ST_PROBS* probs_out,
                           const ST_PROBS* probs,
                           const unsigned int number_streams)
{
    unsigned int i;
    for (i = 0; i < number_streams; ++i)
        probs_out[i] = probs[i];
}

__attribute__ ((always_inline))
inline _MW_STATIC void nu_sum(__global BG_PROB* mu_out,
                              __global ST_PROBS* probs_out,
                              __constant ASTRONOMY_PARAMETERS* ap,
                              __constant STREAM_CONSTANTS* sc,
                              __constant INTEGRAL_AREA* ia,
                              const real irv,
                              const real reff_xr_rp3,
                              __constant R_POINTS* r_pts,
                              real* st_probs,
                              ST_PROBS* probs,
                              const size_t r_step)
{
    unsigned int i;
    BG_PROB mu_result;
    size_t nu_step = get_global_id(1);

    real nu = ia->nu_min + (nu_step * ia->nu_step_size);

    real tmp1 = d2r(90.0 - nu - ia->nu_step_size);
    real tmp2 = d2r(90.0 - nu);

    real id = mw_cos(tmp1) - mw_cos(tmp2);
    nu += 0.5 * ia->nu_step_size;

    mu_result = mu_sum(ap,
                       sc,
                       r_pts,
                       irv,
                       reff_xr_rp3,
                       id,
                       nu,
                       ia->mu_steps,
                       ia->mu_step_size,
                       ia->mu_min,
                       st_probs,
                       probs);

    /* Write results out to buffers */
    mu_out[r_step * ia->nu_steps + nu_step] = mu_result;
    write_st_probs(&probs_out[nu_step * ap->number_streams], probs, ap->number_streams);
}

__kernel void r_sum_kernel(__global BG_PROB* mu_out,
                           __global ST_PROBS* probs_out,

                           __constant ASTRONOMY_PARAMETERS* ap,
                           __constant INTEGRAL_AREA* ia,
                           __constant STREAM_CONSTANTS* sc,
                           __constant R_POINTS* r_pts_all,
                           __local R_POINTS* r_pts)
{
    unsigned int i;
    __constant R_POINTS* r_pts_this;
    __global ST_PROBS* nu_probs;
    size_t r_step = get_global_id(0);

    real st_probs[3]; /* FIXME: hardcoded stream limit */
    ST_PROBS probs[3];

    if (r_step > ia->r_steps)
        return;

    for (i = 0; i < ap->number_streams; ++i)
    {
        probs[i].st_prob_int = 0.0;
        probs[i].st_prob_int_c = 0.0;
    }

    #if 0
    /* Load r_pts into local memory */
    r_pts_this = &r_pts_all[r_step * ap->convolve];
    for (i = 0; i < ap->convolve; ++i)
        r_pts[i] = r_pts_this[i];

    mem_fence(CLK_LOCAL_MEM_FENCE);
    #endif

    r_pts_this = &r_pts_all[r_step * ap->convolve];
    nu_probs = &probs_out[r_step * ia->nu_steps * ap->number_streams];

    R_PRIME rp;
    real reff_xr_rp3;

    rp = calcRPrime(ia, r_step);
    reff_xr_rp3 = calcReffXrRp3(rp.rPrime);

    //nu_sum(mu_out, nu_probs, ap, sc, ia, rp.irv, reff_xr_rp3, r_pts, st_probs, probs, r_step);
    nu_sum(mu_out, nu_probs, ap, sc, ia, rp.irv, reff_xr_rp3, r_pts_this, st_probs, probs, r_step);
}

