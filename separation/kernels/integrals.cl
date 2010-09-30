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
  #error "Bad bad bad bad bad"
#endif /* __FAST_RELAXED_MATH__ */

#include "milkyway_math.h"
#include "separation_types.h"
#include "r_points.h"
#include "milkyway_cl.h"
#include "milkyway_extra.h"

#include "integrals_likelihood.h"
#include "coordinates.h"
#include "integrals_common.h"

#pragma OPENCL EXTENSION cl_amd_printf : enable

__attribute__ ((always_inline))
inline void write_st_probs(__global ST_PROBS* probs_out,
                           const ST_PROBS* probs,
                           const unsigned int number_streams)
{
    unsigned int i;

    for (i = 0; i < number_streams; ++i)
    {
        probs_out[i] = probs[i];

       #if 0
        ST_PROBS arst = { 0xabcd, 0xbeef };
        if (isnan(probs[i].st_prob_int) || isnan(probs[i].st_prob_int_c))
            probs_out[i] = arst;
        else
            probs_out[i] = probs[i];
        #endif
    }

}

__attribute__ ((always_inline))
inline _MW_STATIC void nu_sum(__global BG_PROB* mu_out,
                              __global ST_PROBS* probs_out,
                              __MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                              __MW_CONSTANT STREAM_CONSTANTS* sc,
                              __MW_CONSTANT INTEGRAL_AREA* ia,
                              const real irv,
                              const real reff_xr_rp3,
                              __MW_CONSTANT R_POINTS* r_pts,
                              real* st_probs,
                              ST_PROBS* probs,
                              const size_t r_step)
{
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

__attribute__ ((always_inline))
inline _MW_STATIC void zero_stream_probs(ST_PROBS* probs, const unsigned int nstreams)
{
    unsigned int i;

    for (i = 0; i < nstreams; ++i)
    {
        probs[i].st_prob_int = 0.0;
        probs[i].st_prob_int_c = 0.0;
    }
}


__kernel void r_sum_kernel(__global BG_PROB* mu_out,
                           __global ST_PROBS* probs_out,

                           __MW_CONSTANT ASTRONOMY_PARAMETERS* ap,
                           __MW_CONSTANT INTEGRAL_AREA* ia,
                           __MW_CONSTANT STREAM_CONSTANTS* sc,
                           __MW_CONSTANT STREAM_GAUSS* sg,

{
    __MW_CONSTANT R_POINTS* r_pts_this;
    __global ST_PROBS* nu_probs;
    size_t r_step = get_global_id(0);

    real st_probs[3];   /* FIXME: hardcoded stream limit */
    ST_PROBS probs[3];

    if (r_step > ia->r_steps)
        return;

    zero_stream_probs(probs, ap->number_streams);

    #if 0
    /* Load r_pts into local memory */
    r_pts_this = &r_pts_all[r_step * ap->convolve];
    for (i = 0; i < ap->convolve; ++i)
        r_pts[i] = r_pts_this[i];

    mem_fence(CLK_LOCAL_MEM_FENCE);
    #endif

    R_PRIME rp;
    real reff_xr_rp3;

    rp = calcRPrime(ia, r_step);
    reff_xr_rp3 = calcReffXrRp3(rp.rPrime);

    nu_probs = &probs_out[r_step * ia->nu_steps * ap->number_streams];

    nu_sum(mu_out, nu_probs, ap, sc, ia, rp.irv, reff_xr_rp3, sg, st_probs, probs, r_step);

}

