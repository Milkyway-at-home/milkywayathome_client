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
#include "coordinates.h"
#include "integrals_common.h"

#pragma OPENCL EXTENSION cl_amd_printf : enable

#define FAST_HPROB 1
#define AUX_BG_PROFILE 0
#define ZERO_BG_PROB 0

#if FAST_HPROB
  #define h_prob_f h_prob_fast
#else
  #define h_prob_f h_prob_slow
#endif /* FAST_H_PROB */


__attribute__ ((always_inline))
inline void stream_sums_cl(real* st_probs,
                           __constant STREAM_CONSTANTS* sc,
                           const mwvector xyz,
                           const real qw_r3_N,
                           const unsigned int nstreams)

{
    //unsigned int i;
    //for (i = 0; i < nstreams; ++i)
    st_probs[0] += calc_st_prob_inc(&sc[0], xyz, qw_r3_N);
    st_probs[1] += calc_st_prob_inc(&sc[1], xyz, qw_r3_N);
    st_probs[2] += calc_st_prob_inc(&sc[2], xyz, qw_r3_N);

}

__attribute__ ((always_inline))
inline void mult_probs_cl(real* st_probs, const real V_reff_xr_rp3, const unsigned int n_stream)
{
    //unsigned int i;
    //for (i = 0; i < n_stream; ++i)

    st_probs[0] *= V_reff_xr_rp3;
    st_probs[1] *= V_reff_xr_rp3;
    st_probs[2] *= V_reff_xr_rp3;
}



__attribute__ ((always_inline))
inline real bg_probability(__constant ASTRONOMY_PARAMETERS* ap,
                           __constant STREAM_CONSTANTS* sc,
                           __constant real* sg_dx,
                           __global const R_POINTS* r_pts,
                           const LB_TRIG lbt,
                           const real gPrime,
                           const unsigned int convolve,
                           real* st_probs)
{
    unsigned int i;
    real rg;
    mwvector xyz = ZERO_VECTOR;
    real bg_prob = 0.0;
    R_POINTS r_pt;


  #if ZERO_BG_PROB
    /* if q is 0, there is no probability */
    /* CHECKME: What happens to the st_probs? */
    return -1.0;
  #endif /* ZERO_BG_PROB */


    for (i = 0; i < convolve; ++i)
    {
        r_pt = r_pts[i];
        //r_pt.r_point = gPrime * bg_prob; // attempt to estimate cost of load
        //r_pt.qw_r3_N = lbt.lcos * gPrime;

        xyz = lbr2xyz_2(r_pt.r_point, lbt);

        rg = rg_calc(xyz, ap->q_inv_sqr);

        bg_prob += h_prob_f(ap, r_pt.qw_r3_N, rg);

      #if AUX_BG_PROFILE
        /* Add a quadratic term in g to the Hernquist profile */
        real g = gPrime + sg_dx[i];
        bg_prob += aux_prob(ap, r_pt.qw_r3_N, g);
      #endif /* AUX_BG_PROFILE */

        //stream_sums(st_probs, sc, xyz, r_pt.qw_r3_N, ap->number_streams);
        stream_sums_cl(st_probs, sc, xyz, r_pt.qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

__attribute__ ((always_inline))
real r_calculation(__constant ASTRONOMY_PARAMETERS* ap,
                   __constant STREAM_CONSTANTS* sc,
                   __constant real* sg_dx,
                   __constant R_CONSTS* rcs,
                   __global const R_POINTS* r_pts,
                   const LB_TRIG lbt,
                   const real id,
                   real* st_probs,
                   const unsigned int r_step)
{
    real bg_prob = bg_probability(ap, sc, sg_dx,
                                  &r_pts[r_step * ap->convolve],
                                  lbt,
                                  rcs[r_step].gPrime,
                                  ap->convolve,
                                  st_probs);

    real V = id * rcs[r_step].irv;
    real V_reff_xr_rp3 = V * rcs[r_step].reff_xr_rp3;

    mult_probs_cl(st_probs, V_reff_xr_rp3, ap->number_streams);

    return V_reff_xr_rp3 * bg_prob;
}

__attribute__ ((always_inline))
inline void write_st_probs(__global real* probs_out,
                           const real* st_probs,
                           const unsigned int n_streams)
{
    //unsigned int i;
    //for (i = 0; i < n_streams; ++i)

    probs_out[0] = st_probs[0];
    probs_out[1] = st_probs[1];
    probs_out[2] = st_probs[2];
}

__kernel void mu_sum_kernel(__global real* restrict mu_out,
                            __global real* restrict probs_out,

                            __constant ASTRONOMY_PARAMETERS* ap,
                            __constant INTEGRAL_AREA* ia,
                            __constant STREAM_CONSTANTS* sc,
                            __constant R_CONSTS* rcs,
                            __global const LB_TRIG* lbts,
                            __global const R_POINTS* r_pts,
                            __constant real* sg_dx,
                            const unsigned int nu_step)

{
    NU_ID nuid;
    LB_TRIG lbt;
    real r_result;

    //real st_probs[3];     /* FIXME: hardcoded stream limit */
    real st_probs[3] = { 0.0, 0.0, 0.0 };

    size_t idx;     /* index into stream probs output buffer */
    size_t mu_step = get_global_id(0);
    size_t r_step = get_global_id(1);

    //zero_st_probs(st_probs, ap->number_streams);

    /* Actual calculations */
    nuid = calc_nu_step(ia, nu_step);

    lbt = lbts[nu_step * ia->mu_steps + mu_step];

    r_result = r_calculation(ap, sc, sg_dx, rcs, r_pts, lbt, nuid.id, st_probs, r_step);

    /* Output to buffers */
    mu_out[mu_step * ia->r_steps + r_step] = r_result;
    idx = mu_step * ia->r_steps * ap->number_streams + r_step * ap->number_streams;
    write_st_probs(&probs_out[idx],
                   st_probs,
                   ap->number_streams);

}


