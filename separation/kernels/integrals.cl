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

#define FAST_HPROB 1
#define AUX_BG_PROFILE 0

__attribute__ ((always_inline))
inline real sub_bg_probability(__constant ASTRONOMY_PARAMETERS* ap,
                               __constant STREAM_CONSTANTS* sc,
                               __global const R_POINTS* r_pts,
                               const real gPrime,
                               const LB_TRIG lbt,
                               const unsigned int convolve,
                               real* st_probs)
{
    unsigned int i;
    real h_prob;
    real rg, rs;
    vector xyz;
    real bg_prob = 0.0;
    R_POINTS r_pt;

    for (i = 0; i < convolve; ++i)
    {
        r_pt = r_pts[i];
        lbr2xyz_2(xyz, r_pt.r_point, lbt);

        rg = rg_calc(xyz, ap->q_inv_sqr);
        rs = rg + ap->r0;

      #if FAST_HPROB
        h_prob = h_prob_fast(r_pt.qw_r3_N, rg, rs);

        #if AUX_BG_PROFILE
        /* the Hernquist profile includes a quadratic term in g */
        h_prob += aux_prob(ap, r_pt.qw_r3_N, r_pt.r_in_mag, r_pt.r_in_mag2);
        #endif /* AUX_BG_PROFILE */

        bg_prob += h_prob;

      #else
        bg_prob += h_prob_slow(ap, r_pt.qw_r3_N, rg);
      #endif /* FAST_HPROB */

        stream_sums(st_probs, sc, xyz, r_pt.qw_r3_N, ap->number_streams);
    }

    return bg_prob;
}

__attribute__ ((always_inline))
inline void mult_probs(real* st_probs, const real V_reff_xr_rp3, const unsigned int n_stream)
{
    unsigned int i;

    for (i = 0; i < n_stream; ++i)
        st_probs[i] *= V_reff_xr_rp3;
}

__attribute__ ((always_inline))
inline real bg_probability(__constant ASTRONOMY_PARAMETERS* ap,
                           __constant STREAM_CONSTANTS* sc,
                           __global const R_POINTS* r_pts,
                           const LB_TRIG lbt,
                           const R_CONSTS rc,
                           const real V,
                           real* st_probs)
{
    real bg_prob;

    /* if q is 0, there is no probability */
    if (ap->zero_q)
        return -1.0;

    bg_prob = sub_bg_probability(ap, sc, r_pts, rc.gPrime, lbt, ap->convolve, st_probs);
    bg_prob *= rc.reff_xr_rp3;

    mult_probs(st_probs, V * rc.reff_xr_rp3, ap->number_streams);

    return bg_prob;
}

__attribute__ ((always_inline))
real r_calculation(__constant ASTRONOMY_PARAMETERS* ap,
                   __constant INTEGRAL_AREA* ia,
                   __constant STREAM_CONSTANTS* sc,
                   __constant R_CONSTS* rcs,
                   __global const R_POINTS* r_pts,
                   const LB_TRIG lbt,
                   const real id,
                   real* st_probs,
                   const unsigned int r_step)
{
    real V;
    R_CONSTS rc;

    rc = rcs[r_step];
    V = id * rc.irv;

    return V * bg_probability(ap, sc, &r_pts[r_step * ap->convolve], lbt, rc, V, st_probs);
}

__attribute__ ((always_inline))
inline void write_st_probs(__global real* probs_out,
                           const real* st_probs,
                           const unsigned int n_streams)
{
    unsigned int i;

    for (i = 0; i < n_streams; ++i)
        probs_out[i] = st_probs[i];
}

__kernel void mu_sum_kernel(__global real* mu_out,
                            __global real* probs_out,

                            __constant ASTRONOMY_PARAMETERS* ap,
                            __constant INTEGRAL_AREA* ia,
                            __constant STREAM_CONSTANTS* sc,
                            __constant R_CONSTS* rcs,
                            __global const R_POINTS* r_pts,
                            const unsigned int nu_step)
{
    NU_ID nuid;
    LB lb;
    LB_TRIG lbt;
    real mu, r_result;

    real st_probs[3];     /* FIXME: hardcoded stream limit */

    size_t idx;     /* index into stream probs output buffer */
    size_t mu_step = get_global_id(0);
    size_t r_step = get_global_id(1);


    if (mu_step > ia->mu_steps || r_step > ia->r_steps)
        return;

    zero_st_probs(st_probs, ap->number_streams);

    /* Actual calculations */
    mu = ia->mu_min + (((real) mu_step + 0.5) * ia->mu_step_size);
    nuid = calc_nu_step(ia, nu_step);
    lb = gc2lb(ap->wedge, mu, nuid.nu);
    lbt = lb_trig(lb);

    r_result = r_calculation(ap, ia, sc, rcs, r_pts, lbt, nuid.id, st_probs, r_step);

    /* Output to buffers */
    mu_out[mu_step * ia->r_steps + r_step] = r_result;
    idx = mu_step * ia->r_steps * ap->number_streams + r_step * ap->number_streams;
    write_st_probs(&probs_out[idx],
                   st_probs,
                   ap->number_streams);

}


