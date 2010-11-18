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
#include "milkyway_cl.h"
#include "milkyway_extra.h"
#include "integrals_common.h"

#if FAST_H_PROB
  #define h_prob_f h_prob_fast
#else
  #define h_prob_f h_prob_slow
#endif /* FAST_H_PROB */

/* Compiler apparently can't unroll the loops over the streams if we
 * loop over NStreams (as of Stream SDK 2.2), so we have to manually
 * expand. A macro to expand an arbitrary number of times is also
 * quite terrible. Not unrolling these loops murders the
 * performance. */
#if NSTREAM > 5
  #error "Requested number of streams greater than supported maximum (5)"
#elif NSTREAM <= 0
  #error "Requested number of streams is absurd"
#endif


#if USE_IMAGES

#if 0
/* This breaks Nvidia compiler */
RPoints readImageDouble(uint4 a)
{
    union {
        uint4 i;
        RPoints d;
    } arst;
    arst.i = a;
    return arst.d;
}
#endif

inline RPoints readImageDouble(uint4 a)
{
    union {
        uint2 i[2];
        RPoints d;
    } arst;
    arst.i[0] = a.lo;
    arst.i[1] = a.hi;
    return arst.d;
}

inline RPoints readRPts(__read_only image2d_t r_pts, __constant IntegralArea* ia, int2 i)
{
    const sampler_t sample = CLK_ADDRESS_NONE
                           | CLK_NORMALIZED_COORDS_FALSE
                           | CLK_FILTER_NEAREST;

    return readImageDouble(read_imageui(r_pts, sample, i));
}

#else

inline RPoints readRPts(__global const RPoints* r_pts, __constant IntegralArea* ia, int2 i)
{
    return r_pts[i.y * ia->r_steps + i.x];
}

#endif /* USE_IMAGES */

__attribute__ ((always_inline))
inline void stream_sums_cl(real* st_probs,
                           __constant StreamConstants* sc,
                           const mwvector xyz,
                           const RPoints r_pt)
{
  #if NSTREAM >= 1
    st_probs[0] = mw_mad(QW_R3_N(r_pt), calc_st_prob_inc(sc[0], xyz), st_probs[0]);
  #endif

  #if NSTREAM >= 2
    st_probs[1] = mw_mad(QW_R3_N(r_pt), calc_st_prob_inc(sc[1], xyz), st_probs[1]);
  #endif

  #if NSTREAM >= 3
    st_probs[2] = mw_mad(QW_R3_N(r_pt), calc_st_prob_inc(sc[2], xyz), st_probs[2]);
  #endif

  #if NSTREAM >= 4
    st_probs[3] = mw_mad(QW_R3_N(r_pt), calc_st_prob_inc(sc[3], xyz), st_probs[3]);
  #endif

  #if NSTREAM >= 5
    st_probs[4] = mw_mad(QW_R3_N(r_pt), calc_st_prob_inc(sc[4], xyz), st_probs[4]);
  #endif
}

__attribute__ ((always_inline))
inline void write_mult_st_probs(__global real* probs_out, real V_reff_xr_rp3, const real* st_probs)
{
  #if NSTREAM >= 1
    probs_out[0] += V_reff_xr_rp3 * st_probs[0];
  #endif

  #if NSTREAM >=2
    probs_out[1] += V_reff_xr_rp3 * st_probs[1];
  #endif

  #if NSTREAM >= 3
    probs_out[2] += V_reff_xr_rp3 * st_probs[2];
  #endif

  #if NSTREAM >= 4
    probs_out[3] += V_reff_xr_rp3 * st_probs[3];
  #endif

  #if NSTREAM >= 5
     probs_out[4] += V_reff_xr_rp3 * st_probs[4];
  #endif
}


#define MAX_CONST(n, type) __attribute__((max_constant_size(n * sizeof(type))))

__kernel void mu_sum_kernel(__global real* restrict mu_out,
                            __global real* restrict probs_out,

                            __constant AstronomyParameters* ap MAX_CONST(1, AstronomyParameters),
                            __constant IntegralArea* ia MAX_CONST(1, IntegralArea),
                            __constant StreamConstants* sc MAX_CONST(NSTREAM, StreamConstants),
                            __constant RConsts* rcs MAX_CONST(200, RConsts),
                            __constant real* restrict sg_dx MAX_CONST(200, real),

                          #if USE_IMAGES
                            __read_only image2d_t r_pts,
                          #else
                            __global const RPoints* r_pts,
                          #endif

                            __global const LBTrig* lbts,
                            const real nu_id)
{
    size_t nu_step = get_global_id(2);
    size_t mu_step = get_global_id(1);
    size_t r_step  = get_global_id(0);

    LBTrig lbt = lbts[nu_step * ia->mu_steps + mu_step]; /* 32-byte read */

    real bg_prob = 0.0;
    real st_probs[NSTREAM] = { 0.0 };

    unsigned int i;
    unsigned int convolve = ap->convolve; /* Faster to load this into register first */
    for (i = 0; i < convolve; ++i)
    {
        RPoints r_pt = readRPts(r_pts, ia, (int2) (r_step, i));

        mwvector xyz = lbr2xyz_2(ap, R_POINT(r_pt), lbt);
        real rg = rg_calc(ap, xyz);
        bg_prob += h_prob_f(ap, QW_R3_N(r_pt), rg);

        stream_sums_cl(st_probs, sc, xyz, r_pt);

        /* Moving stream_sums up from here reduces GPR usage by 2, but also
         * for some reason gets slightly slower. */

        /* Using stream_sums_cl twice (while giving a nonsense result)
         * somehow ends up using fewer registers */

      #if AUX_BG_PROFILE
        /* Currently not used */
        /* Add a quadratic term in g to the Hernquist profile */
        real g = GPRIME(rc[r_step]) + sg_dx[i];
        bg_prob += aux_prob(ap, QW_R3_N(r_pt), g);
      #endif /* AUX_BG_PROFILE */
    }

    real V_reff_xr_rp3 = nu_id * IRV_REFF_XR_RP3(rcs[r_step]);
    size_t idx = mu_step * ia->r_steps + r_step; /* Index into output buffers */

    bg_prob *= V_reff_xr_rp3;
    mu_out[idx] += bg_prob;

    write_mult_st_probs(&probs_out[NSTREAM * idx], V_reff_xr_rp3, st_probs);

}


