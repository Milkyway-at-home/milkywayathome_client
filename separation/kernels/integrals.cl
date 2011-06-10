/*
Copyright 2010, 2011 Matthew Arsenault, Travis Desell, Dave Przybylo,
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

#include "separation_kernel_types.h"


#define cube(x) ((x) * (x) * (x))
#define sqr(x) ((x) * (x))

#define USE_CUSTOM_SQRT 1
#define USE_CUSTOM_DIVISION 1

#if !defined(__Cypress__) && !defined(__ATI_RV770__) && !defined(__CPU__)
  #define USE_CUSTOM_DIVISION 1

  /* It's faster to load the stream constants into private memory on
   * Nvidia, which I think there spills into shared stuff. It seems
   * much slower on ATI though. */
  #define LOAD_STREAM_CONSTANTS 1
#endif


#if LOAD_STREAM_CONSTANTS
  #define __SC_CONSTANT __private
#else
  #define __SC_CONSTANT __constant
#endif /* LOAD_STREAM_CONSTANTS */

#if USE_IMAGES

const sampler_t sample = CLK_ADDRESS_NONE
                       | CLK_NORMALIZED_COORDS_FALSE
                       | CLK_FILTER_NEAREST;

inline RPoints readRPts(__read_only image2d_t r_pts, int2 i)
{
    return as_double2((read_imageui(r_pts, sample, i)));
}

#else

inline RPoints readRPts(__global const RPoints* r_pts, int2 i)
{
    return r_pts[i.x * CONVOLVE + i.y];
}

#endif /* USE_IMAGES */

#if LOAD_STREAM_CONSTANTS

inline void set_sc_priv(StreamConstants* sc, __constant StreamConstants* sc_c)
{
    unsigned int i;

    #pragma unroll NSTREAM
    for (i = 0; i < NSTREAM; ++i)
    {
        sc[i] = sc_c[i];
    }
}

#endif /* LOAD_STREAM_CONSTANTS */

#if USE_CUSTOM_DIVISION && DOUBLEPREC

/* TODO: Move these */

double mw_div(double a, double b)  // accurate to 1 ulp, i.e the last bit of the double precision number
{
    // cuts some corners on the numbers range but is significantly faster, employs "faithful rounding"
    double r;
    double y;
    double c;
    y = (double)(1.0f / convert_float_rtn(b)); // 22bit estimate of reciprocal, limits range to float range, but spares the exponent extraction
    c = 1.0 - b * y;
    y += y * c;        // first Newton iteration => 44bit accurate
    r = a * y;         // second iteration works directly on a/b, so it is effectively
    c = a - b * r;     // a Newton-Markstein iteration without guaranteed round to nearest
    return(r + y * c); // should be generally accurate to 1 ulp, i.e "faithful rounding" (or close to it, depending on definition)
}                      // on a GT200 should be round to nearest most of the time, albeit not guaranteed
                       // one would have to add a second Newton iteration before the Markstein rounding step,
                       // but one looses the gained half bit of precision in the following additions, so the added effort doesn't make sense

#else
  #define mw_div(a, b) ((a) / (b))
#endif /* USE_CUSTOM_DIVISION && DOUBLEPREC */

#if USE_CUSTOM_SQRT && DOUBLEPREC

double mw_fsqrt(double y)  // accurate to 1 ulp, i.e the last bit of the double precision number
{
    // cuts some corners on the numbers range but is significantly faster, employs "faithful rounding"
    double x, res;
    x = (double)rsqrt((float)y); // 22bit estimate for reciprocal square root, limits range to float range, but spares the exponent extraction
    x = x * (3.0 - y * (x * x)); // first Newton iteration (44bit accurate)
    res = x * y;  // do final iteration directly on sqrt(y) and not on the inverse
    return(res * (0.75 - 0.0625 * (res * x)));
}   // same precision as division (1 ulp)

#else
  #define mw_fsqrt sqrt
#endif /* USE_CUSTOM_SQRT && DOUBLEPREC */

#if 0
  #define MAX_CONST(n, type) __attribute__((max_constant_size(n * sizeof(type))))
#else
  #define MAX_CONST(n, type)
#endif

#if LOAD_STREAM_CONSTANTS
  #define SC_ARG sc_c
#else
  #define SC_ARG sc
#endif /* LOAD_STREAM_CONSTANTS */


inline real aux_prob(__constant AstronomyParameters* ap,
                     const real qw_r3_N,
                     const real r_in_mag)
{
    real tmp;

    tmp = mad(ap->bg_b, r_in_mag, ap->bg_c); /* bg_b * r_in_mag + bg_c */
    tmp = mad(ap->bg_a, sqr(r_in_mag), tmp); /* bg_a * r_in_mag2 + (bg_b * r_in_mag + bg_c)*/

    return qw_r3_N * tmp;
}


__kernel void mu_sum_kernel(__global real* restrict bgOut,
                            __global real* restrict streamsOut,

                            __constant AstronomyParameters* ap MAX_CONST(1, AstronomyParameters),
                            __constant IntegralArea* ia MAX_CONST(1, IntegralArea),

                            __constant StreamConstants* SC_ARG MAX_CONST(NSTREAM, StreamConstants),
                            __constant RConsts* rcs,
                            __constant real* restrict sg_dx MAX_CONST(256, real),

                          #if USE_IMAGES
                            __read_only image2d_t r_pts,
                          #else
                            __global const __read_only RPoints* r_pts,
                          #endif

                            __global const __read_only LBTrig* lbts,
                            const unsigned int extra,
                            const real nu_id)
{
    size_t nu_step = get_global_id(1);
    size_t mu_step = (get_global_id(0) - extra) % ia->mu_steps;
    size_t r_step  = (get_global_id(0) - extra) / ia->mu_steps;

    if (r_step >= ia->r_steps || mu_step >= ia->mu_steps) /* Avoid out of bounds from roundup */
        return;

    LBTrig lbt = lbts[nu_step * ia->mu_steps + mu_step]; /* 32-byte read */

    real bg_prob = 0.0;
    real st_probs[NSTREAM] = { 0.0 };

  #if LOAD_STREAM_CONSTANTS
    StreamConstants sc[NSTREAM];
    set_sc_priv(sc, sc_c);
  #endif /* LOAD_STREAM_CONSTANTS */

    unsigned int i, j;

    for (i = 0; i < CONVOLVE; ++i)
    {
        RPoints r_pt = readRPts(r_pts, (int2) (r_step, i));

        real x = mad(R_POINT(r_pt), LCOS_BCOS(lbt), (real) -SUN_R0);
        real y = R_POINT(r_pt) * LSIN_BCOS(lbt);
        real z = R_POINT(r_pt) * BSIN(lbt);

        /* sqrt(x^2 + y^2 + q_inv_sqr * z^2) */
        real tmp = x * x;
        tmp = mad(y, y, tmp);           /* x^2 + y^2 */
        tmp = mad((real) Q_INV_SQR, z * z, tmp);   /* (q_invsqr * z^2) + (x^2 + y^2) */

        real rg = mw_fsqrt(tmp);
        real rs = rg + R0;

      #if FAST_H_PROB
        bg_prob += mw_div(QW_R3_N(r_pt), (rg * cube(rs)));
      #else
        bg_prob += mw_div(qw_r3_N, mw_powr(rg, ap->alpha) * mw_powr(rs, ap->alpha_delta3));
      #endif /* FAST_H_PROB */

      #if AUX_BG_PROFILE
        /* Currently not used */
        /* Add a quadratic term in g to the Hernquist profile */
        real g = GPRIME(rcs[r_step]) + sg_dx[i];
        bg_prob += aux_prob(ap, QW_R3_N(r_pt), g);
      #endif /* AUX_BG_PROFILE */

        #pragma unroll NSTREAM
        for (j = 0; j < NSTREAM; ++j)
        {
            real xs = x - X(sc[j].c);
            real ys = y - Y(sc[j].c);
            real zs = z - Z(sc[j].c);

            real dotted = X(sc[j].a) * xs;
            dotted = mad(Y(sc[j].a), ys, dotted);
            dotted = mad(Z(sc[j].a), zs, dotted);

            xs = mad(dotted, -X(sc[j].a), xs);
            ys = mad(dotted, -Y(sc[j].a), ys);
            zs = mad(dotted, -Z(sc[j].a), zs);

            real sqrv = xs * xs;
            sqrv = mad(ys, ys, sqrv);
            sqrv = mad(zs, zs, sqrv);

            real tmp = exp(-sqrv * sc[j].sigma_sq2_inv);

            st_probs[j] = mad(QW_R3_N(r_pt), tmp, st_probs[j]);
        }
    }

    real V_reff_xr_rp3 = nu_id * IRV_REFF_XR_RP3(rcs[r_step]);
    size_t idx = mu_step * ia->r_steps + r_step; /* Index into output buffers */

    bg_prob *= V_reff_xr_rp3;
    bgOut[idx] += bg_prob;

    #pragma unroll NSTREAM
    for (j = 0; j < NSTREAM; ++j)
    {
        streamsOut[NSTREAM * idx + j] += V_reff_xr_rp3 * st_probs[j];
    }
}

