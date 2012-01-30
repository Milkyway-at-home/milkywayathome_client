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

#ifndef DOUBLEPREC
  #error DOUBLEPREC not defined
#endif

#if DOUBLEPREC
  #if defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
  #elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
  #else
    #error No double extension available
  #endif /* defined(cl_khr_fp64) */

  #if defined(__AMD__) || defined(__GPU__)
    /* Older AMD versions did not define __AMD__ but did define __GPU__ */
    #if defined(__ATI_RV770__) || defined(__Cypress__) || defined(__Cayman__)
      #define USE_CUSTOM_DIVISION 1
      #define USE_CUSTOM_SQRT 1
    #elif defined(__Tahiti__)
    /* Using the default division is a bit faster, but using ddiv is
     * slower than the custom division in the IL kernel?  */
      #define USE_CUSTOM_DIVISION 0

      /* This is just barely faster */
      #define USE_CUSTOM_SQRT 1
    #else
      /* Using neither is a safe fallback */
      #define USE_CUSTOM_DIVISION 0
      #define USE_CUSTOM_SQRT 0
    #endif
  #elif defined(__CPU__)
    /* AMD OpenCL's CPU implementation */
    #define USE_CUSTOM_DIVISION 0
    #define USE_CUSTOM_SQRT 0
  #else /* Not AMD OpenCL */
    /* I can't find anything about things Nvidia defines and reasonable
     * guesses don't work, so do what's best for Nvidia */
    #define USE_CUSTOM_DIVISION 1
    #define USE_CUSTOM_SQRT 1
  #endif /* defined(__AMD__) || defined(__GPU__) */
#endif /* DOUBLEPREC */


#define MAX_CONVOLVE 256


#if DOUBLEPREC
typedef double real;
typedef double2 real2;
typedef double4 real4;
#else
typedef float real;
typedef float2 real2;
typedef float4 real4;
#endif /* DOUBLEPREC */


#define cube(x) ((x) * (x) * (x))
#define sqr(x) ((x) * (x))


#if DOUBLEPREC
double mw_div_custom(double a, double b)
{
    // accurate to 1 ulp, i.e the last bit of the double precision number
    // cuts some corners on the numbers range but is significantly faster, employs "faithful rounding"
    // 22bit estimate of reciprocal, limits range to float range, but spares the exponent extraction

    double y = (double)(1.0f / convert_float_rtn(b));
    double c = mad(-b, y, 1.0);
    y = mad(y, c, y);         // first Newton iteration => 44bit accurate
    double r = a * y;         // second iteration works directly on a/b, so it is effectively
    c = mad(-b, r, a);        // a Newton-Markstein iteration without guaranteed round to nearest

    return mad(y, c, r);
    // should be generally accurate to 1 ulp, i.e "faithful rounding" (or close to it, depending on definition)
    // on a GT200 should be round to nearest most of the time, albeit not guaranteed
    // one would have to add a second Newton iteration before the Markstein rounding step,
    // but one looses the gained half bit of precision in the following additions, so the added effort doesn't make sense
}

double mw_sqrt_custom(double y)  // accurate to 1 ulp, i.e the last bit of the double precision number
{
    // cuts some corners on the numbers range but is significantly faster, employs "faithful rounding"
    // 22bit estimate for reciprocal square root, limits range to float range, but spares the exponent extraction
    double x = (double)rsqrt((float)y);
    x = x * mad(-y, x * x, 3.0);         // first Newton iteration (44bit accurate)
    double res = x * y;                  // do final iteration directly on sqrt(y) and not on the inverse
    return res * mad(-0.0625, res * x, 0.75);
    // same precision as division (1 ulp)
}
#endif /* DOUBLEPREC */


#if USE_CUSTOM_DIVISION
  #define mw_div(a, b) mw_div_custom((a), (b))
#else
  #define mw_div(a, b) ((a) / (b))
#endif /* USE_CUSTOM_DIVISION */

#if USE_CUSTOM_SQRT
  #define mw_fsqrt mw_sqrt_custom
#else
  #define mw_fsqrt sqrt
#endif /* USE_CUSTOM_SQRT */


inline real aux_prob(real r_in_mag)
{
    real tmp;

    tmp = mad(BG_B, r_in_mag, BG_C); /* bg_b * r_in_mag + bg_c */
    tmp = mad(BG_A, sqr(r_in_mag), tmp); /* bg_a * r_in_mag2 + (bg_b * r_in_mag + bg_c)*/

    return tmp;
}

typedef struct
{
    real x_a, x_c;
    real y_a, y_c;
    real z_a, z_c;
    real sigma_sq2_inv;
    real _pad;
} SC;

inline real2 kahanSum(real2 running, real term)
{
    real correctedNextTerm = term + running.y;
    real newSum = running.x + correctedNextTerm;
    running.y = correctedNextTerm - (newSum - running.x);
    running.x = newSum;
    return running;
}

/* Be careful of changing the arguments. It will most likely break the IL kernel */
__kernel void probabilities(__global real2* restrict bgOut,
                            __global real2* restrict streamsOut,

                            __global const real2* restrict rConsts,
                            __global const real2* restrict rPts,
                            __global const real2* restrict lTrigBuf,
                            __global const real* restrict bSinBuf,


                            /* Placeholder for IL kernel */
                            __constant real* _ap_consts __attribute__((max_constant_size(16 * sizeof(real)))),

                            __constant SC sc[NSTREAM] __attribute__((max_constant_size(NSTREAM * sizeof(SC)))),
                            __constant real sg_dx[CONVOLVE] __attribute__((max_constant_size(256 * sizeof(real)))),

                            const unsigned int extra,
                            const unsigned int r_steps,
                            const unsigned int mu_steps,
                            const unsigned int nu_steps,
                            const real nu_id,
                            const unsigned int nu_step)
{
    size_t gid = get_global_id(0) - extra;
    size_t mu_step = gid % mu_steps;
    size_t r_step  = gid / mu_steps;

    if (r_step >= r_steps || mu_step >= mu_steps) /* Avoid out of bounds from roundup */
        return;

    size_t idx = mu_step * r_steps + r_step; /* Index into output buffers */

    size_t trigIdx = nu_step * mu_steps + mu_step;
    real2 lTrig = lTrigBuf[trigIdx];
    real bSin = bSinBuf[trigIdx];
    real2 rc = rConsts[r_step];


    real2 bgTmp;
    real2 streamTmp[NSTREAM];

    bgTmp = bgOut[idx];
    #pragma unroll NSTREAM
    for (int j = 0; j < NSTREAM; ++j)
    {
        streamTmp[j] = streamsOut[j * mu_steps * r_steps + idx];
    }

    real bg_prob = 0.0;
    real st_probs[NSTREAM] = { 0.0 };

    for (int i = 0; i < CONVOLVE; ++i)
    {
        real2 rPt = rPts[CONVOLVE * r_step + i];

        real x = mad(rPt.x, lTrig.x, (real) -SUN_R0);
        real y = rPt.x * lTrig.y;
        real z = rPt.x * bSin;

        /* sqrt(x^2 + y^2 + q_inv_sqr * z^2) */
        real tmp = x * x;
        tmp = mad(y, y, tmp);           /* x^2 + y^2 */
        tmp = mad((real) Q_INV_SQR, z * z, tmp);   /* (q_invsqr * z^2) + (x^2 + y^2) */

        real rg = mw_fsqrt(tmp);
        real rs = rg + R0;

        if (FAST_H_PROB)
        {
            bg_prob += mw_div(rPt.y, rg * cube(rs));
        }
        else
        {
            bg_prob += mw_div(rPt.y, powr(rg, ALPHA) * powr(rs, ALPHA_DELTA_3));
        }

        if (AUX_BG_PROFILE)
        {
            /* Currently not used */
            /* Add a quadratic term in g to the Hernquist profile */
            real g = rc.y + sg_dx[i];
            bg_prob = mad(rPt.y, aux_prob(g), bg_prob);
        }

        #pragma unroll NSTREAM
        for (int j = 0; j < NSTREAM; ++j)
        {
            real xs = x - sc[j].x_c;
            real ys = y - sc[j].y_c;
            real zs = z - sc[j].z_c;

            real dotted = sc[j].x_a * xs;
            dotted = mad(sc[j].y_a, ys, dotted);
            dotted = mad(sc[j].z_a, zs, dotted);

            xs = mad(dotted, -sc[j].x_a, xs);
            ys = mad(dotted, -sc[j].y_a, ys);
            zs = mad(dotted, -sc[j].z_a, zs);

            real sqrv = xs * xs;
            sqrv = mad(ys, ys, sqrv);
            sqrv = mad(zs, zs, sqrv);

            real tmp = exp(-sqrv * sc[j].sigma_sq2_inv);

            st_probs[j] = mad(rPt.y, tmp, st_probs[j]);
        }
    }

    real V_reff_xr_rp3 = nu_id * rc.x;
    //uint nElement = mu_steps * r_steps;

    bg_prob *= V_reff_xr_rp3;
    #pragma unroll NSTREAM
    for (uint j = 0; j < NSTREAM; ++j)
    {
        st_probs[j] *= V_reff_xr_rp3;
    }

    bgTmp = kahanSum(bgTmp, bg_prob);
    #pragma unroll NSTREAM
    for (uint j = 0; j < NSTREAM; ++j)
    {
        streamTmp[j] = kahanSum(streamTmp[j], st_probs[j]);
    }


    bgOut[idx] = bgTmp;
    #pragma unroll NSTREAM
    for (uint j = 0; j < NSTREAM; ++j)
    {
        streamsOut[j * mu_steps * r_steps + idx] = streamTmp[j];
    }
}

