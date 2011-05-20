/*
 * Copyright (c) 2008-2010 Travis Desell, Nathan Cole, Boleslaw
 * Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
 * Rensselaer Polytechnic Institute.
 * Copyright (c) 2010-2011 Matthew Arsenault
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
 *
 *
 * This file incorporates work covered by the following copyright and
 * permission notice:
 *
 *      Fast exp(x) computation (with SSE2 optimizations).
 *
 * Copyright (c) 2010, Naoaki Okazaki
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the names of the authors nor the names of its contributors
 *       may be used to endorse or promote products derived from this
 *       software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "evaluation_state.h"
#include "evaluation.h"
#include "coordinates.h"
#include "r_points.h"
#include "milkyway_util.h"
#include "calculated_constants.h"
#include "probabilities.h"

#include <time.h>

#ifdef __SSE3__
  #include <pmmintrin.h>
#elif __SSE2__
  #include <emmintrin.h>
#endif /* __SSE3__ */

#if defined(__GNUC__)
  #pragma GCC diagnostic ignored "-Wunknown-pragmas"
#elif defined(_MSC_VER)
#pragma warning( disable : 4068 )
#endif /* __GNUC__ */



#if defined(__SSE2__) && DOUBLEPREC && !ENABLE_CRLIBM
  #define USE_CUSTOM_MATH 1
#endif


#if defined(__SSE3__)
  #define INIT_PROBABILITIES initProbabilities_SSE3
#elif defined(__SSE2__)
  #define INIT_PROBABILITIES initProbabilities_SSE2
#else
  #define INIT_PROBABILITIES initProbabilities
#endif


#ifdef __SSE2__

typedef   signed char      ssp_s8;
typedef unsigned char      ssp_u8;
typedef   signed short     ssp_s16;
typedef unsigned short     ssp_u16;
typedef   signed int       ssp_s32;
typedef unsigned int       ssp_u32;
typedef float              ssp_f32;
typedef double             ssp_f64;
typedef   signed long long ssp_s64;
typedef unsigned long long ssp_u64;

typedef union
{
    __m128  f;
    __m128d d;
    __m128i i;
    __m64   m64[2];
    ssp_u64 u64[2];
    ssp_s64 s64[2];
    ssp_f64 f64[2];
    ssp_u32 u32[4];
    ssp_s32 s32[4];
    ssp_f32 f32[4];
    ssp_u16 u16[8];
    ssp_s16 s16[8];
    ssp_u8  u8 [16];
    ssp_s8  s8 [16];
} ssp_m128;
#endif /* __SSE2__ */




#ifdef __SSE2__

#define EXP2_TABLE_SIZE_LOG2 12
#define EXP2_TABLE_SIZE (1 << EXP2_TABLE_SIZE_LOG2)
#define EXP2_TABLE_OFFSET (EXP2_TABLE_SIZE/2)
#define EXP2_TABLE_SCALE ((double) ((EXP2_TABLE_SIZE/2)-1))
#define TBL_SIZE_RECIP ((double)(1/EXP2_TABLE_SIZE))

/* 2^x, for x in [-1.0, 1.0[ */
SEPARATION_ALIGN(16) static double exp2_table[2 * EXP2_TABLE_SIZE];

static void initExpTable()
{
    unsigned int i;

    #pragma ivdep
    #pragma vector always
    for (i = 0; i < EXP2_TABLE_SIZE; ++i)
    {
        exp2_table[i] = (double) mw_exp2((i - EXP2_TABLE_OFFSET) / EXP2_TABLE_SCALE);
    }
}

union d2
{
    int i[4];
    unsigned int u[4];
    long int lu[2];
    double d[2];
    __m128d m;
    __m128i mi;
};

#define CONST_128D(var, val)                                    \
    SEPARATION_ALIGN(16) static const double var[2] = {(val), (val)}
#define CONST_128I(var, v1, v2, v3, v4)                                 \
    SEPARATION_ALIGN(16) static const int var[4] = {(v1), (v2), (v3), (v4)}

CONST_128D(log2e, 1.4426950408889634073599);
CONST_128D(maxlog, 7.09782712893383996843e2);   // log(2**1024)
CONST_128D(minlog, -7.08396418532264106224e2);  // log(2**-1022)
CONST_128D(c1, 6.93145751953125E-1);
CONST_128D(c2, 1.42860682030941723212E-6);
CONST_128D(w11, 3.5524625185478232665958141148891055719216674475023e-8);
CONST_128D(w10, 2.5535368519306500343384723775435166753084614063349e-7);
CONST_128D(w9, 2.77750562801295315877005242757916081614772210463065e-6);
CONST_128D(w8, 2.47868893393199945541176652007657202642495832996107e-5);
CONST_128D(w7, 1.98419213985637881240770890090795533564573406893163e-4);
CONST_128D(w6, 1.3888869684178659239014256260881685824525255547326e-3);
CONST_128D(w5, 8.3333337052009872221152811550156335074160546333973e-3);
CONST_128D(w4, 4.1666666621080810610346717440523105184720007971655e-2);
CONST_128D(w3, 0.166666666669960803484477734308515404418108830469798);
CONST_128D(w2, 0.499999999999877094481580370323249951329122224389189);
CONST_128D(w1, 1.0000000000000017952745258419615282194236357388884);
CONST_128D(w0, 0.99999999999999999566016490920259318691496540598896);

static const double MAXLOG =  7.08396418532264106224E2;     /* log 2**1022 */
static const double MINLOG = -7.08396418532264106224E2;     /* log 2**-1022 */
static const double LOG2E  =  1.4426950408889634073599;     /* 1/log(2) */
//static const double INFINITY = 1.79769313486231570815E308;
static const double C1 = 6.93145751953125E-1;
static const double C2 = 1.42860682030941723212E-6;


#define DONEHALF    (_mm_set_pd(0.5, 0.5))
#define DONE        (_mm_set_pd(1.0, 1.0))
#define DONE_NEG    (_mm_set_pd(-1.0,-1.0))
#define DTWO        (_mm_set_pd(2.0, 2.0))
#define DTHREE      (_mm_set_pd(3.0, 3.0))
#define DFOUR       (_mm_set_pd(4.0, 4.0))
#define DFIVE       (_mm_set_pd(5.0, 5.0))
#define LOG2ME      (_mm_set_pd(1.44269504089, 1.44269504089))


#define Vneg(x) _mm_sub_pd(_mm_setzero_pd(), (x))

static inline __m128d hadd_pd_SSE2(__m128d a, __m128d b)
{
    ssp_m128 A, B, C;
    A.d = a;
    C.d = a;
    B.d = b;
    A.f = _mm_movelh_ps(A.f, B.f);
    B.f = _mm_movehl_ps(B.f, C.f);
    A.d = _mm_add_pd(A.d, B.d);
    return A.d;
}

static inline __m128d _mm_fexpd_lut_pd(__m128d x)
{
    __m128i ipart;
    __m128d fpart, expipart, xmm0, y;
    union d2 idx, expfpart;
    const __m128i offset = _mm_setr_epi32(1023, 1023, 0, 0);

    y     = _mm_min_pd(x, _mm_set1_pd(1023.00000));
    y     = _mm_max_pd(x, _mm_set1_pd(-1022.99999));
    xmm0  = _mm_load_pd(log2e);
    y     = _mm_mul_pd(y, xmm0);
    ipart = _mm_cvtpd_epi32(y);
    fpart = _mm_sub_pd(y, _mm_cvtepi32_pd(ipart));

    /* remez */
    expipart = _mm_castsi128_pd(_mm_shuffle_epi32(_mm_slli_epi32(_mm_add_epi32(ipart, offset), 23), _MM_SHUFFLE(1, 3, 0, 2)));
    idx.mi = _mm_add_epi32(_mm_cvtpd_epi32(_mm_mul_pd(fpart, _mm_set1_pd(EXP2_TABLE_SCALE))), _mm_set1_epi32(EXP2_TABLE_OFFSET));

    expfpart.d[0] = exp2_table[idx.u[0]];
    expfpart.d[1] = exp2_table[idx.u[1]];
    /*
      expfpart.d[2] = exp2_table[idx.u[2]];
      expfpart.d[3] = exp2_table[idx.u[3]];
    */

    return _mm_mul_pd(expipart, expfpart.m);
}


static inline __m128d _mm_fexp_poly_pd(__m128d x1) // precise,but slower than _mm_exp_pd
{
    /* remez11_0_log2_sse */
    __m128i k1;
    __m128d p1, a1;
    __m128d xmm0, xmm1;

    const __m128i offset = _mm_setr_epi32(1023, 1023, 0, 0);
    x1 = _mm_min_pd(x1, _mm_set1_pd( 129.00000));
    x1 = _mm_max_pd(x1, _mm_set1_pd(-126.99999));

    xmm0 = _mm_load_pd(log2e);
    xmm1 = _mm_setzero_pd();
    a1 = _mm_mul_pd(x1, xmm0);
    /* k = (int)floor(a); p = (float)k; */
    p1 = _mm_cmplt_pd(a1, xmm1);
    p1 = _mm_and_pd(p1, DONE);
    a1 = _mm_sub_pd(a1, p1);
    k1 = _mm_cvttpd_epi32(a1); // ipart
    p1 = _mm_cvtepi32_pd(k1);
    /* x -= p * log2; */
    xmm0 = _mm_load_pd(c1);
    xmm1 = _mm_load_pd(c2);
    a1 = _mm_mul_pd(p1, xmm0);
    x1 = _mm_sub_pd(x1, a1);
    a1 = _mm_mul_pd(p1, xmm1);
    x1 = _mm_sub_pd(x1, a1);
    /* Compute e^x using a polynomial approximation. */

    xmm0 = _mm_load_pd(w11);
    xmm1 = _mm_load_pd(w10);

    a1 = _mm_mul_pd(x1, xmm0);
    a1 = _mm_add_pd(a1, xmm1);
    xmm0 = _mm_load_pd(w9);
    xmm1 = _mm_load_pd(w8);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm0);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm1);

    xmm0 = _mm_load_pd(w7);
    xmm1 = _mm_load_pd(w6);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm0);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm1);

    xmm0 = _mm_load_pd(w5);
    xmm1 = _mm_load_pd(w4);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm0);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm1);

    xmm0 = _mm_load_pd(w3);
    xmm1 = _mm_load_pd(w2);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm0);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm1);

    xmm0 = _mm_load_pd(w1);
    xmm1 = _mm_load_pd(w0);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm0);
    a1 = _mm_mul_pd(a1, x1);
    a1 = _mm_add_pd(a1, xmm1);
    /* p = 2^k; */
    k1 = _mm_add_epi32(k1, offset);
    k1 = _mm_slli_epi32(k1, 20);
    k1 = _mm_shuffle_epi32(k1, _MM_SHUFFLE(1, 3, 0, 2));
    p1 = _mm_castsi128_pd(k1);

    /* a *= 2^k. */
    return  _mm_mul_pd(a1, p1);
}

static inline __m128d _mm_muladd_pd(__m128d x, __m128d y, __m128d z)
{
    return _mm_add_pd(_mm_mul_pd(x, y), z);
}

static real probabilities_SSE2(const AstronomyParameters* ap,
                               const StreamConstants* sc,
                               const real* RESTRICT sg_dx,
                               const real* RESTRICT r_point,
                               const real* RESTRICT qw_r3_N,
                               LBTrig lbt,
                               real gPrime,
                               real reff_xr_rp3,
                               real* RESTRICT streamTmps)
{
    double bg_prob, dotted, xyz_norm;
    unsigned int i, j, k, convolve, nStreams;
    SEPARATION_ALIGN(16) double psgt[1024], psgf[1024], xyzstr[1024];
    SEPARATION_ALIGN(16) mwvector xyz[256];

    __m128d REF_XR = _mm_set1_pd(reff_xr_rp3);

    //__m128d ALPHAV = _mm_set1_pd(alpha);
    __m128d COSBL    = _mm_set1_pd(lbt.lCosBCos);
    __m128d SINB     = _mm_set1_pd(lbt.bSin);
    __m128d SINCOSBL = _mm_set1_pd(lbt.lSinBCos);
    __m128d SUNR0    = _mm_set1_pd(ap->sun_r0);
    __m128d R0       = _mm_set1_pd(ap->r0);
    __m128d QV_RECIP = _mm_set1_pd(ap->q_inv);
    __m128d RI, QI;
    ssp_m128 xyz0, xyz1, xyz2, tmp0, tmp1, PROD, PBXV, BGP;

    BGP.d = _mm_setzero_pd();

    convolve = ap->convolve;
    nStreams = ap->number_streams;

    #pragma ivdep
    #pragma vector always
    for (i = 0; i < convolve; i += 2)
    {
        RI = _mm_load_pd(&r_point[i]);
        QI = _mm_load_pd(&qw_r3_N[i]);

        xyz0.d = _mm_sub_pd(_mm_mul_pd(RI, COSBL), SUNR0);

        _mm_storel_pd(&X(xyz[i]), xyz0.d);
        _mm_storeh_pd(&X(xyz[i + 1]), xyz0.d);

        xyz1.d = _mm_mul_pd(RI, SINCOSBL);
        _mm_storel_pd(&Y(xyz[i]), xyz1.d);
        _mm_storeh_pd(&Y(xyz[i + 1]), xyz1.d);

        xyz2.d = _mm_mul_pd(RI, SINB);
        tmp0.d = _mm_mul_pd(xyz2.d, QV_RECIP);
        _mm_storel_pd(&Z(xyz[i]), xyz2.d);
        _mm_storeh_pd(&Z(xyz[i + 1]), xyz2.d);

        xyz0.d = _mm_mul_pd(xyz0.d, xyz0.d);
        xyz1.d = _mm_mul_pd(xyz1.d, xyz1.d);
        tmp0.d = _mm_mul_pd(tmp0.d, tmp0.d);
     // PROD.d = _mm_fsqrt_pd(_mm_add_pd(xyz0.d, _mm_add_pd(xyz1.d, tmp0.d))); // slower
        PROD.d = _mm_sqrt_pd(_mm_add_pd(xyz0.d, _mm_add_pd(xyz1.d, tmp0.d)));
        tmp1.d = _mm_add_pd(PROD.d, R0);
     // PBXV.d = _mm_rcp_pd(_mm_mul_pd(PROD.d, _mm_mul_pd(tmp1.d, _mm_mul_pd(tmp1.d,tmp1.d))));  // slower
        PBXV.d = _mm_div_pd(DONE, _mm_mul_pd(PROD.d, _mm_mul_pd(tmp1.d, _mm_mul_pd(tmp1.d, tmp1.d))));
        BGP.d  = _mm_add_pd(BGP.d, _mm_mul_pd(QI, PBXV.d));
    }

    BGP.d = _mm_mul_pd(BGP.d, REF_XR);

  #ifdef __SSE3__
    BGP.d = _mm_hadd_pd(BGP.d, BGP.d);
  #else // SSE2 emu hadd_pd
    #pragma message("USING SSE2 HADD")
    BGP.d = hadd_pd_SSE2(BGP.d, BGP.d);
  #endif
    _mm_store_sd(&bg_prob, BGP.d);

    #pragma ivdep
    #pragma vector always
    for (i = 0; i < nStreams; ++i)
    {
        #pragma ivdep
        #pragma vector always
        for (j = 0; j < convolve; ++j)
        {
            xyzstr[0] = X(xyz[j]) - X(sc[i].c);
            xyzstr[1] = Y(xyz[j]) - Y(sc[i].c);
            xyzstr[2] = Z(xyz[j]) - Z(sc[i].c);

            dotted = X(sc[i].a) * xyzstr[0] + Y(sc[i].a) * xyzstr[1] + Z(sc[i].a) * xyzstr[2];

            xyzstr[0] = xyzstr[0] - dotted * X(sc[i].a);
            xyzstr[1] = xyzstr[1] - dotted * Y(sc[i].a);
            xyzstr[2] = xyzstr[2] - dotted * Z(sc[i].a);

            xyz_norm = xyzstr[0] * xyzstr[0] + xyzstr[1] * xyzstr[1] + xyzstr[2] * xyzstr[2];
            psgt[j] = xyz_norm * sc[i].sigma_sq2_inv;
        }

        for (k = 0; k < convolve; k += 2)
        {
            //_mm_store_pd(&psgf[k], _mm_mul_pd(_mm_load_pd(&qw_r3_N[k]), _mm_fexpd_lut_pd(Vneg(_mm_load_pd(&psgt[k])))));

            // arst
            //_mm_store_pd(&psgf[k], _mm_mul_pd(_mm_load_pd(&qw_r3_N[k]), _mm_exp_pd(Vneg(_mm_load_pd(&psgt[k])))));

            _mm_store_pd(&psgf[k], _mm_mul_pd(_mm_load_pd(&qw_r3_N[k]), _mm_fexp_poly_pd(Vneg(_mm_load_pd(&psgt[k])))));
        }

        streamTmps[i] = 0.0;
        #pragma ivdep
        #pragma vector always
        for (k = 0; k < convolve; ++k)
        {
            streamTmps[i] += psgf[k];
        }
        streamTmps[i] *= reff_xr_rp3;
    }

    return bg_prob;
}

#endif /* __SSE2__ */


static inline mwvector lbr2xyz_2(const AstronomyParameters* ap, real rPoint, LBTrig lbt)
{
    mwvector xyz;

    xyz.x = mw_mad(rPoint, LCOS_BCOS(lbt), ap->m_sun_r0);
    xyz.y = rPoint * LSIN_BCOS(lbt);
    xyz.z = rPoint * BSIN(lbt);

    return xyz;
}

static inline real streamIncrement(const StreamConstants* sc, mwvector xyz)
{
    real xyz_norm, dotted;
    mwvector xyzs;

    xyzs = mw_subv(xyz, sc->c);
    dotted = mw_dotv(sc->a, xyzs);
    mw_incsubv_s(xyzs, sc->a, dotted);

    xyz_norm = mw_sqrv(xyzs);

    return mw_exp(-xyz_norm * sc->sigma_sq2_inv);
}

HOT
static inline void streamSums(real* st_probs,
                              const StreamConstants* sc,
                              const mwvector xyz,
                              const real qw_r3_N,
                              const unsigned int nstreams)
{
    unsigned int i;

    for (i = 0; i < nstreams; ++i)
        st_probs[i] += qw_r3_N * streamIncrement(&sc[i], xyz);
}

HOT
static inline real h_prob_fast(const AstronomyParameters* ap, real qw_r3_N, real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (rg * cube(rs));
}

HOT
static inline real h_prob_slow(const AstronomyParameters* ap, real qw_r3_N, real rg)
{
    const real rs = rg + ap->r0;
    return qw_r3_N / (mw_powr(rg, ap->alpha) * mw_powr(rs, ap->alpha_delta3));
}

HOT
static inline real rg_calc(const AstronomyParameters* ap, const mwvector xyz)
{
    /* sqrt(x^2 + y^2 + q_inv_sqr * z^2) */

    real tmp;

    tmp = sqr(X(xyz));
    tmp = mw_mad(Y(xyz), Y(xyz), tmp);               /* x^2 + y^2 */
    tmp = mw_mad(ap->q_inv_sqr, sqr(Z(xyz)), tmp);   /* (q_invsqr * z^2) + (x^2 + y^2) */

    return mw_sqrt(tmp);
}

static inline void zero_st_probs(real* st_probs, const unsigned int nstream)
{
    unsigned int i;

    for (i = 0; i < nstream; ++i)
        st_probs[i] = 0.0;
}

static inline real aux_prob(const AstronomyParameters* ap,
                            const real qw_r3_N,
                            const real r_in_mag)
{
    return qw_r3_N * (ap->bg_a * sqr(r_in_mag) + ap->bg_b * r_in_mag + ap->bg_c);
}

HOT
static real bg_probability_fast_hprob(const AstronomyParameters* ap,
                                      const StreamConstants* sc,
                                      const real* RESTRICT sg_dx,
                                      const real* RESTRICT r_point,
                                      const real* RESTRICT qw_r3_N,
                                      LBTrig lbt,
                                      real gPrime,
                                      real reff_xr_rp3,
                                      real* RESTRICT streamTmps)
{
    unsigned int i;
    real h_prob, g, rg;
    mwvector xyz;
    real bg_prob = 0.0;
    unsigned int convolve = ap->convolve;
    int aux_bg_profile = ap->aux_bg_profile;

    zero_st_probs(streamTmps, ap->number_streams);
    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_point[i], lbt);
        rg = rg_calc(ap, xyz);

        h_prob = h_prob_fast(ap, qw_r3_N[i], rg);

        /* Add a quadratic term in g to the the Hernquist profile */
        if (aux_bg_profile)
        {
            g = gPrime + sg_dx[i];
            h_prob += aux_prob(ap, qw_r3_N[i], g);
        }

        bg_prob += h_prob;
        streamSums(streamTmps, sc, xyz, qw_r3_N[i], ap->number_streams);
    }

    bg_prob *= reff_xr_rp3;
    for (i = 0; i < ap->number_streams; ++i)
        streamTmps[i] *= reff_xr_rp3;

    return bg_prob;
}

HOT
static real bg_probability_slow_hprob(const AstronomyParameters* ap,
                                      const StreamConstants* sc,
                                      const real* RESTRICT sg_dx,
                                      const real* RESTRICT r_point,
                                      const real* RESTRICT qw_r3_N,
                                      LBTrig lbt,
                                      real gPrime,
                                      real reff_xr_rp3,
                                      real* RESTRICT streamTmps)
{
    unsigned int i;
    real rg, g;
    mwvector xyz;
    real bg_prob = 0.0;
    unsigned int convolve = ap->convolve;
    int aux_bg_profile = ap->aux_bg_profile;

    zero_st_probs(streamTmps, ap->number_streams);

    for (i = 0; i < convolve; ++i)
    {
        xyz = lbr2xyz_2(ap, r_point[i], lbt);

        rg = rg_calc(ap, xyz);

        bg_prob += h_prob_slow(ap, qw_r3_N[i], rg);
        if (aux_bg_profile)
        {
            g = gPrime + sg_dx[i];
            bg_prob += aux_prob(ap, qw_r3_N[i], g);
        }

        streamSums(streamTmps, sc, xyz, qw_r3_N[i], ap->number_streams);
    }

    bg_prob *= reff_xr_rp3;
    for (i = 0; i < ap->number_streams; ++i)
        streamTmps[i] *= reff_xr_rp3;

    return bg_prob;
}


/* We have some more deciding on which function to use for whatever SSE level */

ProbabilityFunc INIT_PROBABILITIES(const AstronomyParameters* ap, int forceNoIntrinsics)
{
  #ifndef __SSE2__
    return ap->fast_h_prob ? bg_probability_fast_hprob : bg_probability_slow_hprob;
  #else

    initExpTable();

    if (ap->fast_h_prob && !ap->aux_bg_profile && !forceNoIntrinsics && USE_CUSTOM_MATH)
    {
        /* TODO: add aux_bg_profile in SSE2 stuff */
        return probabilities_SSE2;
    }
    else if (ap->fast_h_prob)
    {
        return bg_probability_fast_hprob;
    }
    else
    {
        return bg_probability_slow_hprob;
    }

    mw_unreachable();
    return NULL;
  #endif /* __SSE2__ */
}

