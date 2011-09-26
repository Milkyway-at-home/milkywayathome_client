/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newbergb
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 1991-2000 University of Groningen, The Netherlands
 *  Copyright (c) 2001-2009 The GROMACS Development Team
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
 *  This file incorporates work covered by the following copyright and
 *  permission notice:
 *
 *    Fast exp(x) computation (with SSE2 optimizations).
 *
 *  Copyright (c) 2010, Naoaki Okazaki
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * Neither the names of the authors nor the names of its contributors
 *        may be used to endorse or promote products derived from this
 *        software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 *  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _PROBABILITIES_INTRIN_H_
#define _PROBABILITIES_INTRIN_H_


#if defined(__SSE4_1__)
  #include <smmintrin.h>
#elif defined(__SSE3__)
  #include <pmmintrin.h>
#elif defined(__SSE2__)
  #include <emmintrin.h>
#endif /* __SSE3__ */

#if defined(__GNUC__)
  #pragma GCC diagnostic ignored "-Wunknown-pragmas"
#elif defined(_MSC_VER)
  #pragma warning( disable : 4068 )
#endif /* __GNUC__ */




#define SSE2_EXP_LUT 0

#if defined(__SSE2__) && DOUBLEPREC && !ENABLE_CRLIBM
  #define USE_CUSTOM_MATH 1
#endif



#ifdef __cplusplus
extern "C" {
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

#define EXP2_TABLE_SIZE_LOG2 12
#define EXP2_TABLE_SIZE (1 << EXP2_TABLE_SIZE_LOG2)
#define EXP2_TABLE_OFFSET (EXP2_TABLE_SIZE/2)
#define EXP2_TABLE_SCALE ((double) ((EXP2_TABLE_SIZE/2)-1))
#define TBL_SIZE_RECIP ((double)(1/EXP2_TABLE_SIZE))

/* 2^x, for x in [-1.0, 1.0[ */
SEPARATION_ALIGN(16) static double exp2_table[2 * EXP2_TABLE_SIZE];

static inline void initExpTable()
{
    int i;

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
static inline __m128d gmx_mm_exp_pd(__m128d x)
{
    const __m128d argscale = _mm_set1_pd(1.442695040888963387);
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m128d arglimit = _mm_set1_pd(-1022.0/1.442695040888963387);
    const __m128i expbase  = _mm_set1_epi32(1023);

    const __m128d CA0       = _mm_set1_pd(7.0372789822689374920e-9);
    const __m128d CB0       = _mm_set1_pd(89.491964762085371);
    const __m128d CB1       = _mm_set1_pd(-9.7373870675164587);
    const __m128d CC0       = _mm_set1_pd(51.247261867992408);
    const __m128d CC1       = _mm_set1_pd(-0.184020268133945);
    const __m128d CD0       = _mm_set1_pd(36.82070153762337);
    const __m128d CD1       = _mm_set1_pd(5.416849282638991);
    const __m128d CE0       = _mm_set1_pd(30.34003452248759);
    const __m128d CE1       = _mm_set1_pd(8.726173289493301);
    const __m128d CF0       = _mm_set1_pd(27.73526969472330);
    const __m128d CF1       = _mm_set1_pd(10.284755658866532);

    __m128d valuemask;
    __m128i iexppart;
    __m128d fexppart;
    __m128d intpart;
    __m128d z,z2;
    __m128d factB,factC,factD,factE,factF;

      z       = _mm_mul_pd(x,argscale);
    iexppart  = _mm_cvtpd_epi32(z);

#if defined __SSE4_1__
    /* This reduces latency and speeds up the code by roughly 5% when supported */
    intpart   = _mm_round_pd(z,0);
#else
    intpart   = _mm_cvtepi32_pd(iexppart);
#endif
    /* The two lowest elements of iexppart now contains 32-bit numbers with a correctly biased exponent.
     * To be able to shift it into the exponent for a double precision number we first need to
     * shuffle so that the lower half contains the first element, and the upper half the second.
     * This should really be done as a zero-extension, but since the next instructions will shift
     * the registers left by 52 bits it doesn't matter what we put there - it will be shifted out.
     * (thus we just use element 2 from iexppart).
     */
    iexppart  = _mm_shuffle_epi32(iexppart,_MM_SHUFFLE(2,1,2,0));

    /* Do the shift operation on the 64-bit registers */
    iexppart  = _mm_add_epi32(iexppart,expbase);
    iexppart  = _mm_slli_epi64(iexppart,52);
    valuemask = _mm_cmpgt_pd(x,arglimit);

      z         = _mm_sub_pd(z,intpart);
    z2        = _mm_mul_pd(z,z);

    fexppart  = _mm_and_pd(valuemask,_mm_castsi128_pd(iexppart));

    /* Since SSE doubleing-point has relatively high latency it is faster to do
     * factorized polynomial summation with independent terms than using alternating add/multiply, i.e.
     * p(z) = A0 * (B0 + z) * (C0 + C1*z + z^2) * (D0 + D1*z + z^2) * (E0 + E1*z + z^2) * (F0 + F1*z + z^2)
     */

    factB     = _mm_add_pd(CB0,_mm_mul_pd(CB1,z) );
    factB     = _mm_add_pd(factB,z2);
    factC     = _mm_add_pd(CC0,_mm_mul_pd(CC1,z) );
    factC     = _mm_add_pd(factC,z2);
    factD     = _mm_add_pd(CD0,_mm_mul_pd(CD1,z) );
    factD     = _mm_add_pd(factD,z2);
    factE     = _mm_add_pd(CE0,_mm_mul_pd(CE1,z) );
    factE     = _mm_add_pd(factE,z2);
    factF     = _mm_add_pd(CF0,_mm_mul_pd(CF1,z) );
    factF     = _mm_add_pd(factF,z2);

    z         = _mm_mul_pd(CA0,fexppart);
    factB     = _mm_mul_pd(factB,factC);
    factD     = _mm_mul_pd(factD,factE);
    z         = _mm_mul_pd(z,factF);
    factB     = _mm_mul_pd(factB,factD);
    z         = _mm_mul_pd(z,factB);

      /* Currently uses 32 actual (real, not including casts) SSE instructions */
      return  z;
}
static inline __m128d _mm_invsqrt_pd(__m128d x)
{
    const __m128d half  = _mm_set1_pd(0.5);
    const __m128d three = _mm_set1_pd(3.0);

    /* Lookup instruction only exists in single precision, convert back and forth... */
    __m128d xmm0 = _mm_cvtps_pd(_mm_rsqrt_ps( _mm_cvtpd_ps(x)));

    xmm0 = _mm_mul_pd(half,_mm_mul_pd(_mm_sub_pd(three,_mm_mul_pd(_mm_mul_pd(xmm0,xmm0),x)),xmm0));
    return _mm_mul_pd(half,_mm_mul_pd(_mm_sub_pd(three,_mm_mul_pd(_mm_mul_pd(xmm0,xmm0),x)),xmm0));
}

static inline __m128d _mm_fsqrt_pd(__m128d y)  // accurate to 1 ulp, i.e the last bit of the double precision number
{
    // cuts some corners on the numbers range but is significantly faster, employs "faithful rounding"
    __m128d x, res;
	const __m128d C0  = _mm_set1_pd(0.75);
	const __m128d C1  = _mm_set1_pd(0.0625);
    x = _mm_cvtps_pd(_mm_rsqrt_ps( _mm_cvtpd_ps(y))); // 22bit estimate for reciprocal square root, limits range to float range, but spares the exponent extraction
    x = _mm_mul_pd(x,_mm_sub_pd(DTHREE,_mm_mul_pd(y,_mm_mul_pd(x,x)))); // first Newton iteration (44bit accurate)
    res = _mm_mul_pd(x, y);  // do final iteration directly on sqrt(y) and not on the inverse
    return(_mm_mul_pd(res,_mm_sub_pd(C0,_mm_mul_pd(C1,_mm_mul_pd(res , x)))));
}

static inline __m128d _mm_fdiv_pd(__m128d a, __m128d b) // fasted on P4 compared to _mm_div_pd !!!
{
    __m128d r,y,c;

    y = _mm_cvtps_pd(_mm_rcp_ps( _mm_cvtpd_ps(b)));
    c = _mm_sub_pd(DONE,_mm_mul_pd(b,y));
    y = _mm_add_pd(y,_mm_mul_pd(y,c));
    r = _mm_mul_pd(a,y);
    c = _mm_sub_pd(a,_mm_mul_pd(b,r));
    return _mm_add_pd(r,_mm_mul_pd(y,c));
}

static inline __m128d _mm_rcp_pd(__m128d x)
{
    __m128d xmm0 = _mm_cvtps_pd(_mm_rcp_ps( _mm_cvtpd_ps(x)));
    /* Newton Raphson */
    xmm0 = _mm_mul_pd(xmm0,_mm_sub_pd(DTWO,_mm_mul_pd(x,xmm0)));
    return _mm_mul_pd(xmm0,_mm_sub_pd(DTWO,_mm_mul_pd(x,xmm0)));
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

#endif /* __SSE2__ */


#ifdef __cplusplus
}
#endif

#endif /* _PROBABILITIES_INTRIN_H_ */

