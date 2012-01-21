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

#if defined(__AVX__)
  #include <immintrin.h>
#elif defined(__SSE4_1__)
  #include <smmintrin.h>
#elif defined(__SSE3__)
  #include <pmmintrin.h>
#elif defined(__SSE2__)
  #include <emmintrin.h>
#else
  #error No intrinsic enabled
#endif /* __AVX__ */

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
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

#ifdef __AVX__
typedef union
{
    __m256  f;
    __m256d d;
    __m256i i;

    __m128 m128[2];
    __m128 f128[2];
    __m128d d128[2];
   __m128i i128[2];

    __m64   m64[4];
    ssp_u64 u64[4];
    ssp_s64 s64[4];
    ssp_f64 f64[4];
    ssp_u32 u32[8];
    ssp_s32 s32[8];
    ssp_f32 f32[8];
    ssp_u16 u16[16];
    ssp_s16 s16[16];
    ssp_u8  u8 [32];
    ssp_s8  s8 [32];
} ssp_m256;
#endif /* __AVX__ */


#define EXP2_TABLE_SIZE_LOG2 12
#define EXP2_TABLE_SIZE (1 << EXP2_TABLE_SIZE_LOG2)
#define EXP2_TABLE_OFFSET (EXP2_TABLE_SIZE/2)
#define EXP2_TABLE_SCALE ((double) ((EXP2_TABLE_SIZE/2)-1))
#define TBL_SIZE_RECIP ((double)(1/EXP2_TABLE_SIZE))

/* 2^x, for x in [-1.0, 1.0[ */
MW_ALIGN_V(32) static double exp2_table[2 * EXP2_TABLE_SIZE];

static inline void initExpTable(void)
{
    int i;

    #pragma ivdep
    #pragma vector always
    for (i = 0; i < EXP2_TABLE_SIZE; ++i)
    {
        exp2_table[i] = (double) mw_exp2((i - EXP2_TABLE_OFFSET) / EXP2_TABLE_SCALE);
    }
}

//static const double MAXLOG =  7.08396418532264106224E2;     /* log 2**1022 */
//static const double MINLOG = -7.08396418532264106224E2;     /* log 2**-1022 */
//static const double LOG2E  =  1.4426950408889634073599;     /* 1/log(2) */
//static const double INFINITY = 1.79769313486231570815E308;
//static const double C1 = 6.93145751953125E-1;
//static const double C2 = 1.42860682030941723212E-6;

#define DONE256        (_mm256_set_pd(1.0, 1.0, 1.0, 1.0))
#define DTHREE256      (_mm256_set_pd(3.0, 3.0, 3.0, 3.0))

//#define DONEHALF    (_mm_set_pd(0.5, 0.5))
#define DONE        (_mm_set_pd(1.0, 1.0))
//#define DONE_NEG    (_mm_set_pd(-1.0,-1.0))
#define DTWO        (_mm_set_pd(2.0, 2.0))
#define DTHREE      (_mm_set_pd(3.0, 3.0))
//#define DFOUR       (_mm_set_pd(4.0, 4.0))
//#define DFIVE       (_mm_set_pd(5.0, 5.0))
//#define LOG2ME      (_mm_set_pd(1.44269504089, 1.44269504089))


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

    z         = _mm_mul_pd(x,argscale);
    iexppart  = _mm_cvtpd_epi32(z);

  #ifdef __SSE4_1__
    /* This reduces latency and speeds up the code by roughly 5% when supported */
    intpart   = _mm_round_pd(z,0);
  #else
    intpart   = _mm_cvtepi32_pd(iexppart);
  #endif /* __SSE4_1__ */

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



#ifdef __cplusplus
}
#endif

#endif /* _PROBABILITIES_INTRIN_H_ */

