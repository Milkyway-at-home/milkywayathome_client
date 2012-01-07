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

#ifndef __AVX__
  #error AVX not defined
#endif

#include "milkyway_util.h"
#include "probabilities.h"
#include "probabilities_intrin.h"

#define Vneg(x) _mm256_sub_pd(_mm256_setzero_pd(), (x))

#if 0
static inline __m256d gmx_mm256_exp_pd(__m256d x)
{
    const __m256d argscale = _mm256_set1_pd(1.442695040888963387);
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m256d arglimit = _mm256_set1_pd(-1022.0/1.442695040888963387);
    const __m128i expbase  = _mm_set1_epi32(1023);

    const __m256d CA0       = _mm256_set1_pd(7.0372789822689374920e-9);
    const __m256d CB0       = _mm256_set1_pd(89.491964762085371);
    const __m256d CB1       = _mm256_set1_pd(-9.7373870675164587);
    const __m256d CC0       = _mm256_set1_pd(51.247261867992408);
    const __m256d CC1       = _mm256_set1_pd(-0.184020268133945);
    const __m256d CD0       = _mm256_set1_pd(36.82070153762337);
    const __m256d CD1       = _mm256_set1_pd(5.416849282638991);
    const __m256d CE0       = _mm256_set1_pd(30.34003452248759);
    const __m256d CE1       = _mm256_set1_pd(8.726173289493301);
    const __m256d CF0       = _mm256_set1_pd(27.73526969472330);
    const __m256d CF1       = _mm256_set1_pd(10.284755658866532);

    __m256d valuemask;
    __m256d fexppart;
    __m256d intpart;
    __m256d z, z2;
    __m256d factB, factC, factD, factE, factF;

    //__m128i iexppart;
    union
    {
        __m128i i128[2];
        __m128d d128[2];
        __m256d d;
        __m256 m;
        __m256i i;
    } iexppart;

    z          = _mm256_mul_pd(x, argscale);
    iexppart.i128[0] = _mm256_cvtpd_epi32(z);

    intpart = _mm256_round_pd(z, 0);


    /* The four lowest elements of iexppart now contains 32-bit numbers with a correctly biased exponent.
     * To be able to shift it into the exponent for a double precision number we first need to
     * shuffle so that the lower half contains the first element, and the upper half the second.
     * This should really be done as a zero-extension, but since the next instructions will shift
     * the registers left by 52 bits it doesn't matter what we put there - it will be shifted out.
     * (thus we just use element 2 from iexppart).
     */

  #ifdef __AVX2__
    // untested since I don't think anything exists with this yet
    iexppart.m = _mm256_shuffle_epi32(iexppart.m, _MM_SHUFFLE(2, 1, 2, 0));
  #else
    iexppart.m = _mm256_shuffle_ps(iexppart.m, iexppart.m, _MM_SHUFFLE(2, 1, 2, 0));
  #endif


    /* Do the shift operation on the 64-bit registers */
    iexppart.i128[0]  = _mm_add_epi32(iexppart.i128[0], expbase);
    iexppart.i128[1]  = _mm_slli_epi64(iexppart.i128[1], 52);


    valuemask = _mm256_cmpgt_pd(x, arglimit);

    z         = _mm256_sub_pd(z, intpart);
    z2        = _mm256_mul_pd(z, z);

    fexppart  = _mm256_and_pd(valuemask, _mm256_castsi128_pd(iexppart));

    /* Since SSE doubleing-point has relatively high latency it is faster to do
     * factorized polynomial summation with independent terms than using alternating add/multiply, i.e.
     * p(z) = A0 * (B0 + z) * (C0 + C1*z + z^2) * (D0 + D1*z + z^2) * (E0 + E1*z + z^2) * (F0 + F1*z + z^2)
     */

    factB     = _mm256_add_pd(CB0, _mm256_mul_pd(CB1,z));
    factB     = _mm256_add_pd(factB,z2);
    factC     = _mm256_add_pd(CC0, _mm256_mul_pd(CC1,z));
    factC     = _mm256_add_pd(factC,z2);
    factD     = _mm256_add_pd(CD0, _mm256_mul_pd(CD1,z));
    factD     = _mm256_add_pd(factD,z2);
    factE     = _mm256_add_pd(CE0, _mm256_mul_pd(CE1,z));
    factE     = _mm256_add_pd(factE,z2);
    factF     = _mm256_add_pd(CF0, _mm256_mul_pd(CF1,z));
    factF     = _mm256_add_pd(factF,z2);

    z         = _mm256_mul_pd(CA0, fexppart);
    factB     = _mm256_mul_pd(factB, factC);
    factD     = _mm256_mul_pd(factD, factE);
    z         = _mm256_mul_pd(z, factF);
    factB     = _mm256_mul_pd(factB, factD);
    z         = _mm256_mul_pd(z, factB);

    return  z;
}
#else

static inline __m256d gmx_mm256_exp_pd(__m256d x)
{
    union
    {
        __m256d d;
        __m128d d128[2];
    } tmp;

    tmp.d = x;
    tmp.d128[0] = gmx_mm_exp_pd(tmp.d128[0]);
    tmp.d128[1] = gmx_mm_exp_pd(tmp.d128[1]);
    return tmp.d;
}

#endif

static real probabilities_avx(const AstronomyParameters* ap,
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
    int i, j, k, convolve, nStreams;
    MW_ALIGN_V(64) double psgt[256], psgf[256], xyzstr[256];
    MW_ALIGN_V(64) double xs[256], ys[256], zs[256];

    const __m256d REF_XR = _mm256_set1_pd(reff_xr_rp3);

    const __m256d COSBL    = _mm256_set1_pd(lbt.lCosBCos);
    const __m256d SINB     = _mm256_set1_pd(lbt.bSin);
    const __m256d SINCOSBL = _mm256_set1_pd(lbt.lSinBCos);
    const __m256d SUNR0    = _mm256_set1_pd(ap->sun_r0);
    const __m256d R0       = _mm256_set1_pd(ap->r0);
    const __m256d QV_RECIP = _mm256_set1_pd(ap->q_inv);
    __m256d RI, QI;
    ssp_m256 xyz0, xyz1, xyz2, tmp0, tmp1, tmp2, PROD, PBXV, BGP;

    BGP.d = _mm256_setzero_pd();

    convolve = ap->convolve;
    nStreams = ap->number_streams;

    for (i = 0; i < convolve; i += 4)
    {
        RI = _mm256_load_pd(&r_point[i]);
        QI = _mm256_load_pd(&qw_r3_N[i]);

        xyz0.d = _mm256_sub_pd(_mm256_mul_pd(RI, COSBL), SUNR0);
        //xyz0.d = _mm256_fmadd_pd(RI, COSBL, NSUNR0);

        _mm256_store_pd(&xs[i], xyz0.d);

        xyz1.d = _mm256_mul_pd(RI, SINCOSBL);
        _mm256_store_pd(&ys[i], xyz1.d);

        xyz2.d = _mm256_mul_pd(RI, SINB);
        tmp0.d = _mm256_mul_pd(xyz2.d, QV_RECIP);

        _mm256_store_pd(&zs[i], xyz2.d);

        xyz0.d = _mm256_mul_pd(xyz0.d, xyz0.d);
        xyz1.d = _mm256_mul_pd(xyz1.d, xyz1.d);
        tmp0.d = _mm256_mul_pd(tmp0.d, tmp0.d);

        tmp1.d = _mm256_add_pd(xyz0.d, _mm256_add_pd(xyz1.d, tmp0.d));

        PROD.d = _mm256_sqrt_pd(tmp1.d);
        tmp2.d = _mm256_add_pd(PROD.d, R0);

        PBXV.d = _mm256_div_pd(DONE256, _mm256_mul_pd(PROD.d, _mm256_mul_pd(tmp2.d, _mm256_mul_pd(tmp2.d, tmp2.d))));
        BGP.d  = _mm256_add_pd(BGP.d, _mm256_mul_pd(QI, PBXV.d));
    }

    BGP.d = _mm256_mul_pd(BGP.d, REF_XR);
    BGP.d = _mm256_hadd_pd(BGP.d, BGP.d);
    BGP.d128[0] = _mm_add_pd(BGP.d128[0], BGP.d128[1]);
    _mm_store_sd(&bg_prob, BGP.d128[0]);

    for (i = 0; i < nStreams; ++i)
    {
        for (j = 0; j < convolve; ++j)
        {
            xyzstr[0] = xs[j] - X(sc[i].c);
            xyzstr[1] = ys[j] - Y(sc[i].c);
            xyzstr[2] = zs[j] - Z(sc[i].c);

            dotted = X(sc[i].a) * xyzstr[0] + Y(sc[i].a) * xyzstr[1] + Z(sc[i].a) * xyzstr[2];

            xyzstr[0] = xyzstr[0] - dotted * X(sc[i].a);
            xyzstr[1] = xyzstr[1] - dotted * Y(sc[i].a);
            xyzstr[2] = xyzstr[2] - dotted * Z(sc[i].a);

            xyz_norm = xyzstr[0] * xyzstr[0] + xyzstr[1] * xyzstr[1] + xyzstr[2] * xyzstr[2];
            psgt[j] = xyz_norm * sc[i].sigma_sq2_inv;
        }

        for (k = 0; k < convolve; k += 4)
        {
            /* Other maths from SSE2 one when those are figured out */
            #pragma message("Using GMX Polynomal EXP")
            _mm256_store_pd(&psgf[k], _mm256_mul_pd(_mm256_load_pd(&qw_r3_N[k]), gmx_mm256_exp_pd(Vneg(_mm256_load_pd(&psgt[k])))));
		}

        streamTmps[i] = 0.0;
        #pragma ivdep
        #pragma vector always
        for (j = 0; j < convolve; ++j)
        {
            streamTmps[i] += psgf[j];
        }
        streamTmps[i] *= reff_xr_rp3;
    }

    return bg_prob;
}

ProbabilityFunc INIT_PROBABILITIES(void)
{
    assert(mwAllocA32Safe());
    initExpTable();
    return probabilities_avx;
}


