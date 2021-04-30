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
static inline __m256d gmx_mm256_exp2_pd(__m256d x)
{
    /* Lower bound: We do not allow numbers that would lead to an IEEE fp representation exponent smaller than -126. */
    const __m256d arglimit = _mm256_set1_pd(1022.0);
    const __m128i expbase  = _mm_set1_epi32(1023);

    const __m256d P2       = _mm256_set1_pd(2.30933477057345225087e-2);
    const __m256d P1       = _mm256_set1_pd(2.02020656693165307700e1);
    const __m256d P0       = _mm256_set1_pd(1.51390680115615096133e3);
    /* Q2 == 1.0 */
    const __m256d Q1       = _mm256_set1_pd(2.33184211722314911771e2);
    const __m256d Q0       = _mm256_set1_pd(4.36821166879210612817e3);
    const __m256d one      = _mm256_set1_pd(1.0);
    const __m256d two      = _mm256_set1_pd(2.0);

    __m256d       valuemask;
    __m256i       iexppart;
    __m128i       iexppart128a, iexppart128b;
    __m256d       fexppart;
    __m256d       intpart;
    __m256d       z, z2;
    __m256d       PolyP, PolyQ;

    iexppart128a  = _mm256_cvtpd_epi32(x);
    intpart       = _mm256_round_pd(x, _MM_FROUND_TO_NEAREST_INT);

    /* Add exponent bias */
    iexppart128a   = _mm_add_epi32(iexppart128a, expbase);

    /* We now want to shift the exponent 52 positions left, but to achieve this we need
     * to separate the 128-bit register data into two registers (4x64-bit > 128bit)
     * shift them, and then merge into a single __m256d.
     * Elements 0/1 should end up in iexppart128a, and 2/3 in iexppart128b.
     * It doesnt matter what we put in the 2nd/4th position, since that data will be
     * shifted out and replaced with zeros.
     */
    iexppart128b   = _mm_shuffle_epi32(iexppart128a, _MM_SHUFFLE(3, 3, 2, 2));
    iexppart128a   = _mm_shuffle_epi32(iexppart128a, _MM_SHUFFLE(1, 1, 0, 0));

    iexppart128b   = _mm_slli_epi64(iexppart128b, 52);
    iexppart128a   = _mm_slli_epi64(iexppart128a, 52);

    iexppart  = _mm256_castsi128_si256(iexppart128a);
    iexppart  = _mm256_insertf128_si256(iexppart, iexppart128b, 0x1);

    valuemask = _mm256_cmp_pd(arglimit, gmx_mm256_abs_pd(x), _CMP_GE_OQ);
    fexppart  = _mm256_and_pd(valuemask, _mm256_castsi256_pd(iexppart));

    z         = _mm256_sub_pd(x, intpart);

    z2        = _mm256_mul_pd(z, z);

    PolyP     = _mm256_mul_pd(P2, z2);
    PolyP     = _mm256_add_pd(PolyP, P1);
    PolyQ     = _mm256_add_pd(z2, Q1);
    PolyP     = _mm256_mul_pd(PolyP, z2);
    PolyQ     = _mm256_mul_pd(PolyQ, z2);
    PolyP     = _mm256_add_pd(PolyP, P0);
    PolyQ     = _mm256_add_pd(PolyQ, Q0);
    PolyP     = _mm256_mul_pd(PolyP, z);

    z         = _mm256_mul_pd(PolyP, gmx_mm256_inv_pd(_mm256_sub_pd(PolyQ, PolyP)));
    z         = _mm256_add_pd(one, _mm256_mul_pd(two, z));

    z         = _mm256_mul_pd(z, fexppart);

    return z;
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

/******************************************************************************************
 ********************** Hernquist Profile Function ****************************************
 ******************************************************************************************/

static real probabilities_avx_hernquist(const AstronomyParameters* ap,
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
    const __m256d THICKLS  = _mm256_set1_pd(-0.285714286);
    const __m256d THICKHS  = _mm256_set1_pd(-1.428571429);
    /*Calculate halo and disk coeffients
      Can probably write this in a better way */    
    const __m256d THICKCOEF   = _mm256_set1_pd(ap->thick_disk_weight);
    const __m256d BGCOEF      = _mm256_set1_pd(ap->background_weight);


    __m256d RI, QI;
    ssp_m256 xyz0, xyz1, xyz2, tmp0, tmp1, tmp2, PROD, PBXV, BGP;
    //xyz0, 1, 2 = x, y, z
    BGP.d = _mm256_setzero_pd();

    convolve = ap->convolve;
    nStreams = ap->number_streams;

    for (i = 0; i < convolve; i += 4)
    {
    	/* Put r_point and qw_r3_n into RI and QI respectively */
        RI = _mm256_load_pd(&r_point[i]);
        QI = _mm256_load_pd(&qw_r3_N[i]);

        /* Coordinate Transform to Galactic Center XYZ */
        xyz0.d = _mm256_sub_pd(_mm256_mul_pd(RI, COSBL), SUNR0); //X Value
        /* xyz0.d = _mm256_fmadd_pd(RI, COSBL, NSUNR0); */

        _mm256_store_pd(&xs[i], xyz0.d);

        xyz1.d = _mm256_mul_pd(RI, SINCOSBL); /* Y Value */
        _mm256_store_pd(&ys[i], xyz1.d);

        xyz2.d = _mm256_mul_pd(RI, SINB); /* Z Value */
        tmp0.d = _mm256_mul_pd(xyz2.d, QV_RECIP); /* Squashed Z */

        _mm256_store_pd(&zs[i], xyz2.d);

        //Finish Convert to Galactic Center XYZ
        xyz0.d = _mm256_mul_pd(xyz0.d, xyz0.d); /* X^2 */
        xyz1.d = _mm256_mul_pd(xyz1.d, xyz1.d); /* Y^2 */
        tmp0.d = _mm256_mul_pd(tmp0.d, tmp0.d); /* (Z/q)^2 Squishy Squishy */

        tmp1.d = _mm256_add_pd(xyz0.d, _mm256_add_pd(xyz1.d, tmp0.d)); /* Calculate R^2 */

        PROD.d = _mm256_sqrt_pd(tmp1.d); /* Calculate R from R^2 */
        tmp2.d = _mm256_add_pd(PROD.d, R0); /* Calculate R + R0 */

        PBXV.d = _mm256_div_pd(BGCOEF, _mm256_mul_pd(PROD.d, _mm256_mul_pd(tmp2.d, _mm256_mul_pd(tmp2.d, tmp2.d)))); /* Calculate Hernquist Profile */
        BGP.d  = _mm256_add_pd(BGP.d, _mm256_mul_pd(QI, PBXV.d));  /* Add probability for this part to total BGP */
        /* Not sure if _mm256_exp_pd even exists but we don't actually compile the AVX so this will probably need to get fixed before we can try to do that. */
        BGP.d  = _mm256_add_pd(BGP.d, _mm256_mul_pd(_mm256_mul_pd(QI, THICKCOEF), _mm256_exp_pd(_mm256_add_pd(_mm256_mul_pd(CylR, THICKLS), _mm256_mul_pd(_mm256_abs_pd(CylZ), THICKHS)))));
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

/******************************************************************************************
 ********************** Broken Power Law Profile Function *********************************
 ******************************************************************************************/
//Warning Untested
static real probabilities_avx_BPL(const AstronomyParameters* ap,
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
    const __m256d INNER    = _mm256_set1_pd(ap->innerPower); /* Exponent for the inner halo */
    const __m256d OUTER    = _mm256_set1_pd(ap->alpha_delta3); /* Change in exp value from inner to outer */
    __m256d RI, QI;
    ssp_m256 xyz0, xyz1, xyz2, tmp0, tmp1, PROD, PBXV, BGP;
    /* xyz0, 1, 2 = x, y, z */
    BGP.d = _mm256_setzero_pd();

    convolve = ap->convolve;
    nStreams = ap->number_streams;

    for (i = 0; i < convolve; i += 4)
    {
    	/* Put r_point and qw_r3_n into RI and QI respectively */
        RI = _mm256_load_pd(&r_point[i]);
        QI = _mm256_load_pd(&qw_r3_N[i]);

        /* Coordinate Transform to Galactic Center XYZ */
        xyz0.d = _mm256_sub_pd(_mm256_mul_pd(RI, COSBL), SUNR0); /* X Value */

        _mm256_store_pd(&xs[i], xyz0.d);

        xyz1.d = _mm256_mul_pd(RI, SINCOSBL); /* Y Value */
        _mm256_store_pd(&ys[i], xyz1.d);

        xyz2.d = _mm256_mul_pd(RI, SINB); /* Z Value */
        tmp0.d = _mm256_mul_pd(xyz2.d, QV_RECIP); /* Squashed Z */

        _mm256_store_pd(&zs[i], xyz2.d);

        /* Finish Convert to Galactic Center XYZ */
        xyz0.d = _mm256_mul_pd(xyz0.d, xyz0.d); /* X^2 */
        xyz1.d = _mm256_mul_pd(xyz1.d, xyz1.d); /* Y^2 */
        tmp0.d = _mm256_mul_pd(tmp0.d, tmp0.d); /* (Z/q)^2 Squishy Squishy */

        tmp1.d = _mm256_add_pd(xyz0.d, _mm256_add_pd(xyz1.d, tmp0.d)); /* Calculate R^2 */ 

        PROD.d = _mm256_sqrt_pd(tmp1.d); /* Calculate R from R^2 */

        PBXV.d = _mm256_add_pd(INNER.d, _mm256_mul_pd(OUTER.d, _mm256_cmp_pd(PROD.d, R0.d, _CMP_GE_OQ))); /* Determine exponent */
        BGP.d  = _mm256_add_pd(BGP.d, _mm256_mul_pd(QI.d, _mm256_pow_pd(_mm256_div_pd(SUNR0.d, PROD.d), PBXV.d))); /* Calculates Broken Power Law */
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

ProbabilityFunc INIT_PROBABILITIES(const AstronomyParameters * ap)
{
    assert(mwAllocA32Safe());
    initExpTable();
    if(ap->background_profile == FAST_HERNQUIST)
    {
    	return probabilities_avx_hernquist;
    }
    else if (ap->background_profile == BROKEN_POWER_LAW)
    {
    	return probabilities_avx_BPL;
    }
    return probabilities_avx_hernquist;
}


