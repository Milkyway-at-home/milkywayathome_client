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

#include "milkyway_util.h"
#include "probabilities.h"
#include "probabilities_intrin.h"


#ifndef __SSE2__
  #error SSE2 not defined
#endif

#define Vneg(x) _mm_sub_pd(_mm_setzero_pd(), (x))


static inline __m128d _mm_rcp_pd(__m128d x)
{
    __m128d xmm0 = _mm_cvtps_pd(_mm_rcp_ps( _mm_cvtpd_ps(x)));
    /* Newton Raphson */
    xmm0 = _mm_mul_pd(xmm0,_mm_sub_pd(DTWO,_mm_mul_pd(x,xmm0)));
    return _mm_mul_pd(xmm0,_mm_sub_pd(DTWO,_mm_mul_pd(x,xmm0)));
}


#define CONST_128D(var, val)                                    \
    MW_ALIGN_V(16) static const double var[2] = {(val), (val)}
#define CONST_128I(var, v1, v2, v3, v4)                                 \
    MW_ALIGN_V(16) static const int var[4] = {(v1), (v2), (v3), (v4)}


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

static inline __m128d _mm_fdiv_pd(__m128d a, __m128d b) // fasted on P4 compared to _mm_div_pd !!!
{
    __m128d r, y, c;

    y = _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(b)));
    c = _mm_sub_pd(DONE, _mm_mul_pd(b,y));
    y = _mm_add_pd(y, _mm_mul_pd(y,c));
    r = _mm_mul_pd(a,y);
    c = _mm_sub_pd(a, _mm_mul_pd(b,r));
    return _mm_add_pd(r, _mm_mul_pd(y,c));
}

static inline __m128d _mm_fsqrt_pd(__m128d y)  // accurate to 1 ulp, i.e the last bit of the double precision number
{
    // cuts some corners on the numbers range but is significantly faster, employs "faithful rounding"
    __m128d x, res;
	const __m128d C0  = _mm_set1_pd(0.75);
	const __m128d C1  = _mm_set1_pd(0.0625);
    x = _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(y))); // 22bit estimate for reciprocal square root, limits range to float range, but spares the exponent extraction
    x = _mm_mul_pd(x,_mm_sub_pd(DTHREE,_mm_mul_pd(y,_mm_mul_pd(x,x)))); // first Newton iteration (44bit accurate)
    res = _mm_mul_pd(x, y);  // do final iteration directly on sqrt(y) and not on the inverse
    return(_mm_mul_pd(res, _mm_sub_pd(C0, _mm_mul_pd(C1, _mm_mul_pd(res , x)))));
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


union d2
{
    int i[4];
    unsigned int u[4];
    long int lu[2];
    double d[2];
    __m128d m;
    __m128i mi;
};

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

static real probabilities_intrinsics(const AstronomyParameters* ap,
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

    const __m128d REF_XR = _mm_set1_pd(reff_xr_rp3);

    const __m128d COSBL    = _mm_set1_pd(lbt.lCosBCos);
    const __m128d SINB     = _mm_set1_pd(lbt.bSin);
    const __m128d SINCOSBL = _mm_set1_pd(lbt.lSinBCos);
    const __m128d SUNR0    = _mm_set1_pd(ap->sun_r0);
    const __m128d R0       = _mm_set1_pd(ap->r0);
    const __m128d QV_RECIP = _mm_set1_pd(ap->q_inv);
    __m128d RI, QI;
    ssp_m128 xyz0, xyz1, xyz2, tmp0, tmp1, PROD, PBXV, BGP;

    (void) gPrime, (void) sg_dx;

    BGP.d = _mm_setzero_pd();

    convolve = ap->convolve;
    nStreams = ap->number_streams;

    for (i = 0; i < convolve; i += 2)
    {
        RI = _mm_load_pd(&r_point[i]);
        QI = _mm_load_pd(&qw_r3_N[i]);

        xyz0.d = _mm_sub_pd(_mm_mul_pd(RI, COSBL), SUNR0);
        _mm_store_pd(&xs[i], xyz0.d);

        xyz1.d = _mm_mul_pd(RI, SINCOSBL);
        _mm_store_pd(&ys[i], xyz1.d);

        xyz2.d = _mm_mul_pd(RI, SINB);
        tmp0.d = _mm_mul_pd(xyz2.d, QV_RECIP);
        _mm_store_pd(&zs[i], xyz2.d);

        xyz0.d = _mm_mul_pd(xyz0.d, xyz0.d);
        xyz1.d = _mm_mul_pd(xyz1.d, xyz1.d);
        tmp0.d = _mm_mul_pd(tmp0.d, tmp0.d);

		//PROD.d = _mm_sqrt_pd(_mm_add_pd(xyz0.d, _mm_add_pd(xyz1.d, tmp0.d)));
		PROD.d = _mm_fsqrt_pd(_mm_add_pd(xyz0.d, _mm_add_pd(xyz1.d, tmp0.d)));
        tmp1.d = _mm_add_pd(PROD.d, R0);
		//PBXV.d = _mm_fdiv_pd(DONE, _mm_mul_pd(PROD.d, _mm_mul_pd(tmp1.d, _mm_mul_pd(tmp1.d, tmp1.d)))); FASTEST on P4 but slower on Core2 i7 etc!!!
        //PBXV.d = _mm_rcp_pd(_mm_mul_pd(PROD.d, _mm_mul_pd(tmp1.d, _mm_mul_pd(tmp1.d,tmp1.d))));
        PBXV.d = _mm_div_pd(DONE, _mm_mul_pd(PROD.d, _mm_mul_pd(tmp1.d, _mm_mul_pd(tmp1.d, tmp1.d))));
        BGP.d  = _mm_add_pd(BGP.d, _mm_mul_pd(QI, PBXV.d));

    }

    BGP.d = _mm_mul_pd(BGP.d, REF_XR);

  #if defined(__SSE3__) || defined(__SSE4_1__)
	  #pragma message("Using SSE3 Native HADD")
    BGP.d = _mm_hadd_pd(BGP.d, BGP.d);
  #else // SSE2 emu hadd_pd
    #pragma message("Using SSE2 Emulated HADD")
    BGP.d = hadd_pd_SSE2(BGP.d, BGP.d);
  #endif
    _mm_store_sd(&bg_prob, BGP.d);

    for (i = 0; i < nStreams; i++)
    {
        for (j = 0; j < convolve; j++)
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

        for (k = 0; k < convolve; k += 2)
        {
          #if SSE2_EXP_LUT
            #error "SSE2 EXP LUT is broken ATM !!!"
            _mm_store_pd(&psgf[k], _mm_mul_pd(_mm_load_pd(&qw_r3_N[k]), _mm_fexpd_lut_pd(Vneg(_mm_load_pd(&psgt[k])))));
          #elif defined(__INTEL_COMPILER) && defined(__INTEL_MATH) && defined(__SSE2__)
            #pragma message("Using Intel Math EXP")
            _mm_store_pd(&psgf[k], _mm_mul_pd(_mm_load_pd(&qw_r3_N[k]), _mm_exp_pd(Vneg(_mm_load_pd(&psgt[k])))));
          #elif defined(__SSE4_1__) || defined(__SSE3__) || defined(__SSE2__) // SSE2 mod is faster than _mm_fexp_poly_pd ! If SSE4.1 is available it gains an additinal 5% speedup!
            #pragma message("Using GMX Polynomal EXP")
            _mm_store_pd(&psgf[k], _mm_mul_pd(_mm_load_pd(&qw_r3_N[k]), gmx_mm_exp_pd(Vneg(_mm_load_pd(&psgt[k])))));
           /* http://gromacs.sourcearchive.com/documentation/4.5.3-1/gmx__sse2__double_8h-source.html */
          #else
             #pragma message("Using Polynomal EXP")
            _mm_store_pd(&psgf[k], _mm_mul_pd(_mm_load_pd(&qw_r3_N[k]), _mm_fexp_poly_pd(Vneg(_mm_load_pd(&psgt[k])))));
          #endif
		}

        streamTmps[i] = 0.0;
        #pragma ivdep
        #pragma vector always
        for (j = 0; j < convolve; j++)
        {
            streamTmps[i] += psgf[j];
        }
        streamTmps[i] *= reff_xr_rp3;
    }

    return bg_prob;
}

ProbabilityFunc INIT_PROBABILITIES(void)
{
    assert(mwAllocA16Safe());
    initExpTable();

    return probabilities_intrinsics;
}

