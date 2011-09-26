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

#include "evaluation_state.h"
#include "evaluation.h"
#include "coordinates.h"
#include "r_points.h"
#include "milkyway_util.h"
#include "calculated_constants.h"
#include "probabilities.h"
#include "probabilities_intrin.h"

#include <time.h>

#ifdef __SSE2__

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
    SEPARATION_ALIGN(16) double psgt[256], psgf[256], xyzstr[256];
    SEPARATION_ALIGN(16) mwvector xyz[256];

    const __m128d REF_XR = _mm_set1_pd(reff_xr_rp3);

    //__m128d ALPHAV = _mm_set1_pd(alpha);
    const __m128d COSBL    = _mm_set1_pd(lbt.lCosBCos);
    const __m128d SINB     = _mm_set1_pd(lbt.bSin);
    const __m128d SINCOSBL = _mm_set1_pd(lbt.lSinBCos);
    const __m128d SUNR0    = _mm_set1_pd(ap->sun_r0);
    const __m128d R0       = _mm_set1_pd(ap->r0);
    const __m128d QV_RECIP = _mm_set1_pd(ap->q_inv);
    __m128d RI, QI;
    ssp_m128 xyz0, xyz1, xyz2, tmp0, tmp1, PROD, PBXV, BGP;

    BGP.d = _mm_setzero_pd();

    convolve = ap->convolve;
    nStreams = ap->number_streams;

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

#endif /* __SSE2__ */

ProbabilityFunc INIT_PROBABILITIES()
{
  #ifdef __SSE2__
    initExpTable();
    return probabilities_intrinsics;
  #else
    return NULL;
  #endif
}

