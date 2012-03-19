/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newbergb
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
 *
 */

#include "milkyway_util.h"
#include "probabilities.h"
#include "milkyway_simd.h"
#include "separation_constants.h"


#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
  #pragma GCC diagnostic ignored "-Wunknown-pragmas"
#elif defined(_MSC_VER)
#pragma warning( disable : 4068 )
#endif /* __GNUC__ */

/* This is slower with AVX for now */
#ifndef __AVX__
  #define STREAM_VEC 1
#endif

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
    double bg_prob;
    int i, j, k, convolve, nStreams;
    MW_ALIGN_V(16) double psgt[MAX_CONVOLVE], psgf[MAX_CONVOLVE];
    MW_ALIGN_V(16) double  xs[MAX_CONVOLVE], ys[MAX_CONVOLVE], zs[MAX_CONVOLVE];

  #ifdef STREAM_VEC
    __m128d xyzstr[MAX_CONVOLVE];
  #else
    double dotted0, dotted1;
    MW_ALIGN_V(16) double xyz_norm[MAX_CONVOLVE];
    MW_ALIGN_V(16) double xyzstr[MAX_CONVOLVE];
  #endif

    const __m128d REF_XR   = mw_set1_pd(reff_xr_rp3);
    const __m128d COSBL    = mw_set1_pd(lbt.lCosBCos);
    const __m128d SINB     = mw_set1_pd(lbt.bSin);
    const __m128d SINCOSBL = mw_set1_pd(lbt.lSinBCos);
    const __m128d SUNR0    = mw_set1_pd(ap->sun_r0);
    const __m128d R0       = mw_set1_pd(ap->r0);
    const __m128d QV_RECIP = mw_set1_pd(ap->q_inv);

    __m128d RI, tmp0,tmp1, PROD, PBXV, BGP;
    __m128d xyz0, xyz1, xyz2;

    (void) gPrime, (void) sg_dx;

    BGP = mw_setzero_pd();

    convolve = ap->convolve;
    nStreams = ap->number_streams;

    #pragma ivdep
    #pragma vector always
    #pragma vector aligned
    for (i = 0; i < convolve; i += 2)
    {
        RI =  mw_load_pd(&r_point[i]);

        xyz0 = mw_fms_pd(RI, COSBL, SUNR0);
        xyz1 = mw_mul_pd(RI, SINCOSBL);
        xyz2 = mw_mul_pd(RI, SINB);

        mw_store_pd(&xs[i], xyz0);
        mw_store_pd(&ys[i], xyz1);
        mw_store_pd(&zs[i], xyz2);

        tmp0 = mw_mul_pd(xyz2, QV_RECIP);

        xyz0 = mw_mul_pd(xyz0, xyz0);
        xyz1 = mw_mul_pd(xyz1, xyz1);
        tmp0 = mw_mul_pd(tmp0, tmp0);

        PROD = mw_fsqrt_pd(mw_add_pd(xyz0, mw_add_pd(xyz1, tmp0)));
        tmp1 = mw_add_pd(PROD, R0);

        PBXV = mw_div_pd(DONE, mw_mul_pd(PROD, mw_mul_pd(tmp1, mw_mul_pd(tmp1, tmp1))));

        BGP  = mw_fma_pd(mw_load_pd(&qw_r3_N[i]), PBXV, BGP);
    }

    BGP = mw_mul_pd(BGP, REF_XR);

  #if defined  (__SSE3__) || (__SSSE3__) || (__SSE41__) && (__SSE2__)
    #pragma message("Using SSE3 native HADD")
    BGP = mw_hadd_pd(BGP,  BGP);
  #else
    #pragma message("Using SSE2 emulated HADD")
    BGP = hadd_pd_SSE2(BGP,  BGP);
  #endif

    mw_store_sd(&bg_prob, BGP);

    #pragma ivdep
    #pragma vector always
    #pragma vector aligned
    for (i = 0; i < nStreams; i++)
    {
     #ifdef STREAM_VEC
        __m128d Xvec_c, Yvec_c, Zvec_c;
        __m128d Xvec_a, Yvec_a, Zvec_a;
        __m128d x0, x1, x2;
        __m128d sc_inv, dotp0, xyz_normal;

        Xvec_c =  mw_loaddup_pd(&X(sc[i].c));
        Yvec_c =  mw_loaddup_pd(&Y(sc[i].c));
        Zvec_c =  mw_loaddup_pd(&Z(sc[i].c));

        Xvec_a =  mw_loaddup_pd(&X(sc[i].a));
        Yvec_a =  mw_loaddup_pd(&Y(sc[i].a));
        Zvec_a =  mw_loaddup_pd(&Z(sc[i].a));
        sc_inv =  mw_loaddup_pd(&sc[i].sigma_sq2_inv);
     #else
        double Xvec_c, Yvec_c, Zvec_c;
        double Xvec_a, Yvec_a, Zvec_a;
        double sc_inv;

        Xvec_c = X(sc[i].c);
        Yvec_c = Y(sc[i].c);
        Zvec_c = Z(sc[i].c);

        Xvec_a = X(sc[i].a);
        Yvec_a = Y(sc[i].a);
        Zvec_a = Z(sc[i].a);
        sc_inv = sc[i].sigma_sq2_inv;
      #endif

        #pragma ivdep
        #pragma vector always
        #pragma vector aligned
        for (j = 0; j < convolve; j += 2)
        {
         #ifdef STREAM_VEC
         #pragma message("Using hand vectorized code path")
            xyzstr[0] = mw_sub_pd(mw_load_pd(&xs[j]),Xvec_c);
            xyzstr[1] = mw_sub_pd(mw_load_pd(&ys[j]),Yvec_c);
            xyzstr[2] = mw_sub_pd(mw_load_pd(&zs[j]),Zvec_c);

            x0 = mw_mul_pd(Xvec_a,xyzstr[0]);
            x1 = mw_mul_pd(Yvec_a,xyzstr[1]);
            x2 = mw_mul_pd(Zvec_a,xyzstr[2]);

            dotp0 = mw_add_pd(x0, mw_add_pd(x1 ,x2));

            xyzstr[0] = mw_sub_pd(xyzstr[0], mw_mul_pd(dotp0, Xvec_a));
            xyzstr[1] = mw_sub_pd(xyzstr[1], mw_mul_pd(dotp0, Yvec_a));
            xyzstr[2] = mw_sub_pd(xyzstr[2], mw_mul_pd(dotp0, Zvec_a));

            x0 = mw_mul_pd(xyzstr[0],xyzstr[0]);
            x1 = mw_mul_pd(xyzstr[1],xyzstr[1]);
            x2 = mw_mul_pd(xyzstr[2],xyzstr[2]);

            xyz_normal = mw_vneg_pd(mw_add_pd(x0, mw_add_pd(x1,x2) ) );

            mw_store_pd(&psgt[j], mw_mul_pd(xyz_normal, sc_inv));
         #else
            xyzstr[0] = xs[j] - Xvec_c;             xyzstr[1] = ys[j] - Yvec_c;
            xyzstr[2] = zs[j] - Zvec_c;             xyzstr[3] = xs[j+1] - Xvec_c;
            xyzstr[4] = ys[j+1] - Yvec_c;           xyzstr[5] = zs[j+1] - Zvec_c;

            dotted0 = (Xvec_a * xyzstr[0]) + (Yvec_a * xyzstr[1]) + (Zvec_a * xyzstr[2]);
            dotted1 = (Xvec_a * xyzstr[3]) + (Yvec_a * xyzstr[4]) + (Zvec_a * xyzstr[5]);

            xyzstr[0] = xyzstr[0] - dotted0 * Xvec_a;   xyzstr[1] = xyzstr[1] - dotted0 * Yvec_a;
            xyzstr[2] = xyzstr[2] - dotted0 * Zvec_a;   xyzstr[3] = xyzstr[3] - dotted1 * Xvec_a;
            xyzstr[4] = xyzstr[4] - dotted1 * Yvec_a;   xyzstr[5] = xyzstr[5] - dotted1 * Zvec_a;

            xyz_norm[0] = (xyzstr[0] * xyzstr[0]) + (xyzstr[1] * xyzstr[1]) + (xyzstr[2] * xyzstr[2]);
            xyz_norm[1] = (xyzstr[3] * xyzstr[3]) + (xyzstr[4] * xyzstr[4]) + (xyzstr[5] * xyzstr[5]);

            psgt[j] = -xyz_norm[0] * sc_inv;   // [i].sigma_sq2_inv;
            psgt[j+1] = -xyz_norm[1] * sc_inv; //sc[i].sigma_sq2_inv;
         #endif /* STREAM_VEC */
        }

        #pragma ivdep
        #pragma vector always
        #pragma vector aligned
        for (k = 0; k < convolve; k += 4)
        {
            #pragma message("Using GMX Polynomal EXP")
            mw_store_pd(&psgf[k],   mw_mul_pd(mw_load_pd(&qw_r3_N[k]),   gmx_mm_exp_pd(mw_load_pd(&psgt[k]))));
            mw_store_pd(&psgf[k+2], mw_mul_pd(mw_load_pd(&qw_r3_N[k+2]), gmx_mm_exp_pd(mw_load_pd(&psgt[k+2]))));
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
    return probabilities_intrinsics;
}

