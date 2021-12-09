/*
 *  Copyright (c) 2010-2022 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2018-2022 Eric Mendelsohn
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

#if !defined(_MILKYWAY_MATH_H_INSIDE_) && !defined(MILKYWAY_MATH_COMPILATION)
  #error "Only milkyway_math.h can be included directly."
#endif

#ifndef _MILKYWAY_MATH_SUPPLEMENTAL2_H_
#define _MILKYWAY_MATH_SUPPLEMENTAL2_H_

#include "milkyway_extra.h"

#ifndef DOUBLEPREC
  #error DOUBLEPREC not defined
#endif


#ifdef __cplusplus
extern "C" {
#endif

typedef struct MW_ALIGN_TYPE_V(2*realsize)
{
    real sum;
    real correction;
} Kahan;

#define ZERO_KAHAN { ZERO_REAL, ZERO_REAL }

#define CLEAR_KAHAN(k) { (k).sum = ZERO_REAL; (k).correction = ZERO_REAL; }


#define KAHAN_ADD(k, item)                                  \
    {                                                       \
        real _y = mw_sub((item), (k).correction);           \
        real _t = mw_add((k).sum, _y);                      \
        (k).correction = mw_sub(mw_sub(_t, (k).sum), _y);   \
        (k).sum = _t;                                       \
    }

/* inOut += in with 2 Kahan sums correctly combining their error terms */
#define KAHAN_REDUCTION(inOut, in)                                                         \
    {                                                                                      \
        real correctedNextTerm, newSum;                                                    \
                                                                                           \
        correctedNextTerm = mw_add((in).sum, mw_add((in).correction, (inOut).correction)); \
        newSum = mw_add((inOut).sum, correctedNextTerm);                                   \
        (inOut).correction = mw_sub(correctedNextTerm, mw_sub(newSum, (inOut).sum));     \
        (inOut).sum = newSum;                                                              \
    }


/* other useful nonstandard constants */

/* (4 * pi) / 3 */
#define PI_4_3   4.188790204786390984616857844372670512262892532500141
#define PI_2_3   2.094395102393195492308428922186335256131446266250071
#define PI_3_2   4.712388980384689857693965074919254326295754099062659
#define M_2PI    6.2831853071795864769252867665590057683943387987502
#define SQRT_2PI 2.506628274631000502415765284811045253006986740609938


#ifndef M_PI
    # define M_E	2.7182818284590452354	/* e */
    # define M_LOG2E	1.4426950408889634074	/* log_2 e */
    # define M_LOG10E	0.43429448190325182765	/* log_10 e */
    # define M_LN2	0.69314718055994530942	/* log_e 2 */
    # define M_LN10	2.30258509299404568402	/* log_e 10 */
    # define M_PI	3.14159265358979323846	/* pi */
    # define M_PI_2	1.57079632679489661923	/* pi/2 */
    # define M_PI_4	0.78539816339744830962	/* pi/4 */
    # define M_1_PI	0.31830988618379067154	/* 1/pi */
    # define M_2_PI	0.63661977236758134308	/* 2/pi */
    # define M_2_SQRTPI	1.12837916709551257390	/* 2/sqrt(pi) */
    # define M_SQRT2	1.41421356237309504880	/* sqrt(2) */    
    # define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif


/* current hubble parameter and cosmological critical density */
#define hubble      (real_0) 73.8 / 1000.0 //km/s/kpc
#define hubble_gyr  (real_0) hubble * 3.154 *inv_0(3.086) //conversion to 1/gyr -> (km/s/kpc * 3.15576e16s/gyr * 1kpc/3.086e16km)
#define pcrit_exact (real_0) 3.0 * sqr_0(hubble_gyr) / (8.0 * M_PI)
#define pcrit       (real_0) 0.000679087369829744220469326744094105320596648627735869652//precalculated version of pcrit
#define vol_pcrit   (real_0) 0.568910904587397184785763397846734505212216314432372653620//vol_pcrit = 200.0 * pcrit * PI_4_3 

#if !defined(NAN) && defined(_MSC_VER)
static const union
{
    unsigned __int32 _mw_nanbytes[2];
    double _mw_nanval;
} _mw_nanhack = { { 0xffffffff, 0x7fffffff } };
#define NAN (_mw_nanhack._mw_nanval)
#endif /* !defined(NAN) && defined(_MSC_VER) */

#if defined(_WIN32)
  /* MSVC hacks */
  #ifndef INFINITY
    //#pragma message("FIXME: USING MAX_DOUBLE FOR INFINITY for MSVC")
    //#define INFINITY MAX_DOUBLE
    #define INFINITY HUGE_VAL
  #endif /* INFINITY */
#endif /* _WIN32 */

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_MATH_SUPPLEMENTAL2_H_ */
