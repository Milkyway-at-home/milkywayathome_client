/*
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
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

#if !defined(_MILKYWAY_MATH_H_INSIDE_) && !defined(MILKYWAY_MATH_COMPILATION)
  #error "Only milkyway_math.h can be included directly."
#endif

#ifndef _MILKYWAY_MATH_SUPPLEMENTAL_H_
#define _MILKYWAY_MATH_SUPPLEMENTAL_H_

#include "milkyway_extra.h"

#ifndef DOUBLEPREC
  #error DOUBLEPREC not defined
#endif


#ifdef __cplusplus
extern "C" {
#endif

/* simple math macros */
#define sixth(x) ((x) * (x) * (x) * (x) * (x) * (x))
#define fifth(x) ((x) * (x) * (x) * (x) * (x))
#define fourth(x) ((x) * (x) * (x) * (x))
#define cube(x) ((x) * (x) * (x))
#define sqr(x)  ((x) * (x))
#define inv(x)  ((real) 1.0 / (x))

#define fivehalves(x) ( mw_sqrt(fifth(x) ) )
#define threehalves(x) ( mw_sqrt(cube(x)  ) )

#define minusfivehalves(x) (inv(fivehalves(x)))
#define minusthreehalves(x) (inv(threehalves(x)) )
#define minushalf(x) ( inv(mw_sqrt(x)) )


/* TODO: Have fma */
#define mw_fma(a, b, c) (((a) * (b)) + (c))

/* TODO: Have hypot */
#define mw_hypot(x, y) mw_sqrt(sqr(x) + sqr(y))


#if HAVE_FMAX
  #if DOUBLEPREC
    #define mw_fmax fmax
  #else
    #define mw_fmax fmaxf
  #endif /* DOUBLEPREC */
#else
CONST_F ALWAYS_INLINE
static inline real mw_fmax(real a, real b)
{
    return ((a >= b) || isnan(b)) ?  a : b;
}
#endif /* HAVE_FMAX */


#if HAVE_FMIN
  #if DOUBLEPREC
    #define mw_fmin   fmin
  #else
    #define mw_fmin   fminf
  #endif /*  DOUBLEPREC */
#else
CONST_F ALWAYS_INLINE
static inline real mw_fmin(real a, real b)
{
    return ((a <= b) || isnan(b)) ? a : b;
}
#endif /* HAVE_FMIN */


#if HAVE_RSQRT && USE_RSQRT
  /* warning: This loses precision */
#if DOUBLEPREC
  #define mw_rsqrt rsqrt
#else
  #define mw_rsqrt rsqrtf
#endif

#else
  #if DOUBLEPREC
    #define mw_rsqrt(x) (1.0 / sqrt(x))
  #else
    #define mw_rsqrt(x) (1.0f / sqrtf(x))
  #endif /* DOUBLEPREC */
#endif /* HAVE_RSQRT && USE_RSQRT */


#if HAVE_SINCOS
  /* Assuming glibc style, e.g. sincos(x, &sinval, &cosval) */
  #if DOUBLEPREC
    #define mw_sincos sincos
  #else
    #define mw_sincos sincosf
  #endif /* DOUBLEPREC */
#else
    /* TODO: Using real sincos() would be nice in coordinate conversions
   for speed, but it's a glibc extension. It is in opencl though.  We
   could take it from glibc, but it's in assembly and would be kind of
   annoying to add, but probably worth it. */
  #define mw_sincos(x, s, c) { *(s) = mw_sin(x); *(c) = mw_cos(x); }
#endif /* HAVE_SINCOS */




/*  ABS: returns the absolute value of its argument
 *  MAX: returns the argument with the highest value
 *  MIN: returns the argument with the lowest value
 */
#define   ABS(x)       (((x) < 0) ? -(x) : (x))
#define   MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define   MIN(x,y)     (((x) < (y)) ? (x) : (y))

/* Different variants on floating point comparisons */

/* Compare to zero using machine epsilon */
#define mw_cmpzero_machineeps(x) (mw_fabs(x) < REAL_EPSILON)
#define mw_cmpnzero_machineeps(x) (mw_fabs(x) >= REAL_EPSILON)

/* Compare to zero using custom epsilon */
#define mw_cmpzero_eps(x, eps) (mw_fabs(x) < (eps))
#define mw_cmpnzero_eps(x, eps) (mw_fabs(x) >= (eps))

/* Compare to zero using custom epsilon to multiplied by the value */
#define mw_cmpzero_muleps(x, eps) (mw_fabs(x) < (mw_fabs(x) * (eps)))
#define mw_cmpnzero_muleps(x, eps) (mw_fabs(x) >= (mw_fabs(x) * (eps)))

#define mw_cmpf(a, b, eps) (mw_fabs((a) - (b)) < (mw_fmax(mw_fabs(a), mw_fabs(b)) * (eps)))
#define mw_cmpnf(a, b, eps) (mw_fabs((a) - (b)) >= (mw_fmax(mw_fabs(a), mw_fabs(b)) * (eps)))


/* degrees to radians */
#define d2r(x) ((x) * (real) M_PI / (real) 180.0)

/* radians to degrees */
#define r2d(x) ((x) * (real) 180.0 / (real) M_PI)

#define dsign(A,B) ((B) < 0.0 ? -(A) : (A))

typedef struct MW_ALIGN_TYPE_V(16)
{
    real sum;
    real correction;
} Kahan;

#define ZERO_KAHAN { 0.0, 0.0 }

#define CLEAR_KAHAN(k) { (k).sum = 0.0; (k).correction = 0.0; }


#define KAHAN_ADD(k, item)                          \
    {                                               \
        real _y = (item) - (k).correction;          \
        real _t = (k).sum + _y;                     \
        (k).correction = (_t - (k).sum) - _y;       \
        (k).sum = _t;                               \
    }

/* inOut += in with 2 Kahan sums correctly combining their error terms */
#define KAHAN_REDUCTION(inOut, in)                                             \
    {                                                                          \
        real correctedNextTerm, newSum;                                        \
                                                                               \
        correctedNextTerm = (in).sum + ((in).correction + (inOut).correction); \
        newSum = (inOut).sum + correctedNextTerm;                              \
        (inOut).correction = correctedNextTerm - (newSum - (inOut).sum);       \
        (inOut).sum = newSum;                                                  \
    }


/* other useful nonstandard constants */

/* (4 * pi) / 3 */
#define PI_4_3 ((real) 4.188790204786390984616857844372670512262892532500141)
#define PI_2_3 ((real) 2.094395102393195492308428922186335256131446266250071)
#define PI_3_2 ((real) 4.712388980384689857693965074919254326295754099062659)
#define M_2PI ((real) 6.2831853071795864769252867665590057683943387987502)
#define SQRT_2PI ((real) 2.506628274631000502415765284811045253006986740609938)


/* Taken from glibc */
#ifndef M_PI
# define M_E		2.7182818284590452354	/* e */
# define M_LOG2E	1.4426950408889634074	/* log_2 e */
# define M_LOG10E	0.43429448190325182765	/* log_10 e */
# define M_LN2		0.69314718055994530942	/* log_e 2 */
# define M_LN10		2.30258509299404568402	/* log_e 10 */
# define M_PI		3.14159265358979323846	/* pi */
# define M_PI_2		1.57079632679489661923	/* pi/2 */
# define M_PI_4		0.78539816339744830962	/* pi/4 */
# define M_1_PI		0.31830988618379067154	/* 1/pi */
# define M_2_PI		0.63661977236758134308	/* 2/pi */
# define M_2_SQRTPI	1.12837916709551257390	/* 2/sqrt(pi) */
# define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
# define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif /* M_PI */


/* current hubble parameter and cosmological critical density */
#define hubble ((real) 73.8 / 1000.0) //km/s/kpc
#define hubble_gyr ((real) hubble * 3.154 *inv(3.086) ) //conversion to 1/gyr -> (km/s/kpc * 3.15576e16s/gyr * 1kpc/3.086e16km)
#define pcrit_exact ((real) 3.0 * sqr(hubble_gyr) / (8.0 * M_PI) )
#define pcrit       0.000679087369829744220469326744094105320596648627735869652//precalculated version of pcrit
#define vol_pcrit   0.568910904587397184785763397846734505212216314432372653620//vol_pcrit = 200.0 * pcrit * PI_4_3 

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

#endif /* _MILKYWAY_MATH_SUPPLEMENTAL_H_ */
