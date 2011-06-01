/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _MILKYWAY_MATH_H_
#define _MILKYWAY_MATH_H_

#define _MILKYWAY_MATH_H_INSIDE_

#if __OPENCL_VERSION__ && DOUBLEPREC
  #ifdef cl_amd_fp64
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
  #else
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
  #endif /* cl_amd_fp64 */
#endif /* DOUBLEPREC */

#define _USE_MATH_DEFINES

#include "milkyway_config.h"
#include "real.h"
#include "milkyway_vectors.h"
#include "milkyway_math_functions.h"

#ifndef __OPENCL_VERSION__
  #include <float.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

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

#define KAHAN_ADD(k, item)                              \
    {                                                   \
        real _tmp = (k).sum;                            \
        (k).sum += item;                                \
        (k).correction += (item) - ((k).sum - _tmp);    \
    }

#if 0
#define KAHAN_ADD(k, item)                          \
    {                                               \
        real _y = (item) - (k).correction;          \
        real _t = (k).sum + (_y);                   \
        (k).correction = (_t - (k).sum) - (_y);     \
        (k).sum = _t;                               \
    }
#endif

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


#if DOUBLEPREC
  #define REAL_EPSILON DBL_EPSILON
#else
  #define REAL_EPSILON FLT_EPSILON
#endif

#if !defined(NAN) && defined(_MSC_VER) && DOUBLEPREC
  /* CHECKME, also float */
    static const union
    {
        unsigned __int32 _mw_nanbytes[2];
        double _mw_nanval;
    } _mw_nanhack = { { 0xffffffff, 0x7fffffff } };
  #define NAN (_mw_nanhack._mw_nanval)
#endif /* NAN */

#if defined(_WIN32)
/* MSVC hacks */
  #ifndef INFINITY
    //#warning "FIXME: USING MAX_DOUBLE FOR INFINITY for MSVC"
    //#define INFINITY MAX_DOUBLE
    #define INFINITY HUGE_VAL
  #endif /* INFINITY */

#endif /* _WIN32 */

#ifdef __cplusplus
}
#endif

#undef _MILKYWAY_MATH_H_INSIDE_

#endif /* _MILKYWAY_MATH_H_ */

