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
#define sixth_0(x) ((x) * (x) * (x) * (x) * (x) * (x))
#define fifth_0(x) ((x) * (x) * (x) * (x) * (x))
#define fourth_0(x) ((x) * (x) * (x) * (x))
#define cube_0(x) ((x) * (x) * (x))
#define sqr_0(x)  ((x) * (x))
#define inv_0(x)  ((real_0) 1.0 / (x))

#define fivehalves_0(x) ( mw_sqrt_0(fifth_0(x) ) )
#define threehalves_0(x) ( mw_sqrt_0(cube_0(x)  ) )

#define minusfivehalves_0(x) (inv_0(fivehalves_0(x)))
#define minusthreehalves_0(x) (inv_0(threehalves_0(x)) )
#define minushalf_0(x) ( inv_0(mw_sqrt_0(x)) )


/* TODO: Have fma */
#define mw_fma_0(a, b, c) (((a) * (b)) + (c))

/* TODO: Have hypot */
#define mw_hypot_0(x, y) mw_sqrt_0(sqr_0(x) + sqr_0(y))


#if HAVE_FMAX
  #if DOUBLEPREC
    #define mw_fmax_0 fmax
  #else
    #define mw_fmax_0 fmaxf
  #endif /* DOUBLEPREC */
#else
CONST_F ALWAYS_INLINE
static inline real_0 mw_fmax_0(real_0 a, real_0 b)
{
    return ((a >= b) || isnan(b)) ?  a : b;
}
#endif /* HAVE_FMAX */


#if HAVE_FMIN
  #if DOUBLEPREC
    #define mw_fmin_0   fmin
  #else
    #define mw_fmin_0   fminf
  #endif /*  DOUBLEPREC */
#else
CONST_F ALWAYS_INLINE
static inline real_0 mw_fmin_0(real_0 a, real_0 b)
{
    return ((a <= b) || isnan(b)) ? a : b;
}
#endif /* HAVE_FMIN */


#if HAVE_RSQRT && USE_RSQRT
  /* warning: This loses precision */
#if DOUBLEPREC
  #define mw_rsqrt_0 rsqrt
#else
  #define mw_rsqrt_0 rsqrtf
#endif

#else
  #if DOUBLEPREC
    #define mw_rsqrt_0(x) (1.0 / sqrt(x))
  #else
    #define mw_rsqrt_0(x) (1.0f / sqrt(x))
  #endif /* DOUBLEPREC */
#endif /* HAVE_RSQRT && USE_RSQRT */


#if HAVE_SINCOS
  /* Assuming glibc style, e.g. sincos(x, &sinval, &cosval) */
  #if DOUBLEPREC
    #define mw_sincos_0 sincos
  #else
    #define mw_sincos_0 sincosf
  #endif /* DOUBLEPREC */
#else
    /* TODO: Using real sincos() would be nice in coordinate conversions
   for speed, but it's a glibc extension. It is in opencl though.  We
   could take it from glibc, but it's in assembly and would be kind of
   annoying to add, but probably worth it. */
  #define mw_sincos_0(x, s, c) { *(s) = mw_sin_0(x); *(c) = mw_cos_0(x); }
#endif /* HAVE_SINCOS */


/*  ABS: returns the absolute value of its argument
 *  MAX: returns the argument with the highest value
 *  MIN: returns the argument with the lowest value
 */
#define   ABS_0(x)       (((x) < 0) ? -(x) : (x))
#define   MAX_0(x,y)     (((x) > (y)) ? (x) : (y))
#define   MIN_0(x,y)     (((x) < (y)) ? (x) : (y))

/* Different variants on floating point comparisons */

/* Compare to zero using machine epsilon */
#define mw_cmpzero_machineeps_0(x) (mw_fabs_0(x) < REAL_EPSILON_0)
#define mw_cmpnzero_machineeps_0(x) (mw_fabs_0(x) >= REAL_EPSILON_0)

/* Compare to zero using custom epsilon */
#define mw_cmpzero_eps_0(x, eps) (mw_fabs_0(x) < (eps))
#define mw_cmpnzero_eps_0(x, eps) (mw_fabs_0(x) >= (eps))

/* Compare to zero using custom epsilon to multiplied by the value */
#define mw_cmpzero_muleps_0(x, eps) (mw_fabs_0(x) < (mw_fabs_0(x) * (eps)))
#define mw_cmpnzero_muleps_0(x, eps) (mw_fabs_0(x) >= (mw_fabs_0(x) * (eps)))

#define mw_cmpf_0(a, b, eps) (mw_fabs_0((a) - (b)) < (mw_fmax_0(mw_fabs_0(a), mw_fabs_0(b)) * (eps)))
#define mw_cmpnf_0(a, b, eps) (mw_fabs_0((a) - (b)) >= (mw_fmax_0(mw_fabs_0(a), mw_fabs_0(b)) * (eps)))


/* degrees to radians */
#define d2r_0(x) ((x) * (real_0) M_PI / (real_0) 180.0)

/* radians to degrees */
#define r2d_0(x) ((x) * (real_0) 180.0 / (real_0) M_PI)

#define dsign_0(A,B) ((B) < 0.0 ? -(A) : (A))

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_MATH_SUPPLEMENTAL_H_ */

