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

/* TODO: Optionally use native functions */

#if !defined(_MILKYWAY_MATH_H_INSIDE_) && !defined(MILKYWAY_MATH_COMPILATION)
  #error "Only milkyway_math.h can be included directly."
#endif

#ifndef _MILKYWAY_MATH_FUNCTIONS_CL_H_
#define _MILKYWAY_MATH_FUNCTIONS_CL_H_

#ifdef __cplusplus
extern "C" {
#endif

#define mw_acos acos
#define mw_acosh acosh
#define mw_acospi acospi

#if !BROKEN_CL_MATH
  #define mw_asin asin
#else
    #warning "double asin broken on ATI"
    /* Cheap replacement using Taylor expansion that shouldn't actually be used
       asin(x) ~= x + x^3/6 + 3x^5/40 + 5x^7/112 + 35x^9/1152 + 63x^11/2816 + ...
     */
  #define mw_asin(x) ((x) + (cube(x) / 6.0) + 3.0 * (sqr(x) * cube(x) / 40.0)   \
    + (5.0 * (x) * cube(x) * cube(x) / 112.0) + (35.0 * cube(cube(x)) / 1152.0) \
    + (63.0 * sqr(x) * cube(cube(x)) / 2816.0))
#endif /* !BROKEN_CL_MATH */


#define mw_asinh asinh
#define mw_asinpi asinpi

#if !BROKEN_CL_MATH
  #define mw_atan atan
#else
  #define mw_atan(x) ((x) - cube(x) / 3.0 + sqr(x) * cube(x) / 5.0 - cube(x) * cube(x) * (x) / 7.0 + cube(cube(x)) / 9.0)
#endif /* !BROKEN_CL_MATH */

#if !BROKEN_CL_MATH
  #define mw_atan2 atan2
#else
  /* ATI currently missing double version of this.  Wikipedia gives
     this formula, but also an alternate version which avoids
     underflow but is undefined for 0.
   */
  #warning "Using replacement for missing atan2"
  #define mw_atan2(y, x) (((real) 2.0) * mw_atan((y) / (mw_hypot(x, y) + (x))))
#endif /* !BROKEN_CL_MATH */

#define mw_atanh atanh
#define mw_atanpi atanpi
#define mw_atan2pi atan2pi
#define mw_cbrt cbrt
#define mw_ceil ceil
#define mw_copysign copysign

#if !BROKEN_CL_MATH
  #define mw_cos cos
#else
    /* The ATI compiler is the buggiest garbage I've ever used. It
     * makes the compile take 60 seconds when I use the actual cos
     * function. */
  #define mw_cos(x) (1.0 - sqr(x) / 2.0 + sqr(sqr(x)) / 24.0 - sqr(cube(x)) / 720.0 + sqr(sqr(sqr(x))) / 40320.0 - (x) * cube(cube(x)) / 3628800.0)
  #warning "Using cos replacement"
#endif /* !BROKEN_CL_MATH */


#define mw_cosh cosh
#define mw_cospi cospi
#define mw_erfc erfc
#define mw_erf erf
#define mw_exp exp
#define mw_exp2 exp2
#define mw_exp10 exp10
#define mw_expm1 expm1
#define mw_fabs fabs
#define mw_fdim fdim
#define mw_floor floor

#if USE_FMA
  #define mw_fma fma
#else
  #define mw_fma(a, b, c) (((a) * (b)) + (c))
#endif
#define mw_fmax fmax
#define mw_fmin fmin
#define mw_fmod fmod
#define mw_fract fract
#define mw_frexp frexp

#if !BROKEN_CL_MATH
  #define mw_hypot hypot
#else
    #warning "Using replacement for hypot"
  #define mw_hypot(x, y) mw_sqrt(sqr(x) + sqr(y))
#endif /* !BROKEN_CL_MATH */

#define mw_ilogb ilogb
#define mw_ldexp ldexp
#define mw_lgamma lgamma
#define mw_lgamma_r lgamma_r
#define mw_log log
#define mw_log2 log2
#define mw_log10 log10
#define mw_log1p log1p
#define mw_logb logb

#if USE_MAD
  #define mw_mad mad
#else
  #define mw_mad(a, b, c) (((a) * (b)) + (c))
#endif /* USE_MAD */
#define mw_modf modf
#define mw_nan nan
#define mw_nextafter nextafter
#define mw_pow pow
#define mw_pown pown

#if !BROKEN_CL_MATH
  #define mw_powr powr
#else
  #warning "Using pow replacement for powr"
  #define mw_powr pow
#endif /* !BROKEN_CL_MATH */

#define mw_remainder remainder
#define mw_remquo remquo
#define mw_rint rint
#define mw_rootn rootn
#define mw_round round
#define mw_rsqrt rsqrt

#if !BROKEN_CL_MATH
  #define mw_sin sin
#else
  #define mw_sin(x) ((x) - cube(x) / 120.0 - (x) * cube(x) * cube(x) / 5040.0 + cube(cube(x)) / 362880.0)
#endif /* !BROKEN_CL_MATH_ */

#if !BROKEN_CL_MATH
    /* The glibc one uses 2 out arguments. We want to use it when
     * available without CL, so wrap it for consistency */
  #define mw_sincos(x, s, c) (*(s) = sincos((x), (c)))
#else
  #define mw_sincos(x, s, c) { *(s) = mw_sin(x); *(c) = mw_cos(x); }
#endif /* !BROKEN_CL_MATH */

#define mw_sinh sinh
#define mw_sinpi sinpi

#if !BROKEN_CL_MATH
  #define mw_sqrt sqrt
#else
  #if MW_RADEON_4XXX
    /* As of Stream SDK 2.2, double sqrt isn't supported on radeon 4xxx series */
    #warning "Using pow(x, 2) for sqrt replacement"
    #define mw_sqrt(x) mw_pow((x), 2.0)
  #else
    #define mw_sqrt(x) sqrt
  #endif /* Radeon 4xxx */
#endif /* !BROKEN_CL_MATH */

#define mw_tan tan
#define mw_tanh tanh
#define mw_tanpi tanpi
#define mw_tgamma tgamma
#define mw_trunc trunc

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_MATH_FUNCTIONS_CL_H_ */

