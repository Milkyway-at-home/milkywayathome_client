/* ************************************************************************** */
/* REAL.H: include file to support compile-time specification of precision */
/* in floating-point calculations.  If the DOUBLEPREC symbol is defined to */
/* the preprocessor, calculations are done in double precision; otherwise, */
/* they may be done in single precision. */
/* */
/* Rationale: ANSI C enables programmers to write single-precision code, */
/* but does not make it easy to change the precision of code at compile */
/* time, since different functions names are used for floating and double */
/* calculations.  This package introduces the keyword "real", which may be */
/* either float or double, and defines functions which compute with */
/* real-valued numbers. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _REAL_H_
#define _REAL_H_

#include <xmmintrin.h>
#include <emmintrin.h>

#include <math.h>
#include "nbody_types.h"
#include "nbody_config.h"

#if ENABLE_CRLIBM
  #include <crlibm.h>
#endif /* ENABLE_CRLIBM */

/* Real-valued library functions.  Most of these are actually supplied
 * by the standard C libraries.
 */

#ifndef  DOUBLEPREC

#  define rsqrt  sqrtf
#  define rsqr   fsqr
#  define rsin   sinf
#  define rcos   cosf
#  define rtan   tanf
#  define rasin  asinf
#  define racos  acosf
#  define ratan  atanf
#  define ratan2 atan2f
#  define rlog   logf
#  define rlog1p log1pf
#  define rexp   expf
#  define rexpm1 expm1f
#  define rlog10 log10f
#  define rsinh  sinhf
#  define rcosh  coshf
#  define rtanh  tanhf
#  define rpow   powf
#  define rabs   fabsf
#  define rfloor floorf
#  define rceil  ceilf
#  define rmax   fmaxf
#  define rmin   fminf
#  define rfma   fmaf
#  define rdim   fdimf
#  define rrint  rintf

#else /* Use doubles */

#if ENSURE_SSE2_SQRT

/* Temporary workaround for change in inline standard conformance */
#if !defined(__APPLE_CC__) && __GNUC__ <= 4 && __GNUC_MINOR__ < 3
  #define OLD_GCC_EXTERNINLINE extern
#else
  #define OLD_GCC_EXTERNINLINE
#endif

/* Use SSE2 to do double precision square root, ensuring the
    instruction is always used.  TODO: We can SSE2 by hand just about
    everywhere for speed. */
__attribute__ ((always_inline)) OLD_GCC_EXTERNINLINE inline double rsqrt(double x)
{
    double r;
    __m128d xv = _mm_load_sd(&x);
    _mm_store_sd(&r, _mm_sqrt_sd(xv, xv));
    return r;
}

#else

#define rsqrt  sqrt

#endif /* ENSURE_SSE2_SQRT */


/* Same for crlibm or not */
#define rsqr   sqr
#define rabs   fabs
#define rfloor floor
#define rceil  ceil
#define rmax   fmax
#define rmin   fmin
#define rfma   fma
#define rdim   fdim
#define rrint  rint
#define rcbrt  cbrt      /* FIXME: Where is this in crlibm? */
#define ratan2 atan2     /* FIXME: And this? */

#if ENABLE_CRLIBM
  #define rsin   sin_rn
  #define rcos   cos_rn
  #define rtan   tan_rn
  #define rasin  asin_rn
  #define racos  acos_rn
  #define ratan  atan_rn

  #define rlog   log_rn
  #define rlog1p log1p_rn
  #define rexp   exp_rn
  #define rexpm1 expm1_rn
  #define rlog10 log10_rn
  #define rsinh  sinh_rn
  #define rcosh  cosh_rn
  #define rtanh  tanh_rn
  #define rpow   pow_rn   /* FIXME: Don't use, unproven or something */
#else
  #define rsin   sin
  #define rcos   cos
  #define rtan   tan
  #define rasin  asin
  #define racos  acos
  #define ratan  atan
  #define rlog   log
  #define rlog1p log1p
  #define rexp   exp
  #define rexpm1 expm1
  #define rlog10 log10
  #define rsinh  sinh
  #define rcosh  cosh
  #define rtanh  tanh
  #define rpow   pow
#endif /* ENABLE_CRLIBM */

#endif /* DOUBLEPREC */

#endif /* _REAL_H_ */

