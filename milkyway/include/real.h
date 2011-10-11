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

#if !defined(_MILKYWAY_MATH_H_INSIDE_) && !defined(MILKYWAY_MATH_COMPILATION)
  #error "Only milkyway_math.h can be included directly."
#endif


#ifndef _REAL_H_
#define _REAL_H_


/* We need to pull in the __GLIBC__ macros for our lazy check of
 * broken tgmath, and this is the recommended header to use for that
 * purpose
 */
#include <limits.h>

#ifndef DOUBLEPREC
  #error DOUBLEPREC not defined
#endif


#if DOUBLEPREC
typedef double real MW_ALIGN(8);
typedef double double2[2] MW_ALIGN(16);
typedef double double4[4] MW_ALIGN(32);

#else
typedef float real;
typedef float float2[2] MW_ALIGN(8);
typedef float float4[4] MW_ALIGN(16);
#endif /* DOUBLEPREC */


#if DOUBLEPREC
  #define REAL_EPSILON DBL_EPSILON
  #define REAL_MAX DBL_MAX
  #define REAL_MIN DBL_MIN
#else
  #define REAL_EPSILON FLT_EPSILON
  #define REAL_MAX FLT_MAX
  #define REAL_MIN FLT_MIN
#endif


/* FIXME: This happens to work with MSVC with double since the
   functions.  MSVC is stupid and doesn't support C99 or anything,
   so we'll have to do work to make float work to not have horrible
   casting everywhere
*/

#if (__GLIBC__ <= 2 && __GLIBC_MINOR__ < 7 && defined(__GLIBC__)) || defined(__FreeBSD__)
  #define HAVE_BROKEN_TGMATH 1
  #if !DOUBLEPREC
    #warning "This tgmath.h doesn't work, so float doesn't really work"
  #endif
#endif /* __GLIBC__ <= 2 && __GLIBC_MINOR__ < 7 */

#ifndef __cplusplus
  #if HAVE_TGMATH_H && !HAVE_BROKEN_TGMATH
    #include <tgmath.h>
  #else
    #include <math.h>
  #endif /* HAVE_TGMATH_H && !HAVE_BROKEN_TGMATH */
#else
  #include <cmath>
#endif /* __cplusplus */

/* crlibm is a math.h supplement.*/
#if ENABLE_CRLIBM
  #include <crlibm.h>
#endif /* ENABLE_CRLIBM */

#endif /* _REAL_H_ */

