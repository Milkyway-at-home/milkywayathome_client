/*
 *  Copyright (c) 2010, 2011 Matthew Arsenault
 *  Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute
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

#if !defined (_MILKYWAY_MATH_H_INSIDE_) && !defined(MILKYWAY_MATH_COMPILATION)
#error "Only milkyway_math.h can be included directly."
#endif

#ifndef _MILKYWAY_MATH_DOUBLE_H_
#define _MILKYWAY_MATH_DOUBLE_H_

#include "milkyway_extra.h"

#if !DOUBLEPREC
#error Double not enabled for double math
#endif

typedef MW_ALIGN_TYPE_V(8) double real;
typedef MW_ALIGN_TYPE_V(16) double double2[2];
typedef MW_ALIGN_TYPE_V(32) double double4[4];


#define REAL_EPSILON DBL_EPSILON
#define REAL_MAX DBL_MAX
#define REAL_MIN DBL_MIN

#if ENABLE_CRLIBM
  #define mw_sin   sin_rn
  #define mw_cos   cos_rn
  #define mw_tan   tan_rn
  #define mw_asin  asin_rn
  #define mw_acos  acos_rn
  #define mw_atan  atan_rn

  #define mw_log   log_rn
  #define mw_log1p log1p_rn
  #define mw_exp   exp_rn
  #define mw_expm1 expm1_rn
  #define mw_log10 log10_rn
  #define mw_sinh  sinh_rn
  #define mw_cosh  cosh_rn
  #define mw_tanh  tanh_rn
  #define mw_pow   pow_rn
#else
  #define mw_sin   sin
  #define mw_cos   cos
  #define mw_tan   tan
  #define mw_asin  asin
  #define mw_acos  acos
  #define mw_atan  atan
  #define mw_log   log
  #define mw_log1p log1p
  #define mw_exp   exp
  #define mw_expm1 expm1
  #define mw_log10 log10
  #define mw_sinh  sinh
  #define mw_cosh  cosh
  #define mw_tanh  tanh
  #define mw_pow   pow
#endif /* ENABLE_CRLIBM && DOUBLEPREC */

#define mw_abs   fabs

#define mw_acosh acosh
#define mw_acospi(x) (mw_acos(x) / M_PI))
#define mw_asinh asinh
#define mw_asinpi (mw_asin(x) / M_PI)
#define mw_atan2 atan2
#define mw_atanh atanh
#define mw_atanpi(x) (mw_atan(x) / M_PI)
#define mw_atan2pi(x) (mw_atan2(x) / M_PI)
#define mw_cbrt cbrt
#define mw_ceil ceil
#define mw_copysign copysign
#define mw_cospi(x) mw_cos(M_PI * (x))
#define mw_sinpi(x) mw_sin(M_PI * (x))
#define mw_tanpi(x) (mw_tan(M_PI * (x)))
#define mw_erfc erfc
#define mw_erf erf

#if HAVE_EXP2
  #define mw_exp2 exp2
#else
  #define mw_exp2(x) mw_pow(2.0, (x))
#endif /* HAVE_EXP2 */

#if HAVE_EXP10
  #define mw_exp10 exp10
#else
  #define mw_exp10(x) mw_powr((real) 10.0, x)
#endif /* HAVE_EXP10 */

#define mw_fabs fabs
#define mw_fdim fdim
#define mw_floor floor

#define mw_fmod fmod

/* CHECKME: mw_fract */
#define mw_fract(x) mw_fmin((x) - mw_floor(x), 0x1.fffffep-1f)

#define mw_frexp frexp
#define mw_hypot(x, y) mw_sqrt(sqr(x) + sqr(y))
#define mw_ilogb ilogb
#define mw_ldexp ldexp
#define mw_tgamma tgamma
#define mw_tgamma_r tgamma_r
#define mw_lgamma lgamma
#define mw_lgamma_r lgamma_r
#define mw_log2 log2
#define mw_logb logb
#define mw_mad(a, b, c) (((a) * (b)) + (c))
#define mw_modf modf
#define mw_nan nan
#define mw_nextafter nextafter

/* TODO: assertions that these satisfy integer y or x >= 0 */
#define mw_pown(x, iy) pow(x, iy)
#define mw_powr(x, y) pow(x, y)

#define mw_remainder remainder
#define mw_remquo remquo
#define mw_rint rint
#define mw_rootn(x, y) mw_pow((x), 1.0 / (y))
#define mw_round round

#define mw_sqrt sqrt
#define mw_tgamma tgamma
#define mw_trunc trunc


#endif /* _MILKYWAY_MATH_DOUBLE_H_ */
