/*
 *  Copyright (c) 2011 Matthew Arsenault
 *  Copyright (c) 2011 Rensselaer Polytechnic Institute
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

#ifndef _MILKYWAY_MATH_FLOAT_H_
#define _MILKYWAY_MATH_FLOAT_H_

#include "milkyway_extra.h"

#if DOUBLEPREC
#error Double enabled with float math
#endif

#if ENABLE_CRLIBM
#error Crlibm enabled with float math
#endif


typedef MW_ALIGN_TYPE_V(4) float real;
typedef MW_ALIGN_TYPE_V(8) float float2[2];
typedef MW_ALIGN_TYPE_V(16) float float4[4];


#define REAL_EPSILON FLT_EPSILON
#define REAL_MAX FLT_MAX
#define REAL_MIN FLT_MIN


#define mw_sin   sinf
#define mw_cos   cosf
#define mw_tan   tanf
#define mw_asin  asinf
#define mw_acos  acosf
#define mw_atan  atanf
#define mw_log   logf
#define mw_log1p log1pf
#define mw_exp   expf
#define mw_expm1 expm1f
#define mw_log10 log10f
#define mw_sinh  sinhf
#define mw_cosh  coshf
#define mw_tanh  tanhf
#define mw_pow   powf

#define mw_abs   fabsf

#define mw_acosh acoshf
#define mw_acospi(x) (mw_acos(x) / M_PI))
#define mw_asinh asinhf
#define mw_asinpi (mw_asin(x) / M_PI))
#define mw_atan2 atan2f
#define mw_atanh atanhf
#define mw_atanpi(x) (mw_atan(x) / M_PI)
#define mw_atan2pi(x) (mw_atan2(x) / M_PI)
#define mw_cbrt cbrtf
#define mw_ceil ceilf
#define mw_copysign copysignf
#define mw_cospi(x) mw_cos(M_PI * (x))
#define mw_sinpi(x) mw_sin(M_PI * (x))
#define mw_tanpi(x) (mw_tan(M_PI * (x)))
#define mw_erfc erfcf
#define mw_erf erff

#if HAVE_EXP2
  #define mw_exp2 exp2f
#else
  #define mw_exp2(x) mw_pow(2.0, (x))
#endif /* HAVE_EXP2 */

#if HAVE_EXP10
  #define mw_exp10 exp10f
#else
  #define mw_exp10(x) mw_powr((real) 10.0, x)
#endif /* HAVE_EXP10 */

#define mw_fabs fabsf
#define mw_fdim fdimf
#define mw_floor floorf

#define mw_fmod fmodf

/* CHECKME: mw_fract */
#define mw_fract(x) mw_fmin((x) – mw_floor(x), 0x1.fffffep-1f)

#define mw_frexp frexpf
#define mw_ilogb ilogbf
#define mw_ldexp ldexpf
#define mw_tgamma tgammaf
#define mw_tgamma_r tgamma_rf
#define mw_lgamma lgammaf
#define mw_lgamma_r lgamma_rf
#define mw_log2 log2f
#define mw_logb logbf
#define mw_mad(a, b, c) (((a) * (b)) + (c))
#define mw_modf modff
#define mw_nan nanf
#define mw_nextafter nextafterf

/* TODO: assertions that these satisfy integer y or x >= 0 */
#define mw_pown(x, iy) powf(x, iy)
#define mw_powr(x, y) powf(x, y)

#define mw_remainder remainderf
#define mw_remquo remquof
#define mw_rint rintf
#define mw_rootn(x, y) mw_pow((x), 1.0f / (y))
#define mw_round roundf

#define mw_sqrt sqrtf
#define mw_tgamma tgammaf
#define mw_trunc truncf


#endif /* _MILKYWAY_MATH_FLOAT_H_ */

