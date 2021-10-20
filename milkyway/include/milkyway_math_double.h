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

typedef MW_ALIGN_TYPE_V(8) double real_0;
typedef MW_ALIGN_TYPE_V(16) double double2_0[2];
typedef MW_ALIGN_TYPE_V(32) double double4_0[4];

# define M_PI_0 ((real_0) 3.14159265358979323846) /* pi */

#define REAL_EPSILON_0 DBL_EPSILON
#define REAL_MAX_0 DBL_MAX
#define REAL_MIN_0 DBL_MIN

#if ENABLE_CRLIBM
  #define mw_sin_0   sin_rn
  #define mw_cos_0   cos_rn
  #define mw_tan_0   tan_rn
  #define mw_asin_0  asin_rn
  #define mw_acos_0  acos_rn
  #define mw_atan_0  atan_rn

  #define mw_log_0   log_rn
  #define mw_log1p_0 log1p_rn
  #define mw_exp_0   exp_rn
  #define mw_expm1_0 expm1_rn
  #define mw_log10_0 log10_rn
  #define mw_sinh_0  sinh_rn
  #define mw_cosh_0  cosh_rn
  #define mw_tanh_0  tanh_rn
  #define mw_pow_0   pow_rn
#else
  #define mw_sin_0   sin
  #define mw_cos_0   cos
  #define mw_tan_0   tan
  #define mw_asin_0  asin
  #define mw_acos_0  acos
  #define mw_atan_0  atan
  #define mw_log_0   log
  #define mw_log1p_0 log1p
  #define mw_exp_0   exp
  #define mw_expm1_0 expm1
  #define mw_log10_1 log10
  #define mw_sinh_0  sinh
  #define mw_cosh_0  cosh
  #define mw_tanh_0  tanh
  #define mw_pow_0   pow
#endif /* ENABLE_CRLIBM && DOUBLEPREC */

#define mw_abs_0   fabs

#define mw_acosh_0 acosh
#define mw_acospi_0(x) (mw_acos_0(x) / M_PI_0)
#define mw_asinh_0 asinh
#define mw_asinpi_0 (mw_asin_0(x) / M_PI_0)
#define mw_atan2_0 atan2
#define mw_atanh_0 atanh
#define mw_atanpi_0(x) (mw_atan_0(x) / M_PI_0)
#define mw_atan2pi_0(x) (mw_atan2_0(x) / M_PI_0)
#define mw_cbrt_0 cbrt
#define mw_ceil_0 ceil
#define mw_copysign_0 copysign
#define mw_cospi_0(x) mw_cos_0(M_PI_0 * (x))
#define mw_sinpi_0(x) mw_sin_0(M_PI_0 * (x))
#define mw_tanpi_0(x) (mw_tan_0(M_PI_0 * (x)))
#define mw_erfc_0 erfc
#define mw_erf_0 erf

#if HAVE_EXP2
  #define mw_exp2_0 exp2
#else
  #define mw_exp2_0(x) mw_pow_0(2.0, (x))
#endif /* HAVE_EXP2 */

#if HAVE_EXP10
  #define mw_exp10_0 exp10
#else
  #define mw_exp10_0(x) mw_powr_0((real_0) 10.0, x)
#endif /* HAVE_EXP10 */

#define mw_fabs_0 fabs
#define mw_fdim_0 fdim
#define mw_floor_0 floor

#define mw_fmod_0 fmod

/* CHECKME: mw_fract */
#define mw_fract_0(x) mw_fmin_0((x) – mw_floor_0(x), 0x1.fffffep-1f)

#define mw_frexp_0 frexp
#define mw_ilogb_0 ilogb
#define mw_ldexp_0 ldexp
#define mw_tgamma_0 tgamma
#define mw_tgamma_r_0 tgamma_r
#define mw_lgamma_0 lgamma
#define mw_lgamma_r_0 lgamma_r
#define mw_log2_0 log2
#define mw_logb_0 logb
#define mw_mad_0(a, b, c) (((a) * (b)) + (c))
#define mw_modf_0 modf
#define mw_nan_0 nan
#define mw_nextafter_0 nextafter

/* TODO: assertions that these satisfy integer y or x >= 0 */
//#define mw_pown_0(x, iy) pow(x, iy)  NOT USED ANYWHERE IN CODE!
#define mw_powr_0(x, y) pow(x, y)

#define mw_remainder_0 remainder
#define mw_remquo_0 remquo
#define mw_rint_0 rint
#define mw_rootn_0(x, y) mw_pow_0((x), 1.0 / (y))
#define mw_round_0 round

#define mw_sqrt_0 sqrt
#define mw_trunc_0 trunc


#endif /* _MILKYWAY_MATH_DOUBLE_H_ */

