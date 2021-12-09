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


typedef MW_ALIGN_TYPE_V(4) float real_0;

#define REAL_EPSILON FLT_EPSILON
#define REAL_MAX FLT_MAX
#define REAL_MIN FLT_MIN


#define mw_sin_0   sinf
#define mw_cos_0   cosf
#define mw_tan_0   tanf
#define mw_asin_0  asinf
#define mw_acos_0  acosf
#define mw_atan_0  atanf
#define mw_log_0   logf
#define mw_log1p_0 log1pf
#define mw_exp_0   expf
#define mw_expm1_0 expm1f
#define mw_log10_0 log10f
#define mw_sinh_0  sinhf
#define mw_cosh_0  coshf
#define mw_tanh_0  tanhf
#define mw_pow_0   powf

#define mw_abs_0   fabsf

#define mw_acosh_0 acoshf
#define mw_acospi_0(x) (mw_acos_0(x) / M_PI)
#define mw_asinh_0 asinhf
#define mw_asinpi_0(x) (mw_asin_0(x) / M_PI)
#define mw_atan2_0 atan2f
#define mw_atanh_0 atanhf
#define mw_atanpi_0(x) (mw_atan_0(x) / M_PI)
#define mw_atan2pi_0(x,y) (mw_atan2_0(x,y) / M_PI)
#define mw_cbrt_0 cbrtf
#define mw_ceil_0 ceilf
#define mw_copysign_0 copysignf
#define mw_cospi_0(x) mw_cos_0(M_PI * (x))
#define mw_sinpi_0(x) mw_sin_0(M_PI * (x))
#define mw_tanpi_0(x) (mw_tan_0(M_PI * (x)))
#define mw_erfc_0 erfcf
#define mw_erf_0 erff

#if HAVE_EXP2
  #define mw_exp2_0 exp2f
#else
  #define mw_exp2_0(x) mw_pow_0(2.0, (x))
#endif /* HAVE_EXP2 */

#if HAVE_EXP10
  #define mw_exp10_0 exp10f
#else
  #define mw_exp10_0(x) mw_powr_0((real_0) 10.0, x)
#endif /* HAVE_EXP10 */

#define mw_fabs_0 fabsf
#define mw_fdim_0 fdimf
#define mw_floor_0 floorf

#define mw_fmod_0 fmodf

/* CHECKME: mw_fract */
//#define mw_fract_0(x) mw_fmin_0((x) – mw_floor_0(x), 0x1.fffffep-1f)  NOT USED ANYWHERE IN CODE!

//#define mw_frexp_0 frexpf                                       NOT USED ANYWHERE IN CODE!
//#define mw_ilogb_0 ilogbf                                       NOT USED ANYWHERE IN CODE!
#define mw_ldexp_0 ldexpf
#define mw_tgamma_0 tgammaf
#define mw_tgamma_r_0 tgamma_rf
#define mw_lgamma_0 lgammaf
#define mw_lgamma_r_0 lgamma_rf
#define mw_log2_0 log2f
#define mw_logb_0 logbf
#define mw_mad_0(a, b, c) (((a) * (b)) + (c))
//#define mw_modf_0 modff                                         NOT USED ANYWHERE IN CODE!
//#define mw_nan_0 nanf                                           NOT USED ANYWHERE IN CODE!
#define mw_nextafter_0 nextafterf

/* TODO: assertions that these satisfy integer y or x >= 0 */
//#define mw_pown_0(x, iy) powf(x, iy)                            NOT USED ANYWHERE IN CODE!
#define mw_powr_0(x, y) powf(x, y)

//#define mw_remainder_0 remainderf                               NOT USED ANYWHERE IN CODE!
//#define mw_remquo_0 remquof                                     NOT USED ANYWHERE IN CODE!
//#define mw_rint_0 rintf                                         NOT USED ANYWHERE IN CODE!
//#define mw_rootn_0(x, y) mw_pow_0((x), 1.0f / (y))              NOT USED ANYWHERE IN CODE!
#define mw_round_0 roundf

#define mw_sqrt_0 sqrtf
//#define mw_trunc_0 truncf                                       NOT USED ANYWHERE IN CODE!


#endif /* _MILKYWAY_MATH_FLOAT_H_ */

