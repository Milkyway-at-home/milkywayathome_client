/*
 *  Copyright (c) 2010-2022 Rensselaer Polytechnic Institute
 *  Copyright (c) 2018-2022 Eric Mendelsohn
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

#ifndef _MILKYWAY_MATH_AUTODIFF_H_
#define _MILKYWAY_MATH_AUTODIFF_H_

#include "milkyway_extra.h"

#ifndef DOUBLEPREC
  #error DOUBLEPREC not defined
#endif

#if AUTODIFF /*Define types outside of main()*/

    #define NumberOfModelParameters 21     /*Change this number to add more space for parameters to differentiate over*/

    typedef struct MW_ALIGN_TYPE_V(sizeof(real_0)*(1 + NumberOfModelParameters + NumberOfModelParameters*NumberOfModelParameters))
    {
        real_0 value;
        real_0 gradient[NumberOfModelParameters];
        real_0 hessian[NumberOfModelParameters][NumberOfModelParameters];
    } real;

    real_0 ZERO_GRADIENT[NumberOfModelParameters] = {0.0};
    real_0 ZERO_HESSIAN[NumberOfModelParameters][NumberOfModelParameters] = {0.0};
    #define ZERO_REAL {0.0, ZERO_GRADIENT, ZERO_HESSIAN}

#else

    typedef MW_ALIGN_TYPE_V(sizeof(real_0)) real_0 real;
    #define ZERO_REAL 0.0

#endif       /*Define types outside of main()*/


#if DOUBLEPREC
    typedef MW_ALIGN_TYPE_V(2*sizeof(real)) real double2[2];
    typedef MW_ALIGN_TYPE_V(4*sizeof(real)) real double4[4];
#else
    typedef MW_ALIGN_TYPE_V(2*sizeof(real)) real float2[2];
    typedef MW_ALIGN_TYPE_V(4*sizeof(real)) real float4[4];
#endif


#ifdef __cplusplus
extern "C" {
#endif

#if AUTODIFF /*Math functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_real_const(real_0 a)
    {
        real result;
	result.value = a;
        result.gradient = ZERO_GRADIENT;
        result.hessian = ZERO_HESSIAN;
        return result;
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_real_var(real_0 a, int n)
    {
        real result;
	result.value = a;
        result.gradient = ZERO_GRADIENT;
        result.hessian = ZERO_HESSIAN;
        result.gradient[n] = 1.0;
        return result;
    }
    #define INIT_REAL_CONST(r,x) {(r).value = (x); (r).gradient = ZERO_GRADIENT; (r).hessian = ZERO_HESSIAN;}
    #define INIT_REAL_VAR(r,x,n) {INIT_REAL_CONST((r),(x)); (r).gradient[n] = 1.0}

#else
    CONST_F ALWAYS_INLINE
    static inline real mw_real_const(real a)
    {
        return a;
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_real_var(real a, int n)
    {
        return a;
    }

    #define INIT_REAL_CONST(r,x) {(r) = (x);}
    #define INIT_REAL_VAR(r,x,n) {INIT_REAL_CONST((r),(x));}

    CONST_F ALWAYS_INLINE
    static inline real mw_add(real a, real b)
    {
        return (a+b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sub(real a, real b)
    {
        return (a-b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mul(real a, real b)
    {
        return (a*b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_div(real a, real b)
    {
        return (a/b);
    }

    #define mw_sin   mw_sin_0
    #define mw_cos   mw_cos_0
    #define mw_tan   mw_tan_0
    #define mw_asin  mw_asin_0
    #define mw_acos  mw_acos_0
    #define mw_atan  mw_atan_0

    #define mw_log   mw_log_0
    #define mw_log1p mw_log1p_0
    #define mw_exp   mw_exp_0
    #define mw_expm1 mw_expm1_0
    #define mw_log10 mw_log10_0
    #define mw_sinh  mw_sinh_0
    #define mw_cosh  mw_cosh_0
    #define mw_tanh  mw_tanh_0
    #define mw_pow   mw_pow_0

    #define mw_abs   mw_fabs

    #define mw_acosh    mw_acosh_0
    #define mw_acospi   mw_acospi_0
    #define mw_asinh    mw_asinh_0
    #define mw_asinpi   mw_asinpi_0
    #define mw_atan2    mw_atan2_0
    #define mw_atanh    mw_atanh_0
    #define mw_atanpi   mw_atanpi_0
    #define mw_atan2pi  mw_atan2pi_0
    #define mw_cbrt     mw_cbrt_0
    #define mw_ceil     mw_ceil_0
    #define mw_copysign mw_copysign_0
    #define mw_cospi    mw_cospi_0
    #define mw_sinpi    mw_sinpi_0
    #define mw_tanpi    mw_tanpi_0
    #define mw_erfc     mw_erfc_0
    #define mw_erf      mw_erf_0

    #define mw_exp2     mw_exp2_0
    #define mw_exp10    mw_exp10_0

    #define mw_fabs     mw_fabs_0
    #define mw_fdim     mw_fdim_0
    #define mw_floor    mw_floor_0

    #define mw_fmod     mw_fmod_0

    #define mw_fract    mw_fract_0

    #define mw_frexp     mw_frexp_0
    #define mw_ilogb     mw_ilogb_0
    #define mw_ldexp     mw_ldexp_0
    #define mw_tgamma    mw_tgamma_0
    #define mw_tgamma_r  mw_tgamma_r_0
    #define mw_lgamma    mw_lgamma_0
    #define mw_lgamma_r  mw_lgamma_r_0
    #define mw_log2      mw_log2_0
    #define mw_logb      mw_logb_0
    #define mw_mad       mw_mad_0
    #define mw_modf      mw_modf_0
    #define mw_nan       mw_nan_0
    #define mw_nextafter mw_nextafter_0

    #define mw_powr      mw_powr_0

    #define mw_remainder mw_remainder_0
    #define mw_remquo    mw_remquo_0
    #define mw_rint      mw_rint_0
    #define mw_rootn     mw_rootn_0
    #define mw_round     mw_round_0

    #define mw_sqrt      mw_sqrt_0
    #define mw_trunc     mw_trunc_0

    #define sixth        sixth_0
    #define fifth        fifth_0
    #define fourth       fourth_0
    #define cube         cube_0
    #define sqr          sqr_0
    #define inv          inv_0

    #define fivehalves   fivehalves_0
    #define threehalves  threehalves_0

    #define minusfivehalves  minusfivehalves_0
    #define minusthreehalves minusthreehalves_0
    #define minushalf        minushalf_0

    #define mw_fma       mw_fma_0
    #define mw_hypot     mw_hypot_0

    #define mw_fmax   mw_fmax_0
    #define mw_fmin   mw_fmin_0
    #define mw_rsqrt  mw_rsqrt_0
    #define mw_sincos mw_sincos_0

    #define   ABS     ABS_0
    #define   MAX     MAX_0
    #define   MIN     MIN_0

    #define mw_cmpzero_machineeps  mw_cmpzero_machineeps_0
    #define mw_cmpnzero_machineeps mw_cmpnzero_machineeps_0

    #define mw_cmpzero_eps     mw_cmpzero_eps_0
    #define mw_cmpnzero_eps    mw_cmpnzero_eps_0

    #define mw_cmpzero_muleps  mw_cmpzero_muleps_0
    #define mw_cmpnzero_muleps mw_cmpnzero_muleps_0

    #define mw_cmpf    mw_cmpf_0
    #define mw_cmpnf   mw_cmpnf_0

    #define d2r   d2r_0
    #define r2d   r2d_0
    #define dsign dsign_0


#endif /*Math Functions*/


#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_MATH_AUTODIFF_H_ */
