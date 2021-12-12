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

    typedef struct MW_ALIGN_TYPE
    {
        real_0   value;
        real_0*  gradient;
        real_0** hessian;
    } real;

    #define ZERO_REAL (real){ 0.0 , NULL , NULL } //NULL doesn't work here for some reason...

#else

    typedef MW_ALIGN_TYPE_V(sizeof(real_0)) real_0 real;
    #define ZERO_REAL 0.0

#endif       /*Define types outside of main()*/

#define realsize (int) mw_exp2_0(mw_ceil_0(mw_log2_0((real_0) sizeof(real))))

#if DOUBLEPREC
    typedef MW_ALIGN_TYPE_V(2*realsize) real double2[2];
    typedef MW_ALIGN_TYPE_V(4*realsize) real double4[4];
#else
    typedef MW_ALIGN_TYPE_V(2*realsize) real float2[2];
    typedef MW_ALIGN_TYPE_V(4*realsize) real float4[4];
#endif


#ifdef __cplusplus
extern "C" {
#endif

#if AUTODIFF /*Math functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_real_const(real_0 a)
    {
        int i,j;
        real result;
	result.value = a;
        for (i=0;i<NumberOfModelParameters;i++)
        {
            result.gradient[i] = 0.0;
        }
        for (i=0;i<NumberOfModelParameters;i++)
        {
            for (j=i;j<NumberOfModelParameters;j++)
            {
                result.hessian[i][j] = 0.0;
                result.hessian[j][i] = 0.0;
            }
        }
        return result;
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_real_var(real_0 a, int n)
    {
        real result = mw_real_const(a);
        result.gradient[n] = 1.0;
        return result;
    }
    #define INIT_REAL_CONST(r,x) {(r).value = (x);}
    #define INIT_REAL_VAR(r,x,n) {INIT_REAL_CONST((r),(x)); (r).gradient[n] = 1.0}



    /*When initializing a real parameter that we want to differentiate over,
      we must specify which rows and columns of the gradient and hessian attributes
      carry the derivatives with respect to that parameter. These are the parameters
      we currently differentiate over and their assigned column number n.

      +----PARAMETER---------------------------n----------+
           Backwards Evolution Time            0
           Time Ratio                          1
           Baryonic Plummer Radius             2
           Radius Ratio                        3
           Baryonic Mass                       4
           Mass Ratio                          5
           Orbital b coord                     6
           Orbital r coord                     7
           Orbital vx coord                    8
           Orbital vy coord                    9
           Orbital vz coord                   10
           MW Bulge Mass                      11
           MW Bulge Radius                    12
           MW Disk Mass                       13
           MW Disk Scale Length               14
           MW Disk Scale Height               15
           MW Halo Mass (or Scale Velocity)   16
           MW Halo Scale Radius               17
           MW Halo Flattening Parameter (Z)   18
           LMC Mass                           19
           LMC Scale Radius                   20

      FIXME: Most of our Halo models use more parameters that we are not listed here.
      As such, only a Logarithmic Halo and a single Miyamoto-Nagai disk can be used
      with AUTODIFF. To use other models, you will need to either reassign some numbers
      or add more parameters to the list. */


    CONST_F ALWAYS_INLINE
    static inline real_0 showRealValue(real a)
    {
        return a.value;
    }

    CONST_F ALWAYS_INLINE
    static inline void setRealValue(real* a, real_0 b)
    {
        (*a).value = b;
    }

    CONST_F ALWAYS_INLINE
    static int equalReal(const real* a, const real* b)
    {
        int i,j;
        int equalValue = (a->value == b->value);
        int equalGradient = 1;
        int equalHessian = 1;
        for (i=0;i<NumberOfModelParameters;i++)
        {
            equalGradient &= (a->gradient[i] == b->gradient[i]);
        }
        for (i=0;i<NumberOfModelParameters;i++)
        {
            for (j=i;j<NumberOfModelParameters;j++)
            {
                equalHessian &= (a->hessian[i][j] == b->hessian[i][j]);
            }
        }
        return (equalValue && equalGradient && equalHessian);
    }

    /*Defined basic derivative chain rule propagation here.*/
    CONST_F ALWAYS_INLINE
    static inline real mw_AUTODIFF(real x, real y, real_0 z, real_0 dz_dx, real_0 dz_dy, real_0 d2z_dx2, real_0 d2z_dy2, real_0 d2z_dxdy)
    {
        int i,j;
        real_0 x_grad_i, y_grad_i, x_grad_j, y_grad_j, x_hess, y_hess;

        real result = mw_real_const(z);

        for (i=0;i<NumberOfModelParameters;i++)
        {
            if(!x.gradient)  //For case when gradient is NULL (as is for ZERO_REAL)
            {
                x_grad_i = 0.0;
            }
            else
            {
                x_grad_i = x.gradient[i];
            }
            if(!y.gradient)
            {
                y_grad_i = 0.0;
            }
            else
            {
                y_grad_i = y.gradient[i];
            }
            result.gradient[i] = dz_dx*x_grad_i + dz_dy*y_grad_i;
        }

        for (i=0;i<NumberOfModelParameters;i++)
        {
            for (j=i;j<NumberOfModelParameters;j++)
            {
                if(!x.gradient)  //For case when gradient is NULL (as is for ZERO_REAL)
                {
                    x_grad_i = 0.0;
                    x_grad_j = 0.0;
                    x_hess   = 0.0;
                }
                else
                {
                    x_grad_i = x.gradient[i];
                    x_grad_j = x.gradient[j];
                    x_hess   = x.hessian[i][j];
                }
                if(!y.gradient)
                {
                    y_grad_i = 0.0;
                    y_grad_j = 0.0;
                    y_hess   = 0.0;
                }
                else
                {
                    y_grad_i = y.gradient[i];
                    y_grad_j = y.gradient[j];
                    y_hess   = y.hessian[i][j];
                }
                result.hessian[i][j] = dz_dx*x_hess + dz_dy*y_hess + d2z_dx2*x_grad_i*x_grad_j + d2z_dy2*y_grad_i*y_grad_j + d2z_dxdy*(x_grad_i*y_grad_j + y_grad_i*x_grad_j);
                result.hessian[j][i] = result.hessian[i][j];
            }
        }
        return result;
    }

    /*Basic arithmetic operators*/
    CONST_F ALWAYS_INLINE
    static inline real mw_add(real a, real b)
    {
        real_0 z        = a.value + b.value;
        real_0 dz_da    = 1.0;
        real_0 dz_db    = 1.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sub(real a, real b)
    {
        real_0 z        = a.value - b.value;
        real_0 dz_da    = 1.0;
        real_0 dz_db    = -1.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mul(real a, real b)
    {
        real_0 z        = a.value * b.value;
        real_0 dz_da    = b.value;
        real_0 dz_db    = a.value;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 1.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mul_s(real a, real_0 b)
    {
        real_0 z        = a.value * b;
        real_0 dz_da    = b;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_div(real a, real b)
    {
        real_0 z        = a.value / b.value;
        real_0 dz_da    = 1.0 / b.value;
        real_0 dz_db    = -a.value / sqr_0(b.value);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 2.0 * a.value / cube_0(b.value);
        real_0 d2z_dadb = -1.0 / sqr_0(b.value);

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_neg(real a)
    {
        real_0 z        = -a.value;
        real_0 dz_da    = -1.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Trigonometric Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_sin(real a)
    {
        real_0 z        = mw_sin_0(a.value);
        real_0 dz_da    = mw_cos_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -mw_sin_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sinpi(real a)
    {
        real_0 z        = mw_sinpi_0(a.value);
        real_0 dz_da    = M_PI * mw_cospi_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -M_PI * M_PI * mw_sinpi_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cos(real a)
    {
        real_0 z        = mw_cos_0(a.value);
        real_0 dz_da    = -mw_sin_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -mw_cos_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cospi(real a)
    {
        real_0 z        = mw_cospi_0(a.value);
        real_0 dz_da    = -M_PI * mw_sinpi_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -M_PI * M_PI * mw_cospi_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_tan(real a)
    {
        real_0 z        = mw_tan_0(a.value);
        real_0 dz_da    = 1.0 / sqr_0(mw_cos_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0 * mw_tan_0(a.value) / sqr_0(mw_cos_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_tanpi(real a)
    {
        real_0 z        = mw_tanpi_0(a.value);
        real_0 dz_da    = M_PI / sqr_0(mw_cospi_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0 * M_PI * M_PI * mw_tanpi_0(a.value) / sqr_0(mw_cospi_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sincos(real a, real* sinval, real* cosval)
    {
        *sinval = mw_sin(a);
        *cosval = mw_cos(a);
    }


    /*Inverse Trigonometric Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_asin(real a)
    {
        real_0 z        = mw_asin_0(a.value);
        real_0 dz_da    = 1.0 / mw_sqrt_0(1.0 - sqr_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = a.value / threehalves_0(1.0 - sqr_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    static inline real mw_asinpi(real a)
    {
        real_0 z        = mw_asinpi_0(a.value);
        real_0 dz_da    = 1.0 / mw_sqrt_0(1.0 - sqr_0(a.value)) / M_PI;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = a.value / threehalves_0(1.0 - sqr_0(a.value)) / M_PI;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acos(real a)
    {
        real_0 z        = mw_acos_0(a.value);
        real_0 dz_da    = -1.0 / mw_sqrt_0(1.0 - sqr_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -a.value / threehalves_0(1.0 - sqr_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acospi(real a)
    {
        real_0 z        = mw_acospi_0(a.value);
        real_0 dz_da    = -1.0 / mw_sqrt_0(1.0 - sqr_0(a.value)) / M_PI;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -a.value / threehalves_0(1.0 - sqr_0(a.value)) / M_PI;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan(real a)
    {
        real_0 z        = mw_atan_0(a.value);
        real_0 dz_da    = 1.0 / (1.0 + sqr_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -2.0 * a.value / sqr_0(1.0 + sqr_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atanpi(real a)
    {
        real_0 z        = mw_atanpi_0(a.value);
        real_0 dz_da    = 1.0 / (1.0 + sqr_0(a.value)) / M_PI;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -2.0 * a.value / sqr_0(1.0 + sqr_0(a.value)) / M_PI;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan2(real a, real b)
    {
        real_0 z        = mw_atan2_0(a.value, b.value);
        real_0 dz_da    = b.value / (sqr_0(b.value) + sqr_0(a.value));
        real_0 dz_db    = -a.value / (sqr_0(b.value) + sqr_0(a.value));
        real_0 d2z_da2  = -2.0 * a.value * b.value / sqr_0(sqr_0(b.value) + sqr_0(a.value));
        real_0 d2z_db2  = 2.0 * a.value * b.value / sqr_0(sqr_0(b.value) + sqr_0(a.value));
        real_0 d2z_dadb = (sqr_0(a.value) - sqr_0(b.value)) / sqr_0(sqr_0(b.value) + sqr_0(a.value));

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan2pi(real a, real b)
    {
        real_0 z        = mw_atan2pi_0(a.value, b.value);
        real_0 dz_da    = b.value / (sqr_0(b.value) + sqr_0(a.value)) / M_PI;
        real_0 dz_db    = -a.value / (sqr_0(b.value) + sqr_0(a.value)) / M_PI;
        real_0 d2z_da2  = -2.0 * a.value * b.value / sqr_0(sqr_0(b.value) + sqr_0(a.value)) / M_PI;
        real_0 d2z_db2  = 2.0 * a.value * b.value / sqr_0(sqr_0(b.value) + sqr_0(a.value)) / M_PI;
        real_0 d2z_dadb = (sqr_0(a.value) - sqr_0(b.value)) / sqr_0(sqr_0(b.value) + sqr_0(a.value)) / M_PI;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Exponential and Logarithmic Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_log(real a)
    {
        real_0 z        = mw_log_0(a.value);
        real_0 dz_da    = 1.0 / a.value;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_log2(real a)
    {
        real_0 z        = mw_log2_0(a.value);
        real_0 dz_da    = 1.0 / a.value / mw_log_0(2.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(a.value) / mw_log_0(2.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_logb(real a)         //FIXME: Should change 2.0 in here to FLT_RADIX later
    {
        real_0 z        = mw_log10_0(a.value);
        real_0 dz_da    = 1.0 / a.value / mw_log_0(2.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(a.value) / mw_log_0(2.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_log10(real a)
    {
        real_0 z        = mw_log10_0(a.value);
        real_0 dz_da    = 1.0 / a.value / mw_log_0(10.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(a.value) / mw_log_0(10.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_log1p(real a)
    {
        real_0 z        = mw_log1p_0(a.value);
        real_0 dz_da    = 1.0 / (1.0 + a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(1.0 + a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_exp(real a)
    {
        real_0 z        = mw_exp_0(a.value);
        real_0 dz_da    = mw_exp_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_exp_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_expm1(real a)
    {
        real_0 z        = mw_expm1_0(a.value);
        real_0 dz_da    = mw_exp_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_exp_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_exp2(real a)
    {
        real_0 z        = mw_exp2_0(a.value);
        real_0 dz_da    = mw_exp2_0(a.value) * mw_log_0(2.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_exp2_0(a.value) * mw_log_0(2.0) * mw_log_0(2.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_exp10(real a)
    {
        real_0 z        = mw_exp10_0(a.value);
        real_0 dz_da    = mw_exp10_0(a.value) * mw_log_0(10.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_exp10_0(a.value) * mw_log_0(10.0) * mw_log_0(10.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Hyperbolic Trigonometric Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_sinh(real a)
    {
        real_0 z        = mw_sinh_0(a.value);
        real_0 dz_da    = mw_cosh_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_sinh_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cosh(real a)
    {
        real_0 z        = mw_cosh_0(a.value);
        real_0 dz_da    = mw_sinh_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_cosh_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_tanh(real a)
    {
        real_0 z        = mw_tanh_0(a.value);
        real_0 dz_da    = 1.0 / sqr_0(mw_cosh_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -2.0 * mw_tanh_0(a.value) / sqr_0(mw_cosh_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Inverse Hyperbolic Trigonometric Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_asinh(real a)
    {
        real_0 z        = mw_asinh_0(a.value);
        real_0 dz_da    = 1.0 / mw_sqrt_0(sqr_0(a.value) + 1.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -a.value / threehalves_0(sqr_0(a.value) + 1.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acosh(real a)
    {
        real_0 z        = mw_acosh_0(a.value);
        real_0 dz_da    = 1.0 / mw_sqrt_0(sqr_0(a.value) - 1.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -a.value / threehalves_0(sqr_0(a.value) - 1.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atanh(real a)
    {
        real_0 z        = mw_atanh_0(a.value);
        real_0 dz_da    = 1.0 / (1.0 - sqr_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0 * a.value / sqr_0(1.0 - sqr_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Polynomial and Power Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_pow(real a, real b)
    {
        real_0 z        = mw_pow_0(a.value, b.value);
        real_0 dz_da    = b.value * mw_pow_0(a.value, b.value - 1.0);
        real_0 dz_db    = mw_pow_0(a.value, b.value) * mw_log_0(a.value);
        real_0 d2z_da2  = b.value * (b.value - 1.0) * mw_pow_0(a.value, b.value - 2.0);
        real_0 d2z_db2  = mw_pow_0(a.value, b.value) * sqr_0(mw_log_0(a.value));
        real_0 d2z_dadb = mw_pow_0(a.value, b.value - 1.0) * (1.0 + b.value*mw_log_0(a.value));

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_powr(real a, real b)
    {
        real_0 z        = mw_powr_0(a.value, b.value);
        real_0 dz_da    = b.value * mw_pow_0(a.value, b.value - 1.0);
        real_0 dz_db    = mw_pow_0(a.value, b.value) * mw_log_0(a.value);
        real_0 d2z_da2  = b.value * (b.value - 1.0) * mw_pow_0(a.value, b.value - 2.0);
        real_0 d2z_db2  = mw_pow_0(a.value, b.value) * sqr_0(mw_log_0(a.value));
        real_0 d2z_dadb = mw_pow_0(a.value, b.value - 1.0) * (1.0 + b.value*mw_log_0(a.value));

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real sqr(real a)
    {
        real_0 z        = sqr_0(a.value);
        real_0 dz_da    = 2.0 * a.value;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real cube(real a)
    {
        real_0 z        = cube_0(a.value);
        real_0 dz_da    = 3.0 * sqr_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 6.0 * a.value;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real fourth(real a)
    {
        real_0 z        = fourth_0(a.value);
        real_0 dz_da    = 4.0 * cube_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 12.0 * sqr_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real fifth(real a)
    {
        real_0 z        = fifth_0(a.value);
        real_0 dz_da    = 5.0 * fourth_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 20.0 * cube_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real sixth(real a)
    {
        real_0 z        = sixth_0(a.value);
        real_0 dz_da    = 6.0 * fifth_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 30.0 * fourth_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real inv(real a)
    {
        real_0 z        = inv_0(a.value);
        real_0 dz_da    = -1.0 / sqr_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0 / cube_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cbrt(real a)
    {
        real_0 z        = mw_cbrt_0(a.value);
        real_0 dz_da    = 1.0 / 3.0 / mw_cbrt_0(sqr_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -2.0 / 9.0 / mw_cbrt_0(fifth_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sqrt(real a)
    {
        real_0 z        = mw_sqrt_0(a.value);
        real_0 dz_da    = 1.0 / 2.0 / mw_sqrt_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / 4.0 / cube_0(mw_sqrt_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_rsqrt(real a)
    {
        real_0 z        = mw_rsqrt_0(a.value);
        real_0 dz_da    = -1.0 / 2.0 / cube_0(mw_sqrt_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 3.0 / 4.0 / fifth_0(mw_sqrt_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real minushalf(real a)
    {
        real_0 z        = minushalf_0(a.value);
        real_0 dz_da    = -1.0 / 2.0 / cube_0(mw_sqrt_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 3.0 / 4.0 / fifth_0(mw_sqrt_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real threehalves(real a)
    {
        real_0 z        = threehalves_0(a.value);
        real_0 dz_da    = 3.0 / 2.0 * mw_sqrt_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 3.0 / 4.0 / mw_sqrt_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real fivehalves(real a)
    {
        real_0 z        = fivehalves_0(a.value);
        real_0 dz_da    = 5.0 / 2.0 * threehalves_0(a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 15.0 / 4.0 * mw_sqrt_0(a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real minusthreehalves(real a)
    {
        real_0 z        = minusthreehalves_0(a.value);
        real_0 dz_da    = -3.0 / 2.0 / mw_sqrt_0(fifth_0(a.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 15.0 / 4.0 / mw_sqrt_0(sixth_0(a.value)*a.value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real minusfivehalves(real a)
    {
        real_0 z        = minusfivehalves_0(a.value);
        real_0 dz_da    = -5.0 / 2.0 / mw_sqrt_0(sixth_0(a.value)*a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 35.0 / 4.0 / mw_sqrt_0(sixth_0(a.value)*cube_0(a.value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Non-Elementary Integral Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_erf(real a)
    {
        real_0 z        = mw_erf_0(a.value);
        real_0 dz_da    = 2.0 * mw_exp_0(-sqr_0(a.value)) / mw_sqrt_0(M_PI);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -4.0 * a.value * mw_exp_0(-sqr_0(a.value)) / mw_sqrt_0(M_PI);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_erfc(real a)
    {
        real_0 z        = mw_erfc_0(a.value);
        real_0 dz_da    = -2.0 * mw_exp_0(-sqr_0(a.value)) / mw_sqrt_0(M_PI);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 4.0 * a.value * mw_exp_0(-sqr_0(a.value)) / mw_sqrt_0(M_PI);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Miscellaneous Computer Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_abs(real a)
    {
        real_0 z        = mw_abs_0(a.value);
        real_0 dz_da    = (a.value > 0) - (0 > a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fabs(real a)
    {
        real_0 z        = mw_fabs_0(a.value);
        real_0 dz_da    = (a.value > 0) - (0 > a.value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_ceil(real a)
    {
        real_0 z        = mw_ceil_0(a.value);
        real_0 dz_da    = 0.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_floor(real a)
    {
        real_0 z        = mw_floor_0(a.value);
        real_0 dz_da    = 0.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_copysign(real a, real b)
    {
        real_0 z        = mw_copysign_0(a.value, b.value);
        real_0 dz_da    = ((a.value > 0) - (0 > a.value)) * ((b.value > 0) - (0 > b.value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmax(real a, real b)
    {
        real_0 z        = mw_fmax_0(a.value, b.value);
        real_0 dz_da    = (((a.value) > (b.value)) ? 1.0 : 0.0);
        real_0 dz_db    = (((a.value) > (b.value)) ? 0.0 : 1.0);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmin(real a, real b)
    {
        real_0 z        = mw_fmin_0(a.value, b.value);
        real_0 dz_da    = (((a.value) < (b.value)) ? 1.0 : 0.0);
        real_0 dz_db    = (((a.value) < (b.value)) ? 0.0 : 1.0);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fdim(real a, real b)
    {
        real_0 z        = mw_fdim_0(a.value, b.value);
        real_0 dz_da    = (((a.value-b.value) > 0) ? 1.0 : 0.0);
        real_0 dz_db    = (((a.value-b.value) > 0) ? -1.0 : 0.0);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmod(real a, real b)
    {
        real_0 z        = mw_fmod_0(a.value, b.value);
        real_0 dz_da    = 1.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mad(real a, real b, real c)
    {
        return mw_add(mw_mul(a,b),c);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_nextafter(real a, real b)
    {
        real_0 z        = mw_nextafter_0(a.value, b.value);
        real_0 dz_da    = 1.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_round(real a)
    {
        real_0 z        = mw_round_0(a.value);
        real_0 dz_da    = 1.0;      //Done this way to preserve derivative information
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_hypot(real a, real b)
    {
        real_0 z        = mw_hypot_0(a.value, b.value);
        real_0 dz_da    = a.value / mw_hypot_0(a.value, b.value);
        real_0 dz_db    = b.value / mw_hypot_0(a.value, b.value);
        real_0 d2z_da2  = sqr_0(b.value) / cube_0(mw_hypot_0(a.value, b.value));
        real_0 d2z_db2  = sqr_0(a.value) / cube_0(mw_hypot_0(a.value, b.value));
        real_0 d2z_dadb = -a.value * b.value / cube_0(mw_hypot_0(a.value, b.value));

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real d2r(real a)
    {
        real_0 z        = d2r_0(a.value);
        real_0 dz_da    = M_PI / 180.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real r2d(real a)
    {
        real_0 z        = r2d_0(a.value);
        real_0 dz_da    = 180.0 / M_PI;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, ZERO_REAL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


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
    static inline real_0 showRealValue(real a)
    {
        return a;
    }

    CONST_F ALWAYS_INLINE
    static inline void setRealValue(real* a, real b)
    {
        (*a) = b;
    }

    CONST_F ALWAYS_INLINE
    static int equalReal(const real* a, const real* b)
    {
        return (*a == *b);
    }

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
    static inline real mw_mul_s(real a, real b)
    {
        return (a*b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_div(real a, real b)
    {
        return (a/b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_neg(real a)
    {
        return (-a);
    }

    #define mw_sin      mw_sin_0
    #define mw_sinpi    mw_sinpi_0
    #define mw_cos      mw_cos_0
    #define mw_cospi    mw_cospi_0
    #define mw_tan      mw_tan_0
    #define mw_tanpi    mw_tanpi_0
    #define mw_sincos   mw_sincos_0

    #define mw_asin     mw_asin_0
    #define mw_asinpi   mw_asinpi_0
    #define mw_acos     mw_acos_0
    #define mw_acospi   mw_acospi_0
    #define mw_atan     mw_atan_0
    #define mw_atanpi   mw_atanpi_0
    #define mw_atan2    mw_atan2_0
    #define mw_atan2pi  mw_atan2pi_0

    #define mw_log      mw_log_0
    #define mw_log2     mw_log2_0
    #define mw_logb     mw_logb_0
    #define mw_log10    mw_log10_0
    #define mw_log1p    mw_log1p_0
    #define mw_exp      mw_exp_0
    #define mw_expm1    mw_expm1_0
    #define mw_exp2     mw_exp2_0
    #define mw_exp10    mw_exp10_0

    #define mw_sinh     mw_sinh_0
    #define mw_cosh     mw_cosh_0
    #define mw_tanh     mw_tanh_0
    #define mw_asinh    mw_asinh_0
    #define mw_acosh    mw_acosh_0
    #define mw_atanh    mw_atanh_0

    #define mw_pow           mw_pow_0
    #define mw_powr          mw_powr_0
    #define sqr              sqr_0
    #define cube             cube_0
    #define fourth           fourth_0
    #define fifth            fifth_0
    #define sixth            sixth_0
    #define inv              inv_0
    #define mw_cbrt          mw_cbrt_0
    #define mw_sqrt          mw_sqrt_0
    #define mw_rsqrt         mw_rsqrt_0
    #define minushalf        minushalf_0
    #define threehalves      threehalves_0
    #define fivehalves       fivehalves_0
    #define minusthreehalves minusthreehalves_0
    #define minusfivehalves  minusfivehalves_0

    #define mw_erfc      mw_erfc_0
    #define mw_erf       mw_erf_0

    #define mw_abs      mw_abs_0
    #define mw_fabs     mw_fabs_0
    #define mw_ceil     mw_ceil_0
    #define mw_floor    mw_floor_0
    #define mw_copysign mw_copysign_0
    #define mw_fmax     mw_fmax_0
    #define mw_fmin     mw_fmin_0
    #define mw_fdim     mw_fdim_0
    #define mw_fmod     mw_fmod_0
    #define mw_mad       mw_mad_0
    #define mw_nextafter mw_nextafter_0
    #define mw_round     mw_round_0
    #define mw_hypot     mw_hypot_0
    #define d2r   d2r_0
    #define r2d   r2d_0


#endif /*Math Functions*/


#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_MATH_AUTODIFF_H_ */
