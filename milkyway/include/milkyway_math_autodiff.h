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

#include <stdio.h>
#include "milkyway_extra.h"
#include "milkyway_alloc.h"

#ifndef DOUBLEPREC
  #error DOUBLEPREC not defined
#endif

#if AUTODIFF /*Define types outside of main()*/

    #define NumberOfModelParameters 20     /*Change this number to add more space for parameters to differentiate over*/
    #define HessianLength (int) (NumberOfModelParameters * (NumberOfModelParameters + 1)/2)
    #define origRealSize (int) ((NumberOfModelParameters/2.0 + 1) * (NumberOfModelParameters + 1) + 2)

    /* This code sets the size of the real struct to the next highest power of 2 using bitwise operators */
    #define adjustedRealSize ((((((origRealSize|(origRealSize>>1))|((origRealSize|(origRealSize>>1))>>2))|(((origRealSize|(origRealSize>>1))|((origRealSize|(origRealSize>>1))>>2))>>4))|((((origRealSize|(origRealSize>>1))|((origRealSize|(origRealSize>>1))>>2))|(((origRealSize|(origRealSize>>1))|((origRealSize|(origRealSize>>1))>>2))>>4))>>8))|(((((origRealSize|(origRealSize>>1))|((origRealSize|(origRealSize>>1))>>2))|(((origRealSize|(origRealSize>>1))|((origRealSize|(origRealSize>>1))>>2))>>4))|((((origRealSize|(origRealSize>>1))|((origRealSize|(origRealSize>>1))>>2))|(((origRealSize|(origRealSize>>1))|((origRealSize|(origRealSize>>1))>>2))>>4))>>8))>>16)) + 1)

    typedef struct MW_ALIGN_TYPE_V(adjustedRealSize*sizeof(real_0))
    {
        real_0 value;
        real_0 gradient[NumberOfModelParameters];
        real_0 hessian[HessianLength];
        real_0 lnfactor_gradient;      //Stores the natural log of highest gradient value to work in log-space
        real_0 lnfactor_hessian;       //Stores the natural log of highest hessian value to work in log-space
        real_0 buffer[adjustedRealSize - HessianLength - NumberOfModelParameters - 3]; //Buffer to make 'real's memory a power of 2
    } real;

    #define ZERO_REAL (real) { 0.0 , {0.0} , {0.0} , 0.0 , 0.0 , {0.0} }

#else

    typedef MW_ALIGN_TYPE_V(sizeof(real_0)) real_0 real;
    #define ZERO_REAL (0.0)

    #if DOUBLEPREC
        typedef MW_ALIGN_TYPE_V(2*sizeof(real_0)) real double2[2];
        typedef MW_ALIGN_TYPE_V(4*sizeof(real_0)) real double4[4];
    #else
        typedef MW_ALIGN_TYPE_V(2*sizeof(real_0)) real float2[2];
        typedef MW_ALIGN_TYPE_V(4*sizeof(real_0)) real float4[4];
    #endif

#endif       /*Define types outside of main()*/

    /*When initializing a real parameter that we want to differentiate over,
      we must specify which rows and columns of the gradient and hessian attributes
      carry the derivatives with respect to that parameter. These are the parameters
      we currently differentiate over and their assigned column number n.*/

    #define BACKWARDS_TIME_POS 0
    #define BARYON_RADIUS_POS 1
    #define RADIUS_RATIO_POS 2
    #define BARYON_MASS_POS 3
    #define MASS_RATIO_POS 4
    #define B_COORD_POS 5
    #define R_COORD_POS 6
    #define VX_COORD_POS 7
    #define VY_COORD_POS 8
    #define VZ_COORD_POS 9
    #define BULGE_MASS_POS 10
    #define BULGE_RADIUS_POS 11
    #define DISK_MASS_POS 12
    #define DISK_LENGTH_POS 13
    #define DISK_HEIGHT_POS 14
    #define HALO_MASS_POS 15
    #define HALO_RADIUS_POS 16
    #define HALO_ZFLATTEN_POS 17
    #define LMC_MASS_POS 18
    #define LMC_RADIUS_POS 19


#ifdef __cplusplus
extern "C" {
#endif

static inline add_logspace_0(real_0 a, real_0 b)
{
    L_g = (a<b)*b + (a>=b)*a;
    L_l = (b<a)*b + (b>=a)*a;
    return L_g + mw_log_0(1.0 + mw_exp(L_l - L_g));
}

static inline sub_logspace_0(real_0 a, real_0 b)  //Always subtracts the smaller value from the larger. Must keep note of sign separately.
{
    L_g = (a<b)*b + (a>=b)*a;
    L_l = (b<a)*b + (b>=a)*a;
    if (L_g == L_l)
    {
        return -REAL_MAX;
    }
    return L_g + mw_log_0(1.0 - mw_exp(L_l - L_g));
}

#if AUTODIFF /*Math functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_real_const(real_0 a)
    {
        real result = ZERO_REAL;
        result.value = a;
        return result;
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_real_var(real_0 a, int n)
    {
        real result = mw_real_const(a);
        if (n < NumberOfModelParameters)
        {
            result.gradient[n] = 1.0;
        }
        return result;
    }

    CONST_F ALWAYS_INLINE
    static inline real_0 showRealValue(real* a)
    {
        return a->value;
    }

    CONST_F ALWAYS_INLINE
    static inline void setRealValue(real* a, real_0 b)
    {
        (*a).value = b;
    }

    CONST_F ALWAYS_INLINE
    static inline void setRealGradient(real* a, real_0 b, int i)
    {
        (*a).gradient[i] = b;
    }

    CONST_F ALWAYS_INLINE
    static inline void setRealHessian(real* a, real_0 b, int i, int j)
    {
        int eff_i, eff_j, k;
        if(i<j)
        {
            eff_i = j;
            eff_j = i;
        }
        else
        {
            eff_i = i;
            eff_j = j;
        }
        k = (int) (eff_i*(eff_i+1)/2 + eff_j);
        (*a).hessian[k] = b;
    }

    CONST_F ALWAYS_INLINE
    static int equalReal(const real* a, const real* b)
    {
        int i,j,k;
        int equalValue = ((a->value) == (b->value));
        int equalGradient = 1;
        int equalHessian = 1;
        for (i=0;i<NumberOfModelParameters;i++)
        {
            equalGradient &= (a->gradient[i] == b->gradient[i]);
        }
        for (i=0;i<NumberOfModelParameters;i++)
        {
            for (j=0;j<i+1;j++)
            {
                k = (int) (i*(i+1)/2 + j);
                equalHessian &= (a->hessian[k] == b->hessian[k]);
            }
        }
        return (equalValue && equalGradient && equalHessian);
    }

    /*Defined basic derivative chain rule propagation here.*/
    CONST_F ALWAYS_INLINE
    static inline real mw_AUTODIFF(real* x, real* y, real_0 z, real_0 dz_dx, real_0 dz_dy, real_0 d2z_dx2, real_0 d2z_dy2, real_0 d2z_dxdy)
    {
        int i,j,k;
        real_0 x_grad_i, y_grad_i, x_grad_j, y_grad_j, x_hess, y_hess;
        int sgn_grad[NumberOfModelParameters] = {0};
        int sgn_hess[HessianLength] = {0};

        real result = mw_real_const(z);

        for (i=0;i<NumberOfModelParameters;i++)
        {
            x_grad_i = x->gradient[i];
            x_lngrad_i = x->lnfactor_gradient[i];
            if (!y)
            {
                y_grad_i = 0.0;
                y_lngrad_i = 0.0;
            }
            else
            {
                y_grad_i = y->gradient[i];
                y_lngrad_i = y->lnfactor_gradient[i];
            }

            result.gradient[i] = dz_dx*x_grad_i + dz_dy*y_grad_i;
        }

        for (i=0;i<NumberOfModelParameters;i++)
        {
            for (j=0;j<i+1;j++)
            {
                k = (int) (i*(i+1)/2 + j);
                x_grad_i = x->gradient[i];
                x_grad_j = x->gradient[j];
                x_hess   = x->hessian[k];

                if (!y)
                {
                    y_grad_i = 0.0;
                    y_grad_j = 0.0;
                    y_hess   = 0.0;
                }
                else
                {
                    y_grad_i = y->gradient[i];
                    y_grad_j = y->gradient[j];
                    y_hess   = y->hessian[k];
                }
                
                result.hessian[k] = dz_dx*x_hess + dz_dy*y_hess + d2z_dx2*x_grad_i*x_grad_j + d2z_dy2*y_grad_i*y_grad_j + d2z_dxdy*(x_grad_i*y_grad_j + y_grad_i*x_grad_j);
            }
        }
        return result;
    }

    /*Basic arithmetic operators*/
    CONST_F ALWAYS_INLINE
    static inline real mw_add(real* a, real* b)
    {
        real_0 z        = (a->value) + (b->value);
        real_0 dz_da    = 1.0;
        real_0 dz_db    = 1.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_add_s(real* a, real_0 b)
    {
        real_0 z        = (a->value) + b;
        real_0 dz_da    = 1.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sub(real* a, real* b)
    {
        real_0 z        = (a->value) - (b->value);
        real_0 dz_da    = 1.0;
        real_0 dz_db    = -1.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mul(real* a, real* b)
    {
        real_0 z        = (a->value) * (b->value);
        real_0 dz_da    = (b->value);
        real_0 dz_db    = (a->value);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 1.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mul_s(real* a, real_0 b)
    {
        real_0 z        = (a->value) * b;
        real_0 dz_da    = b;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_div(real* a, real* b)
    {
        real_0 z        = (a->value) / (b->value);
        real_0 dz_da    = 1.0 / (b->value);
        real_0 dz_db    = -(a->value) / sqr_0(b->value);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 2.0 * (a->value) / cube_0(b->value);
        real_0 d2z_dadb = -1.0 / sqr_0(b->value);

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_neg(real* a)
    {
        real_0 z        = -1.0 * (a->value);
        real_0 dz_da    = -1.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Trigonometric Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_sin(real* a)
    {
        real_0 z        = mw_sin_0(a->value);
        real_0 dz_da    = mw_cos_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -mw_sin_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sinpi(real* a)
    {
        real_0 z        = mw_sinpi_0(a->value);
        real_0 dz_da    = M_PI * mw_cospi_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -M_PI * M_PI * mw_sinpi_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cos(real* a)
    {
        real_0 z        = mw_cos_0(a->value);
        real_0 dz_da    = -mw_sin_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -mw_cos_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cospi(real* a)
    {
        real_0 z        = mw_cospi_0(a->value);
        real_0 dz_da    = -M_PI * mw_sinpi_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -M_PI * M_PI * mw_cospi_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_tan(real* a)
    {
        real_0 z        = mw_tan_0(a->value);
        real_0 dz_da    = 1.0 / sqr_0(mw_cos_0(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0 * mw_tan_0(a->value) / sqr_0(mw_cos_0(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_tanpi(real* a)
    {
        real_0 z        = mw_tanpi_0(a->value);
        real_0 dz_da    = M_PI / sqr_0(mw_cospi_0(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0 * M_PI * M_PI * mw_tanpi_0(a->value) / sqr_0(mw_cospi_0(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline void mw_sincos(real* a, real* sinval, real* cosval)
    {
        *sinval = mw_sin(a);
        *cosval = mw_cos(a);
    }


    /*Inverse Trigonometric Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_asin(real* a)
    {
        real_0 z        = mw_asin_0(a->value);
        real_0 dz_da    = 1.0 / mw_sqrt_0(1.0 - sqr_0(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = (a->value) / threehalves_0(1.0 - sqr_0(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    static inline real mw_asinpi(real* a)
    {
        real_0 z        = mw_asinpi_0(a->value);
        real_0 dz_da    = 1.0 / mw_sqrt_0(1.0 - sqr_0(a->value)) / M_PI;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = (a->value) / threehalves_0(1.0 - sqr_0(a->value)) / M_PI;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acos(real* a)
    {
        real_0 z        = mw_acos_0(a->value);
        real_0 dz_da    = -1.0 / mw_sqrt_0(1.0 - sqr_0(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -(a->value) / threehalves_0(1.0 - sqr_0(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acospi(real* a)
    {
        real_0 z        = mw_acospi_0(a->value);
        real_0 dz_da    = -1.0 / mw_sqrt_0(1.0 - sqr_0(a->value)) / M_PI;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -(a->value) / threehalves_0(1.0 - sqr_0(a->value)) / M_PI;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan(real* a)
    {
        real_0 z        = mw_atan_0(a->value);
        real_0 dz_da    = 1.0 / (1.0 + sqr_0(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -2.0 * (a->value) / sqr_0(1.0 + sqr_0(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atanpi(real* a)
    {
        real_0 z        = mw_atanpi_0(a->value);
        real_0 dz_da    = 1.0 / (1.0 + sqr_0(a->value)) / M_PI;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -2.0 * (a->value) / sqr_0(1.0 + sqr_0(a->value)) / M_PI;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan2(real* a, real* b)
    {
        real_0 z        = mw_atan2_0(a->value, b->value);
        real_0 dz_da    = (b->value) / (sqr_0(b->value) + sqr_0(a->value));
        real_0 dz_db    = -(a->value) / (sqr_0(b->value) + sqr_0(a->value));
        real_0 d2z_da2  = -2.0 * (a->value) * (b->value) / sqr_0(sqr_0(b->value) + sqr_0(a->value));
        real_0 d2z_db2  = 2.0 * (a->value) * (b->value) / sqr_0(sqr_0(b->value) + sqr_0(a->value));
        real_0 d2z_dadb = (sqr_0(a->value) - sqr_0(b->value)) / sqr_0(sqr_0(b->value) + sqr_0(a->value));

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan2pi(real* a, real* b)
    {
        real_0 z        = mw_atan2pi_0(a->value, b->value);
        real_0 dz_da    = (b->value) / (sqr_0(b->value) + sqr_0(a->value)) / M_PI;
        real_0 dz_db    = -(a->value) / (sqr_0(b->value) + sqr_0(a->value)) / M_PI;
        real_0 d2z_da2  = -2.0 * (a->value) * (b->value) / sqr_0(sqr_0(b->value) + sqr_0(a->value)) / M_PI;
        real_0 d2z_db2  = 2.0 * (a->value) * (b->value) / sqr_0(sqr_0(b->value) + sqr_0(a->value)) / M_PI;
        real_0 d2z_dadb = (sqr_0(a->value) - sqr_0(b->value)) / sqr_0(sqr_0(b->value) + sqr_0(a->value)) / M_PI;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Exponential and Logarithmic Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_log(real* a)
    {
        real_0 z        = mw_log_0(a->value);
        real_0 dz_da    = 1.0 / (a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_log2(real* a)
    {
        real_0 z        = mw_log2_0(a->value);
        real_0 dz_da    = 1.0 / (a->value) / mw_log_0(2.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(a->value) / mw_log_0(2.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_logb(real* a)
    {
        real_0 z        = mw_log10_0(a->value);
        real_0 dz_da    = 0.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_log10(real* a)
    {
        real_0 z        = mw_log10_0(a->value);
        real_0 dz_da    = 1.0 / (a->value) / mw_log_0(10.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(a->value) / mw_log_0(10.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_log1p(real* a)
    {
        real_0 z        = mw_log1p_0(a->value);
        real_0 dz_da    = 1.0 / (1.0 + (a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -1.0 / sqr_0(1.0 + (a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_exp(real* a)
    {
        real_0 z        = mw_exp_0(a->value);
        real_0 dz_da    = mw_exp_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_exp_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_expm1(real* a)
    {
        real_0 z        = mw_expm1_0(a->value);
        real_0 dz_da    = mw_exp_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_exp_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_exp2(real* a)
    {
        real_0 z        = mw_exp2_0(a->value);
        real_0 dz_da    = mw_exp2_0(a->value) * mw_log_0(2.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_exp2_0(a->value) * mw_log_0(2.0) * mw_log_0(2.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_exp10(real* a)
    {
        real_0 z        = mw_exp10_0(a->value);
        real_0 dz_da    = mw_exp10_0(a->value) * mw_log_0(10.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_exp10_0(a->value) * mw_log_0(10.0) * mw_log_0(10.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Hyperbolic Trigonometric Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_sinh(real* a)
    {
        real_0 z        = mw_sinh_0(a->value);
        real_0 dz_da    = mw_cosh_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_sinh_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cosh(real* a)
    {
        real_0 z        = mw_cosh_0(a->value);
        real_0 dz_da    = mw_sinh_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = mw_cosh_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    /*Inverse Hyperbolic Trigonometric Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_asinh(real* a)
    {
        real_0 z        = mw_asinh_0(a->value);
        real_0 dz_da    = 1.0 / mw_sqrt_0(sqr_0(a->value) + 1.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -(a->value) / threehalves_0(sqr_0(a->value) + 1.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acosh(real* a)
    {
        real_0 z        = mw_acosh_0(a->value);
        real_0 dz_da    = 1.0 / mw_sqrt_0(sqr_0(a->value) - 1.0);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -(a->value) / threehalves_0(sqr_0(a->value) - 1.0);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    /*Polynomial and Power Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_pow(real* a, real* b)
    {
        real_0 dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb;
        real_0 z        = mw_pow_0(a->value, b->value);
        if ((a->value) > 0)
        {
            dz_da    = (b->value) * mw_pow_0(a->value, (b->value) - 1.0);
            d2z_da2  = (b->value) * ((b->value) - 1.0) * mw_pow_0(a->value, (b->value) - 2.0);
            dz_db    = mw_pow_0(a->value, b->value) * mw_log_0(a->value);
            d2z_db2  = mw_pow_0(a->value, b->value) * sqr_0(mw_log_0(a->value));
            d2z_dadb = mw_pow_0(a->value, (b->value) - 1.0) * (1.0 + (b->value)*mw_log_0(a->value));
        }
        else
        {
            if ((a->value) != 0.0)
            {
                dz_da    = (b->value) * mw_pow_0(a->value, (b->value) - 1.0);
                d2z_da2  = (b->value) * ((b->value) - 1.0) * mw_pow_0(a->value, (b->value) - 2.0);
            }
            else
            {
                dz_da    = 0.0;
                d2z_da2  = 0.0;
            }
            dz_db    = 0.0;
            d2z_db2  = 0.0;
            d2z_dadb = 0.0;
        }

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_powr(real* a, real* b)
    {
        real_0 dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb;
        real_0 z        = mw_pow_0(a->value, b->value);
        if ((a->value) > 0)
        {
            dz_da    = (b->value) * mw_pow_0(a->value, (b->value) - 1.0);
            d2z_da2  = (b->value) * ((b->value) - 1.0) * mw_pow_0(a->value, (b->value) - 2.0);
            dz_db    = mw_pow_0(a->value, b->value) * mw_log_0(a->value);
            d2z_db2  = mw_pow_0(a->value, b->value) * sqr_0(mw_log_0(a->value));
            d2z_dadb = mw_pow_0(a->value, (b->value) - 1.0) * (1.0 + (b->value)*mw_log_0(a->value));
        }
        else
        {
            if ((a->value) != 0.0)
            {
                dz_da    = (b->value) * mw_pow_0(a->value, (b->value) - 1.0);
                d2z_da2  = (b->value) * ((b->value) - 1.0) * mw_pow_0(a->value, (b->value) - 2.0);
            }
            else
            {
                dz_da    = 0.0;
                d2z_da2  = 0.0;
            }
            dz_db    = 0.0;
            d2z_db2  = 0.0;
            d2z_dadb = 0.0;
        }

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real sqr(real* a)
    {
        real_0 z        = sqr_0(a->value);
        real_0 dz_da    = 2.0 * (a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real cube(real* a)
    {
        real_0 z        = cube_0(a->value);
        real_0 dz_da    = 3.0 * sqr_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 6.0 * (a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real fourth(real* a)
    {
        real_0 z        = fourth_0(a->value);
        real_0 dz_da    = 4.0 * cube_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 12.0 * sqr_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real fifth(real* a)
    {
        real_0 z        = fifth_0(a->value);
        real_0 dz_da    = 5.0 * fourth_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 20.0 * cube_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real sixth(real* a)
    {
        real_0 z        = sixth_0(a->value);
        real_0 dz_da    = 6.0 * fifth_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 30.0 * fourth_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real inv(real* a)
    {
        real_0 z        = inv_0(a->value);
        real_0 dz_da    = -1.0 / sqr_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 2.0 / cube_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cbrt(real* a)
    {
        real_0 dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb;
        real_0 z        = mw_cbrt_0(a->value);
        if (z != 0)
        {
            dz_da    = 1.0 / 3.0 / mw_cbrt_0(sqr_0(a->value));
            d2z_da2  = -2.0 / 9.0 / mw_cbrt_0(fifth_0(a->value));
        }
        else
        {
            dz_da    = 0.0;
            d2z_da2  = 0.0;
        }
        dz_db    = 0.0;
        d2z_db2  = 0.0;
        d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sqrt(real* a)
    {
        real_0 dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb;
        real_0 z        = mw_sqrt_0(a->value);
        if (z > 0)
        {
            dz_da    = 1.0 / 2.0 / mw_sqrt_0(a->value);
            d2z_da2  = -1.0 / 4.0 / cube_0(mw_sqrt_0(a->value));
        }
        else
        {
            dz_da    = 0.0;
            d2z_da2  = 0.0;
        }
        dz_db    = 0.0;
        d2z_db2  = 0.0;
        d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_rsqrt(real* a)
    {
        real_0 z        = mw_rsqrt_0(a->value);
        real_0 dz_da    = -1.0 / 2.0 / cube_0(mw_sqrt_0(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 3.0 / 4.0 / fifth_0(mw_sqrt_0(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real minushalf(real* a)
    {
        real_0 z        = minushalf_0(a->value);
        real_0 dz_da    = -1.0 / 2.0 / cube_0(mw_sqrt_0(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 3.0 / 4.0 / fifth_0(mw_sqrt_0(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real threehalves(real* a)
    {
        real_0 z        = threehalves_0(a->value);
        real_0 dz_da    = 3.0 / 2.0 * mw_sqrt_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 3.0 / 4.0 / mw_sqrt_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real fivehalves(real* a)
    {
        real_0 z        = fivehalves_0(a->value);
        real_0 dz_da    = 5.0 / 2.0 * threehalves_0(a->value);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 15.0 / 4.0 * mw_sqrt_0(a->value);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real minusthreehalves(real* a)
    {
        real_0 z        = minusthreehalves_0(a->value);
        real_0 dz_da    = -3.0 / 2.0 / mw_sqrt_0(fifth_0(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 15.0 / 4.0 / mw_sqrt_0(sixth_0(a->value)*(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real minusfivehalves(real* a)
    {
        real_0 z        = minusfivehalves_0(a->value);
        real_0 dz_da    = -5.0 / 2.0 / mw_sqrt_0(sixth_0(a->value)*(a->value));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 35.0 / 4.0 / mw_sqrt_0(sixth_0(a->value)*cube_0(a->value));
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Non-Elementary Integral Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_erf(real* a)
    {
        real_0 z        = mw_erf_0(a->value);
        real_0 dz_da    = 2.0 * mw_exp_0(-sqr_0(a->value)) / mw_sqrt_0(M_PI);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = -4.0 * (a->value) * mw_exp_0(-sqr_0(a->value)) / mw_sqrt_0(M_PI);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_erfc(real* a)
    {
        real_0 z        = mw_erfc_0(a->value);
        real_0 dz_da    = -2.0 * mw_exp_0(-sqr_0(a->value)) / mw_sqrt_0(M_PI);
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 4.0 * (a->value) * mw_exp_0(-sqr_0(a->value)) / mw_sqrt_0(M_PI);
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }


    /*Miscellaneous Computer Functions*/
    CONST_F ALWAYS_INLINE
    static inline real mw_abs(real* a)
    {
        real_0 z        = mw_abs_0(a->value);
        real_0 dz_da    = ((a->value)/REAL_EPSILON) / mw_hypot_0(1.0, ((a->value)/REAL_EPSILON));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 1.0 / cube_0(mw_hypot_0(1.0, ((a->value)/REAL_EPSILON))) / REAL_EPSILON;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fabs(real* a)
    {
        real_0 z        = mw_fabs_0(a->value);
        real_0 dz_da    = ((a->value)/1.0e-9) / mw_hypot_0(1.0, ((a->value)/1.0e-9));
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 1.0 / cube_0(mw_hypot_0(1.0, ((a->value)/1.0e-9))) / 1.0e-9;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_ceil(real* a)
    {
        real_0 z        = mw_ceil_0(a->value);
        real_0 dz_da    = 0.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_floor(real* a)
    {
        real_0 z        = mw_floor_0(a->value);
        real_0 dz_da    = 0.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_copysign(real* a, real* b)
    {
        real_0 z        = mw_copysign_0(a->value, b->value);
        real_0 dz_da    = (2.0 / (1.0 + mw_exp_0(-(a->value)/REAL_EPSILON)) - 1.0) * (2.0 / (1.0 + mw_exp_0(-(b->value)/REAL_EPSILON)) - 1.0);
        real_0 dz_db    = mw_abs_0(a->value) / (mw_cosh_0((b->value)/REAL_EPSILON) + 1.0) / REAL_EPSILON ;
        real_0 d2z_da2  = (1.0 / (mw_cosh_0((a->value)/REAL_EPSILON) + 1.0) / REAL_EPSILON) * (2.0 / (1.0 + mw_exp_0(-(b->value)/REAL_EPSILON)) - 1.0);
        real_0 d2z_db2  = -mw_abs_0(a->value) * mw_sinh_0((b->value)/REAL_EPSILON) / sqr_0(mw_cosh_0((b->value)/REAL_EPSILON) + 1.0) / REAL_EPSILON / REAL_EPSILON;
        real_0 d2z_dadb = (2.0 / (1.0 + mw_exp_0(-(a->value)/REAL_EPSILON)) - 1.0) / (mw_cosh_0((b->value)/REAL_EPSILON) + 1.0) / REAL_EPSILON;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmax(real* a, real* b)
    {
        real_0 z        = mw_fmax_0(a->value, b->value);
        real_0 dz_da    = (((a->value) > (b->value)) ? 1.0 : 0.0);
        real_0 dz_db    = (((a->value) > (b->value)) ? 0.0 : 1.0);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmin(real* a, real* b)
    {
        real_0 z        = mw_fmin_0(a->value, b->value);
        real_0 dz_da    = (((a->value) < (b->value)) ? 1.0 : 0.0);
        real_0 dz_db    = (((a->value) < (b->value)) ? 0.0 : 1.0);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fdim(real* a, real* b)
    {
        real_0 z        = mw_fdim_0(a->value, b->value);
        real_0 dz_da    = ((((a->value)-(b->value)) > 0) ? 1.0 : 0.0);
        real_0 dz_db    = ((((a->value)-(b->value)) > 0) ? -1.0 : 0.0);
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmod(real* a, real* b)
    {
        real_0 z        = mw_fmod_0(a->value, b->value);
        real_0 dz_da    = 1.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mad(real* a, real* b, real* c)
    {
        real ab = mw_mul(a,b);
        return mw_add(&ab,c);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_nextafter(real* a, real* b)
    {
        real_0 z        = mw_nextafter_0(a->value, b->value);
        real_0 dz_da    = 1.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_round(real* a)
    {
        real_0 z        = mw_round_0(a->value);
        real_0 dz_da    = 1.0;      //Done this way to preserve derivative information
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_hypot(real* a, real* b)
    {
        real_0 dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb;
        real_0 z = mw_sqrt_0(sqr_0(a->value) + sqr_0(b->value));
        if (z == 0.0)
        {
            dz_da    = 0.0;
            dz_db    = 0.0;
            d2z_da2  = 0.0;
            d2z_db2  = 0.0;
            d2z_dadb = 0.0;
        }
        else
        {
            dz_da    = (a->value) / mw_sqrt_0(sqr_0(a->value) + sqr_0(b->value));
            dz_db    = (b->value) / mw_sqrt_0(sqr_0(a->value) + sqr_0(b->value));
            d2z_da2  = sqr_0(b->value) / cube_0(mw_sqrt_0(sqr_0(a->value) + sqr_0(b->value)));
            d2z_db2  = sqr_0(a->value) / cube_0(mw_sqrt_0(sqr_0(a->value) + sqr_0(b->value)));
            d2z_dadb = -(a->value) * (b->value) / cube_0(mw_sqrt_0(sqr_0(a->value) + sqr_0(b->value)));
        }

        return mw_AUTODIFF(a, b, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real d2r(real* a)
    {
        real_0 z        = d2r_0(a->value);
        real_0 dz_da    = M_PI / 180.0;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
    }

    CONST_F ALWAYS_INLINE
    static inline real r2d(real* a)
    {
        real_0 z        = r2d_0(a->value);
        real_0 dz_da    = 180.0 / M_PI;
        real_0 dz_db    = 0.0;
        real_0 d2z_da2  = 0.0;
        real_0 d2z_db2  = 0.0;
        real_0 d2z_dadb = 0.0;

        return mw_AUTODIFF(a, NULL, z, dz_da, dz_db, d2z_da2, d2z_db2, d2z_dadb);
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

//    #define INIT_REAL_CONST(r,x) {(r) = (x);}
//    #define INIT_REAL_VAR(r,x,n) {INIT_REAL_CONST((r),(x));}

    CONST_F ALWAYS_INLINE
    static inline real showRealValue(real* a)
    {
        return *a;
    }

    CONST_F ALWAYS_INLINE
    static inline void setRealValue(real* a, real b)
    {
        (*a) = b;
    }

    CONST_F ALWAYS_INLINE
    static inline void setRealGradient(real* a, real b, int i)
    {
        (*a) = b;
    }

    CONST_F ALWAYS_INLINE
    static inline void setRealHessian(real* a, real b, int i, int j)
    {
        (*a) = b;
    }

    CONST_F ALWAYS_INLINE
    static int equalReal(const real* a, const real* b)
    {
        return ((*a) == (*b));
    }

    //Arithmetic Operators
    CONST_F ALWAYS_INLINE
    static inline real mw_add(real* a, real* b)
    {
        return ((*a) + (*b));
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_add_s(real* a, real b)
    {
        return ((*a) + b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sub(real* a, real* b)
    {
        return ((*a) - (*b));
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mul(real* a, real* b)
    {
        return ((*a) * (*b));
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mul_s(real* a, real b)
    {
        return ((*a) * (b));
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_div(real* a, real* b)
    {
        return ((*a)/(*b));
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_neg(real* a)
    {
        return -(*a);
    }

    //Trigonometric Functions
    CONST_F ALWAYS_INLINE
    static inline real mw_sin(real* a)
    {
        return mw_sin_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sinpi(real* a)
    {
        return mw_sinpi_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cos(real* a)
    {
        return mw_cos_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cospi(real* a)
    {
        return mw_cospi_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_tan(real* a)
    {
        return mw_tan_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_tanpi(real* a)
    {
        return mw_tanpi_0(*a);
    }

    #define mw_sincos   mw_sincos_0

    //Inverse Trigonometric Functions
    CONST_F ALWAYS_INLINE
    static inline real mw_asin(real* a)
    {
        return mw_asin_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_asinpi(real* a)
    {
        return mw_asinpi_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acos(real* a)
    {
        return mw_acos_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acospi(real* a)
    {
        return mw_acospi_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan(real* a)
    {
        return mw_atan_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atanpi(real* a)
    {
        return mw_atanpi_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan2(real* a, real* b)
    {
        return mw_atan2_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_atan2pi(real* a, real* b)
    {
        return mw_atan2pi_0(*a, *b);
    }

    //Exponential and Logarithmic Functions
    CONST_F ALWAYS_INLINE
    static inline real mw_log(real* a)
    {
        return mw_log_0(*a);
    }  

    CONST_F ALWAYS_INLINE
    static inline real mw_log2(real* a)
    {
        return mw_log2_0(*a);
    }  

    CONST_F ALWAYS_INLINE
    static inline real mw_logb(real* a)
    {
        return mw_logb_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_log10(real* a)
    {
        return mw_log10_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_log1p(real* a)
    {
        return mw_log1p_0(*a);
    } 

    CONST_F ALWAYS_INLINE
    static inline real mw_exp(real* a)
    {
        return mw_exp_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_expm1(real* a)
    {
        return mw_expm1_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_exp2(real* a)
    {
        return mw_exp2_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_exp10(real* a)
    {
        return mw_exp10_0(*a);
    } 

    //Hyperbolic Trig Functions
    CONST_F ALWAYS_INLINE
    static inline real mw_sinh(real* a)
    {
        return mw_sinh_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cosh(real* a)
    {
        return mw_cosh_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_asinh(real* a)
    {
        return mw_asinh_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_acosh(real* a)
    {
        return mw_acosh_0(*a);
    }

    // Polynomial and Root Functions
    CONST_F ALWAYS_INLINE
    static inline real mw_pow(real* a, real* b)
    {
        return mw_pow_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_powr(real* a, real* b)
    {
        return mw_powr_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real sqr(real* a)
    {
        return sqr_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real cube(real* a)
    {
        return cube_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real fourth(real* a)
    {
        return fourth_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real fifth(real* a)
    {
        return fifth_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real sixth(real* a)
    {
        return sixth_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real inv(real* a)
    {
        return inv_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_cbrt(real* a)
    {
        return mw_cbrt_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_sqrt(real* a)
    {
        return mw_sqrt_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_rsqrt(real* a)
    {
        return mw_rsqrt_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real minushalf(real* a)
    {
        return minushalf_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real threehalves(real* a)
    {
        return threehalves_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real fivehalves(real* a)
    {
        return fivehalves_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real minusthreehalves(real* a)
    {
        return minusthreehalves_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real minusfivehalves(real* a)
    {
        return minusfivehalves_0(*a);
    }

    //Error Functions
    CONST_F ALWAYS_INLINE
    static inline real mw_erfc(real* a)
    {
        return mw_erfc_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_erf(real* a)
    {
        return mw_erf_0(*a);
    }

    //Other Computer Functions
    CONST_F ALWAYS_INLINE
    static inline real mw_abs(real* a)
    {
        return mw_abs_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fabs(real* a)
    {
        return mw_fabs_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_ceil(real* a)
    {
        return mw_ceil_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_floor(real* a)
    {
        return mw_floor_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_copysign(real* a, real* b)
    {
        return mw_copysign_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmax(real* a, real* b)
    {
        return mw_fmax_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmin(real* a, real* b)
    {
        return mw_fmin_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fdim(real* a, real* b)
    {
        return mw_fdim_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_fmod(real* a, real* b)
    {
        return mw_fmod_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_mad(real* a, real* b, real* c)
    {
        return mw_mad_0(*a, *b, *c);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_nextafter(real* a, real* b)
    {
        return mw_nextafter_0(*a, *b);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_round(real* a)
    {
        return mw_round_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real mw_hypot(real* a, real* b)
    {
        return mw_sqrt_0(sqr_0(*a) + sqr_0(*b));
    }

    CONST_F ALWAYS_INLINE
    static inline real d2r(real* a)
    {
        return d2r_0(*a);
    }

    CONST_F ALWAYS_INLINE
    static inline real r2d(real* a)
    {
        return r2d_0(*a);
    }

#endif /*Math Functions*/


#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_MATH_AUTODIFF_H_ */
