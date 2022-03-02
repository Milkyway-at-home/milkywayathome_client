/* This program tests the AUTODIFF derivative propagation */

#include "milkyway_util.h"
#include "nbody_bessel.h"
#include <stdio.h> 
#include <stdlib.h>
#include <time.h>

#define derv_thresh (0.05)
#define small_derv (0.001)

static inline real_0 add_func(real_0 a, real_0 b)
{
    return a+b;
}

static inline real_0 sub_func(real_0 a, real_0 b)
{
    return a-b;
}

static inline real_0 mul_func(real_0 a, real_0 b)
{
    return a*b;
}

static inline real_0 div_func(real_0 a, real_0 b)
{
    return a/b;
}

static inline real_0 neg_func(real_0 a)
{
    return -a;
}

static inline real_0 sinpi_func(real_0 a)
{
    return mw_sinpi_0(a);
}

static inline real_0 cospi_func(real_0 a)
{
    return mw_cospi_0(a);
}

static inline real_0 tanpi_func(real_0 a)
{
    return mw_tanpi_0(a);
}

static inline real_0 asinpi_func(real_0 a)
{
    return mw_asinpi_0(a);
}

static inline real_0 acospi_func(real_0 a)
{
    return mw_acospi_0(a);
}

static inline real_0 atanpi_func(real_0 a)
{
    return mw_atanpi_0(a);
}

static inline real_0 atan2pi_func(real_0 a, real_0 b)
{
    return mw_atan2pi_0(a, b);
}

static inline real_0 powr_func(real_0 a, real_0 b)
{
    return mw_powr_0(a, b);
}

static inline real_0 sqr_func(real_0 a)
{
    return sqr_0(a);
}

static inline real_0 cube_func(real_0 a)
{
    return cube_0(a);
}

static inline real_0 fourth_func(real_0 a)
{
    return fourth_0(a);
}

static inline real_0 fifth_func(real_0 a)
{
    return fifth_0(a);
}

static inline real_0 sixth_func(real_0 a)
{
    return sixth_0(a);
}

static inline real_0 inv_func(real_0 a)
{
    return inv_0(a);
}

static inline real_0 rsqrt_func(real_0 a)
{
    return mw_rsqrt_0(a);
}

static inline real_0 minushalf_func(real_0 a)
{
    return minushalf_0(a);
}

static inline real_0 threehalves_func(real_0 a)
{
    return threehalves_0(a);
}

static inline real_0 fivehalves_func(real_0 a)
{
    return fivehalves_0(a);
}

static inline real_0 minusthreehalves_func(real_0 a)
{
    return minusthreehalves_0(a);
}

static inline real_0 minusfivehalves_func(real_0 a)
{
    return minusfivehalves_0(a);
}

static inline real_0 hypot_func(real_0 a, real_0 b)
{
    return mw_sqrt_0(sqr_0(a) + sqr_0(b));
}

static inline real_0 d2r_func(real_0 a)
{
    return d2r_0(a);
}

static inline real_0 r2d_func(real_0 a)
{
    return r2d_0(a);
}

static inline int testSingleVar(real (*func)(real*), real_0 (*func0)(real_0), real_0 min_a, real_0 max_a)
{
    int failed = 0;
    real_0 a_val;
    real a, result;
    real_0 h = 0.0001;
    real_0 deriv, check;
    real_0 p1, p2, p3, p4, p5, denom;

    srand(time(0));

    for (int i = 0; i < 100000; i++)
    {
        a_val = ((real_0)rand()/(real_0)RAND_MAX) * (max_a - min_a) + min_a;
        a = mw_real_var(a_val, 0);
        result = (*func)(&a);

        p1 =   1.0 * (*func0)( (a_val - 2.0 * h) );
        p2 = - 8.0 * (*func0)( (a_val - h) );
        p3 = - 1.0 * (*func0)( (a_val + 2.0 * h) );
        p4 =   8.0 * (*func0)( (a_val + h) );
        denom = inv_0( 12.0 * h);
        deriv = (p1 + p2 + p3 + p4) * denom;

        if(mw_abs_0(result.gradient[0]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.gradient[0] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD GRADIENT! At (a=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, result.gradient[0], deriv, check);
        }

        p1 = - 1.0 * (*func0)( (a_val + 2.0 * h) );
        p2 =  16.0 * (*func0)( (a_val + h) );
        p3 = -30.0 * (*func0)( (a_val) );
        p4 =  16.0 * (*func0)( (a_val - h) );
        p5 = - 1.0 * (*func0)( (a_val - 2.0 * h) );
        denom = inv_0( 12.0 * h * h);
        deriv = (p1 + p2 + p3 + p4 + p5) * denom;

        if(mw_abs_0(result.hessian[0]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.hessian[0] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD HESSIAN! At (a=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, result.hessian[0], deriv, check);
        }
    }
    return failed;
}

static inline int testSingleVar_s(real (*func)(real*, real_0), real_0 (*func0)(real_0, real_0), real_0 min_a, real_0 max_a, real_0 min_b, real_0 max_b)
{
    int failed = 0;
    real_0 a_val, b_val;
    real a, result;
    real_0 h = 0.0001;
    real_0 deriv, check;
    real_0 p1, p2, p3, p4, p5, denom;

    srand(time(0));

    for (int i = 0; i < 100000; i++)
    {
        a_val = ((real_0)rand()/(real_0)RAND_MAX) * (max_a - min_a) + min_a;
        b_val = ((real_0)rand()/(real_0)RAND_MAX) * (max_b - min_b) + min_b;
        a = mw_real_var(a_val, 0);
        result = (*func)(&a, b_val);

        p1 =   1.0 * (*func0)( (a_val - 2.0 * h), b_val );
        p2 = - 8.0 * (*func0)( (a_val - h), b_val );
        p3 = - 1.0 * (*func0)( (a_val + 2.0 * h), b_val );
        p4 =   8.0 * (*func0)( (a_val + h), b_val );
        denom = inv_0( 12.0 * h);
        deriv = (p1 + p2 + p3 + p4) * denom;

        if(mw_abs_0(result.gradient[0]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.gradient[0] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD GRADIENT! At (a=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, result.gradient[0], deriv, check);
        }

        p1 = - 1.0 * (*func0)( (a_val + 2.0 * h), b_val );
        p2 =  16.0 * (*func0)( (a_val + h), b_val );
        p3 = -30.0 * (*func0)( (a_val), b_val );
        p4 =  16.0 * (*func0)( (a_val - h), b_val );
        p5 = - 1.0 * (*func0)( (a_val - 2.0 * h), b_val );
        denom = inv_0( 12.0 * h * h);
        deriv = (p1 + p2 + p3 + p4 + p5) * denom;

        if(mw_abs_0(result.hessian[0]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.hessian[0] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD HESSIAN! At (a=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, result.hessian[0], deriv, check);
        }
    }
    return failed;
}

static inline int testDoubleVar(real (*func)(real*, real*), real_0 (*func0)(real_0, real_0), real_0 min_a, real_0 max_a, real_0 min_b, real_0 max_b)
{
    int failed = 0;
    real_0 a_val, b_val;
    real a, b, result;
    real_0 h = 0.0001;
    real_0 deriv, check;
    real_0 p1, p2, p3, p4, p5, denom;
    real_0 p11, p12, p13, p14, p21, p22, p23, p24, p31, p32, p33, p34, p41, p42, p43, p44;

    srand(time(0));

    for (int i = 0; i < 100000; i++)
    {
        a_val = ((real_0)rand()/(real_0)RAND_MAX) * (max_a - min_a) + min_a;
        b_val = ((real_0)rand()/(real_0)RAND_MAX) * (max_b - min_b) + min_b;
        a = mw_real_var(a_val, 0);
        b = mw_real_var(b_val, 1);
        result = (*func)(&a, &b);

        p1 =   1.0 * (*func0)( (a_val - 2.0 * h), b_val );
        p2 = - 8.0 * (*func0)( (a_val - h), b_val );
        p3 = - 1.0 * (*func0)( (a_val + 2.0 * h), b_val );
        p4 =   8.0 * (*func0)( (a_val + h), b_val );
        denom = inv_0( 12.0 * h);
        deriv = (p1 + p2 + p3 + p4) * denom;

        if(mw_abs_0(result.gradient[0]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.gradient[0] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD GRADIENT(a)! At (a=%.15f, b=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, b_val, result.gradient[0], deriv, check);
        }

        p1 =   1.0 * (*func0)( a_val, (b_val - 2.0 * h) );
        p2 = - 8.0 * (*func0)( a_val, (b_val - h) );
        p3 = - 1.0 * (*func0)( a_val, (b_val + 2.0 * h) );
        p4 =   8.0 * (*func0)( a_val, (b_val + h) );
        denom = inv_0( 12.0 * h);
        deriv = (p1 + p2 + p3 + p4) * denom;

        if(mw_abs_0(result.gradient[1]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.gradient[1] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD GRADIENT(b)! At (a=%.15f, b=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, b_val, result.gradient[1], deriv, check);
        }
        p1 = - 1.0 * (*func0)( (a_val + 2.0 * h), b_val );
        p2 =  16.0 * (*func0)( (a_val + h), b_val );
        p3 = -30.0 * (*func0)( a_val, b_val );
        p4 =  16.0 * (*func0)( (a_val - h), b_val );
        p5 = - 1.0 * (*func0)( (a_val - 2.0 * h), b_val );
        denom = inv_0( 12.0 * h * h);
        deriv = (p1 + p2 + p3 + p4 + p5) * denom;

        if(mw_abs_0(result.hessian[0]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.hessian[0] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD HESSIAN(a,a)! At (a=%.15f, b=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, b_val, result.hessian[0], deriv, check);
        }

        p1 = - 1.0 * (*func0)( a_val, (b_val + 2.0 * h) );
        p2 =  16.0 * (*func0)( a_val, (b_val + h) );
        p3 = -30.0 * (*func0)( a_val, b_val );
        p4 =  16.0 * (*func0)( a_val, (b_val - h) );
        p5 = - 1.0 * (*func0)( a_val, (b_val - 2.0 * h) );
        denom = inv_0( 12.0 * h * h);
        deriv = (p1 + p2 + p3 + p4 + p5) * denom;

        if(mw_abs_0(result.hessian[2]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.hessian[2] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD HESSIAN(b,b)! At (a=%.15f, b=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, b_val, result.hessian[2], deriv, check);
        }

        p11 =   1.0 * (*func0)( (a_val - 2.0 * h), (b_val - 2.0 * h) );
        p12 = - 8.0 * (*func0)( (a_val - 2.0 * h), (b_val - h) );
        p13 = - 1.0 * (*func0)( (a_val - 2.0 * h), (b_val + 2.0 * h) );
        p14 =   8.0 * (*func0)( (a_val - 2.0 * h), (b_val + h) );
        p21 = - 8.0 * (*func0)( (a_val - h), (b_val - 2.0 * h) );
        p22 =  64.0 * (*func0)( (a_val - h), (b_val - h) );
        p23 =   8.0 * (*func0)( (a_val - h), (b_val + 2.0 * h) );
        p24 = -64.0 * (*func0)( (a_val - h), (b_val + h) );
        p31 = - 1.0 * (*func0)( (a_val + 2.0 * h), (b_val - 2.0 * h) );
        p32 =   8.0 * (*func0)( (a_val + 2.0 * h), (b_val - h) );
        p33 =   1.0 * (*func0)( (a_val + 2.0 * h), (b_val + 2.0 * h) );
        p34 = - 8.0 * (*func0)( (a_val + 2.0 * h), (b_val + h) );
        p41 =   8.0 * (*func0)( (a_val + h), (b_val - 2.0 * h) );
        p42 = -64.0 * (*func0)( (a_val + h), (b_val - h) );
        p43 = - 8.0 * (*func0)( (a_val + h), (b_val + 2.0 * h) );
        p44 =  64.0 * (*func0)( (a_val + h), (b_val + h) );
        denom = inv_0( 144.0 * h * h);
        deriv = (p11 + p12 + p13 + p14 + p21 + p22 + p23 + p24 + p31 + p32 + p33 + p34 + p41 + p42 + p43 + p44) * denom;

        if(mw_abs_0(result.hessian[1]) < small_derv)
        {
            check = mw_abs_0(deriv);
        }
        else
        {
            check = mw_abs_0(1.0 - result.hessian[1] / deriv);
        }
        if(check > derv_thresh)
        {
            failed += 1;
            printf("\t BAD HESSIAN(a,b)! At (a=%.15f, b=%.15f), Autodiff gave %.15f while Numerical gave %.15f (check=%.15f)\n", a_val, b_val, result.hessian[1], deriv, check);
        }
    }
    return failed;
}

int main()
{
    int failed = 0;
    int test = 0;

    //Arithmetic Operators
    test = testDoubleVar(mw_add, add_func, -100.0, 100.0, -100.0, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_add)\n");
    }
    failed += test;

    test = testSingleVar_s(mw_add_s, add_func, -100.0, 100.0, -100.0, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_add_s)\n");
    }
    failed += test;

    test = testDoubleVar(mw_sub, sub_func, -100.0, 100.0, -100.0, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_sub)\n");
    }
    failed += test;

    test = testDoubleVar(mw_mul, mul_func, -100.0, 100.0, -100.0, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_mul)\n");
    }
    failed += test;

    test = testSingleVar_s(mw_mul_s, mul_func, -100.0, 100.0, -100.0, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_mul_s)\n");
    }
    failed += test;

    test = testDoubleVar(mw_div, div_func, -100.0, 100.0, 1.0, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_div)\n");
    }
    failed += test;

    test = testSingleVar(mw_neg, neg_func, -100.0, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_neg)\n");
    }
    failed += test;

    //Trigonometric Functions
    test = testSingleVar(mw_sin, mw_sin_0, 0.0, 2*M_PI);
    if(test != 0)
    {
        printf("    Test failed! (mw_sin)\n");
    }
    failed += test;

    test = testSingleVar(mw_sinpi, sinpi_func, 0.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_sinpi)\n");
    }
    failed += test;

    test = testSingleVar(mw_cos, mw_cos_0, 0.0, 2*M_PI);
    if(test != 0)
    {
        printf("    Test failed! (mw_cos)\n");
    }
    failed += test;

    test = testSingleVar(mw_cospi, cospi_func, 0.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_cospi)\n");
    }
    failed += test;

    test = testSingleVar(mw_tan, mw_tan_0, -M_PI/2.0*0.95, M_PI/2.0*0.95);
    if(test != 0)
    {
        printf("    Test failed! (mw_tan)\n");
    }
    failed += test;

    test = testSingleVar(mw_tanpi, tanpi_func, -0.5*0.95, 0.5*0.95);
    if(test != 0)
    {
        printf("    Test failed! (mw_tanpi)\n");
    }
    failed += test;

    //Inverse Trigonometric Functions
    test = testSingleVar(mw_asin, mw_asin_0, 0.0, 0.95);
    if(test != 0)
    {
        printf("    Test failed! (mw_asin)\n");
    }
    failed += test;

    test = testSingleVar(mw_asinpi, asinpi_func, 0.0, 0.95);
    if(test != 0)
    {
        printf("    Test failed! (mw_asin)\n");
    }
    failed += test;

    test = testSingleVar(mw_acos, mw_acos_0, 0.0, 0.95);
    if(test != 0)
    {
        printf("    Test failed! (mw_acos)\n");
    }
    failed += test;

    test = testSingleVar(mw_acospi, acospi_func, 0.0, 0.95);
    if(test != 0)
    {
        printf("    Test failed! (mw_acospi)\n");
    }
    failed += test;

    test = testSingleVar(mw_atan, mw_atan_0, -1.0, 1.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_atan)\n");
    }
    failed += test;

    test = testSingleVar(mw_atanpi, atanpi_func, -1.0, 1.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_atanpi)\n");
    }
    failed += test;

    test = testDoubleVar(mw_atan2, mw_atan2_0, 0.05, 1.0, 0.05, 1.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_atan2)\n");
    }
    failed += test;

    test = testDoubleVar(mw_atan2pi, atan2pi_func, 0.05, 1.0, 0.05, 1.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_atan2pi)\n");
    }
    failed += test;

    //Exponential and Logarithmic Functions
    test = testSingleVar(mw_log, mw_log_0, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_log)\n");
    }
    failed += test;

    test = testSingleVar(mw_log2, mw_log2_0, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_log2)\n");
    }
    failed += test;  

    test = testSingleVar(mw_log10, mw_log10_0, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_log10)\n");
    }
    failed += test;

    test = testSingleVar(mw_log1p, mw_log1p_0, -0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_log1p)\n");
    }
    failed += test;

    test = testSingleVar(mw_exp, mw_exp_0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_exp)\n");
    }
    failed += test; 

    test = testSingleVar(mw_expm1, mw_expm1_0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_expm1)\n");
    }
    failed += test;

    test = testSingleVar(mw_exp2, mw_exp2_0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_exp2)\n");
    }
    failed += test;

    test = testSingleVar(mw_exp10, mw_exp10_0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_exp10)\n");
    }
    failed += test;

    //Hyperbolic Trig Functions
    test = testSingleVar(mw_sinh, mw_sinh_0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_sinh)\n");
    }
    failed += test;

    test = testSingleVar(mw_cosh, mw_cosh_0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_cosh)\n");
    }
    failed += test;

    test = testSingleVar(mw_tanh, mw_tanh_0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_tanh)\n");
    }
    failed += test;

    test = testSingleVar(mw_asinh, mw_asinh_0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_asinh)\n");
    }
    failed += test;

    test = testSingleVar(mw_acosh, mw_acosh_0, 1.01, 30.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_acosh)\n");
    }
    failed += test;

    test = testSingleVar(mw_atanh, mw_atanh_0, -0.95, 0.95);
    if(test != 0)
    {
        printf("    Test failed! (mw_atanh)\n");
    }
    failed += test;

    // Polynomial and Root Functions
    test = testDoubleVar(mw_pow, mw_pow_0, 0.0, 10.0, 0.1, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_pow)\n");
    }
    failed += test;

    test = testDoubleVar(mw_powr, powr_func, 0.0, 10.0, 0.1, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_powr)\n");
    }
    failed += test;

    test = testSingleVar(sqr, sqr_func, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (sqr)\n");
    }
    failed += test;

    test = testSingleVar(cube, cube_func, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (cube)\n");
    }
    failed += test;

    test = testSingleVar(fourth, fourth_func, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (fourth)\n");
    }
    failed += test;

    test = testSingleVar(fifth, fifth_func, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (fifth)\n");
    }
    failed += test;

    test = testSingleVar(sixth, sixth_func, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (sixth)\n");
    }
    failed += test;

    test = testSingleVar(inv, inv_func, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (inv)\n");
    }
    failed += test;

    test = testSingleVar(mw_cbrt, mw_cbrt_0, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_cbrt)\n");
    }
    failed += test;

    test = testSingleVar(mw_sqrt, mw_sqrt_0, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_sqrt)\n");
    }
    failed += test;

    test = testSingleVar(mw_rsqrt, rsqrt_func, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_rsqrt)\n");
    }
    failed += test;

    test = testSingleVar(minushalf, minushalf_func, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (minushalf)\n");
    }
    failed += test;

    test = testSingleVar(threehalves, threehalves_func, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (threehalves)\n");
    }
    failed += test;

    test = testSingleVar(fivehalves, fivehalves_func, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (fivehalves)\n");
    }
    failed += test;

    test = testSingleVar(minusthreehalves, minusthreehalves_func, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (minusthreehalves)\n");
    }
    failed += test;

    test = testSingleVar(minusfivehalves, minusfivehalves_func, 0.5, 100.0);
    if(test != 0)
    {
        printf("    Test failed! (minusfivehalves)\n");
    }
    failed += test;

    //Error Functions
    test = testSingleVar(mw_erfc, mw_erfc_0, -1.0, 1.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_erfc)\n");
    }
    failed += test;

    test = testSingleVar(mw_erf, mw_erf_0, -1.0, 1.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_erf)\n");
    }
    failed += test;

    //Bessel Functions
    test = testSingleVar(mw_besselJ0, besselJ0, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_besselJ0)\n");
    }
    failed += test;

    test = testSingleVar(mw_besselJ1, besselJ1, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_besselJ1)\n");
    }
    failed += test;

    test = testSingleVar(mw_besselJ2, besselJ2, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_besselJ2)\n");
    }
    failed += test;

    test = testSingleVar(mw_besselI0, besselI0, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_besselI0)\n");
    }
    failed += test;

    test = testSingleVar(mw_besselI1, besselI1, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_besselI1)\n");
    }
    failed += test;

    test = testSingleVar(mw_besselK0, besselK0, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_besselK0)\n");
    }
    failed += test;

    test = testSingleVar(mw_besselK1, besselK1, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_besselK1)\n");
    }
    failed += test;

    test = testSingleVar(mw_besselK2, besselK2, 0.5, 10.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_besselK2)\n");
    }
    failed += test;

    //Other Computer Functions
    test = testDoubleVar(mw_hypot, hypot_func, -2.0, 2.0, -2.0, 2.0);
    if(test != 0)
    {
        printf("    Test failed! (mw_hypot)\n");
    }
    failed += test;

    test = testSingleVar(d2r, d2r_func, 0.0, 360.0);
    if(test != 0)
    {
        printf("    Test failed! (d2r)\n");
    }
    failed += test;

    test = testSingleVar(r2d, r2d_func, 0.0, 2*M_PI);
    if(test != 0)
    {
        printf("    Test failed! (r2d)\n");
    }
    failed += test;


    if(failed == 0)
    {
        printf("All tests successful!\n");
    }

    return failed;
}
