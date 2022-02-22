/*
 * Copyright (c) 2010, 2011 Matthew Arsenault
 * Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute.
 * Copyright (c) 2016-2018 Siddhartha Shelton
 * Copyright (c) 2018 Eric Mendelsohn
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_priv.h"
#include "nbody_potential.h"
#include "nbody_potential_types.h"
#include "milkyway_util.h"
#include "nbody_caustic.h"
#include "nbody_bessel.h"

#include <time.h>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

/*Methods to be called by potentials*/

static inline real_0 lnfact(int n)
{
     int counter;
     real_0 result = 0.0;
     if (n > 0)
     {
          for (counter = n; counter >= 1; counter--)
          {
               result += mw_log_0((real_0) counter);
          }
     }
     /*mw_printf("ln(%u!) = %.15f \n",n,result);*/
     return result;
}

static inline real_0 binom(int n, int k)
{
    return mw_exp_0(lnfact(n)-lnfact(k)-lnfact(n-k));
}

static inline real leg_pol(real* x, int l)
{
    real sum = ZERO_REAL;
    real tmp1;
    for (int m = 0; m < mw_floor_0(l/2)+1; m++)
    {
        tmp1 = mw_real_const(l-2*m);
        tmp1 = mw_pow(&x, &tmp1);
        tmp1 = mw_mul_s(&tmp1, mw_pow_0(-1, m)/mw_pow_0(2,l)*binom(l,m)*binom(2*(l-m),l));
        sum = mw_add(&sum, &tmp1);
    }
    /*mw_printf("P(%.15f, %u) = %.15f \n",x,l,sum);*/
    return sum;
}

static inline real_0 leg_pol_zeros(int l)
{
    if (l % 2 == 0)
    {
        return mw_pow_0(-1,l/2)*mw_exp_0(lnfact(l)-2*lnfact(l/2) - l*mw_log_0(2));
    }
    else
    {
        return 0.0;
    }
}

static inline real leg_pol_derv(real* x, int l)
{
    real sum = ZERO_REAL;
    real tmp1;
    for (int m = 0; m < mw_floor_0((l-1)/2)+1; m++)
    {
        tmp1 = mw_real_const(l - 2*m - 1);
        tmp1 = mw_pow(&x, &tmp1);
        tmp1 = mw_mul_s(&tmp1, mw_pow_0(-1, m)/mw_pow_0(2,l)*binom(l,m)*binom(2*(l-m),l)*(l - 2*m));
        sum = mw_add(&sum, &tmp1);
    }
    /*mw_printf("P'(%.15f, %u) = %.15f \n",x,l,sum);*/
    return sum;
}

static inline real lower_gamma(int n, real* x)
{
    real sum = ZERO_REAL;
    real tmp1;
    for (int k = 0; k < n; k++)
    {
       tmp1 = mw_real_const(k);
       tmp1 = mw_pow(&x, &tmp1);
       tmp1 = mw_mul_s(&tmp1, inv_0(mw_exp_0(lnfact(k))));
       sum = mw_add(&sum, &tmp1);
    }
    /*mw_printf("g(%u, %.15f) = %.15f \n",n,x,mw_exp(lnfact(n-1))*(1-mw_exp(-x)*sum));*/
    tmp1 = mw_neg(&x);
    tmp1 = mw_exp(&tmp1);
    tmp1 = mw_mul(&tmp1, &sum);
    tmp1 = mw_neg(&tmp1);
    tmp1 = mw_add_s(&tmp1, 1.0);
    return mw_mul_s(&tmp1, mw_exp_0(lnfact(n-1)));
}

static inline real G_cont_frac(real* x, int n, int k)
{
    real tmp1;
    if (k > 15)
    {
        return mw_add_s(x, n+k);
    }
    else
    {
        tmp1 = G_cont_frac(x,n,k+1);
        tmp1 = inv(&tmp1);
        tmp1 = mw_mul_s(&tmp1, (1+k));
        tmp1 = mw_add_s(&tmp1, 1.0);
        tmp1 = inv(&tmp1);
        tmp1 = mw_mul_s(&tmp1, (n+k));
        return mw_add(x, &tmp1);
    }
}

static inline real GenExpIntegral(int n, real* x)
{
    real tmp1, tmp2;
    tmp1 = G_cont_frac(x,n,0);
    tmp2 = mw_exp(x);
    tmp1 = mw_mul(&tmp1, &tmp2);
    return inv(&tmp1);
}

//static inline real RExpIntegrand (real k, real R, real Rd, real z, real zd)
//{
//    real val = k*mw_cos(k*z)*(aExp(k,R,Rd)*besselK1(k*R) - bExp(k,R,Rd)*besselI1(k*R))/(sqr(zd*k) + 1);
    //mw_printf("RExp(%.15f,%.15f,%.15f,%.15f,%.15f) = %.15f\n",k,R,Rd,z,zd,val);
//    return val;
//}

//static inline real ZExpIntegrand (real k, real R, real Rd, real z, real zd)
//{
//    real val = k*mw_sin(k*z)*(aExp(k,R,Rd)*besselK0(k*R) + bExp(k,R,Rd)*besselI0(k*R))/(sqr(zd*k) + 1);
    //mw_printf("ZExp = %.15f\n",val);
//    return val;
//}

static inline real RSechIntegrand (real* k, real* R, real* Rd, real* z, real* zd)
{
    real tmp1, tmp2, tmp3;

    tmp1 = aExp(k,R,Rd);
    tmp2 = mw_mul(k,R);
    tmp2 = mw_besselK1(&tmp2);
    tmp1 = mw_mul(&tmp1,&tmp2);
    tmp2 = bExp(k,R,Rd);
    tmp3 = mw_mul(k,R);
    tmp3 = mw_besselI1(&tmp3);
    tmp2 = mw_mul(&tmp2,&tmp3);
    real ExpStuff = mw_sub(&tmp1, &tmp2);

    tmp1 = mw_mul(k,z);
    tmp1 = mw_cos(&tmp1);
    tmp2 = mw_mul(k,zd);
    tmp2 = mw_mul_s(&tmp2,M_PI/2.0);
    tmp2 = mw_sinh(&tmp2);
    real TrigStuff = mw_div(&tmp1, &tmp2);

    tmp1 = sqr(k);
    tmp1 = mw_mul(&tmp1,zd);
    tmp1 = mw_mul(&tmp1,&TrigStuff);
    return mw_mul(&tmp1,&ExpStuff);
    //mw_printf("RSech(%.15f,%.15f,%.15f,%.15f,%.15f) = %.15f\n",k,R,Rd,z,zd,val);
}

static inline real ZSechIntegrand (real* k, real* R, real* Rd, real* z, real* zd)
{
    real tmp1, tmp2, tmp3;

    tmp1 = aExp(k,R,Rd);
    tmp2 = mw_mul(k,R);
    tmp2 = mw_besselK0(&tmp2);
    tmp1 = mw_mul(&tmp1,&tmp2);
    tmp2 = bExp(k,R,Rd);
    tmp3 = mw_mul(k,R);
    tmp3 = mw_besselI0(&tmp3);
    tmp2 = mw_mul(&tmp2,&tmp3);
    real ExpStuff = mw_sub(&tmp1, &tmp2);

    tmp1 = mw_mul(k,z);
    tmp1 = mw_sin(&tmp1);
    tmp2 = mw_mul(k,zd);
    tmp2 = mw_mul_s(&tmp2,M_PI/2.0);
    tmp2 = mw_sinh(&tmp2);
    real TrigStuff = mw_div(&tmp1, &tmp2);

    tmp1 = sqr(k);
    tmp1 = mw_mul(&tmp1,zd);
    tmp1 = mw_mul(&tmp1,&TrigStuff);
    return mw_mul(&tmp1,&ExpStuff);
    //mw_printf("ZSech = %.15f\n",val);
}

/**********************************************************************************************************************************************************************************************/

mwvector pointAccel(const mwvector* pos, const mwvector* pos1, const real* mass)
{
    mwvector acc;
    mwvector v = mw_subv(pos1, pos);
    real dist = mw_length(&v);
    real tmp = cube(&dist);
    tmp = mw_div(mass, &tmp);

    acc.x = mw_mul(&v.x, &tmp);
    acc.y = mw_mul(&v.y, &tmp);
    acc.z = mw_mul(&v.z, &tmp);

    return acc;
}

mwvector plummerAccel(const mwvector* pos, const mwvector* pos1, const real* mass, const real* scale)
{
    mwvector acc;
    mwvector v = mw_subv(pos1, pos);
    real dist = mw_length(&v);
    real tmp = mw_hypot(scale, &dist);
    tmp = cube(&tmp);
    tmp = mw_div(mass, &tmp);

    acc.x = mw_mul(&v.x, &tmp);
    acc.y = mw_mul(&v.y, &tmp);
    acc.z = mw_mul(&v.z, &tmp);

    return acc;
}


/*spherical bulge potentials*/

static inline mwvector hernquistSphericalAccel(const Spherical* sph, mwvector* pos, real* r)
{
    mwvector acc;
    const real M = mw_real_var(sph->mass, BULGE_MASS_POS);
    const real a = mw_real_var(sph->scale, BULGE_RADIUS_POS);
    real tmp = mw_add(&a, r);
    real neg_r = mw_neg(r);

    tmp = sqr(&tmp);
    tmp = mw_mul(&neg_r, &tmp);
    tmp = mw_div(&M, &tmp);
    //mw_printf("TMP = %.15f\n", showRealValue(&tmp));

    acc.x = mw_mul(&pos->x, &tmp);
    acc.y = mw_mul(&pos->y, &tmp);
    acc.z = mw_mul(&pos->z, &tmp);
    //mw_printf("ACC = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acc)), showRealValue(&Y(&acc)), showRealValue(&Z(&acc)));

    return acc;
}

static inline mwvector plummerSphericalAccel(const Spherical* sph, mwvector* pos, real* r)
{
    mwvector acc;
    real M = mw_real_var(sph->mass, BULGE_MASS_POS);
    real a = mw_real_var(sph->scale, BULGE_RADIUS_POS);
    real tmp = mw_hypot(&a,r);

    tmp = cube(&tmp);
    tmp = mw_div(&M, &tmp);
    tmp = mw_neg(&tmp);

    acc.x = mw_mul(&pos->x, &tmp);
    acc.y = mw_mul(&pos->y, &tmp);
    acc.z = mw_mul(&pos->z, &tmp);

    return acc;
}

/* Disk Potentials */

static inline mwvector miyamotoNagaiDiskAccel(const Disk* disk, mwvector* pos)
{
    mwvector acc;
    real tmp1;
    const real M   = mw_real_var(disk->mass, DISK_MASS_POS);
    const real a   = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);
    const real b   = mw_real_var(disk->scaleHeight, DISK_HEIGHT_POS);
    const real zp  = mw_hypot(&Z(pos), &b);
    const real azp = mw_add(&a, &zp);

    tmp1 = mw_hypot(&X(pos), &Y(pos));
    const real rp  = mw_hypot(&tmp1, &azp);
    const real rth = cube(&rp);  /* rp ^ (3) */

    tmp1 = mw_mul(&M, &X(pos));
    tmp1 = mw_div(&tmp1, &rth);
    X(&acc) = mw_neg(&tmp1);

    tmp1 = mw_mul(&M, &Y(pos));
    tmp1 = mw_div(&tmp1, &rth);
    Y(&acc) = mw_neg(&tmp1);

    tmp1 =  mw_mul(&zp, &rth);
    tmp1 = mw_div(&azp, &tmp1);
    tmp1 = mw_mul(&Z(pos), &tmp1);
    tmp1 = mw_mul(&M, &tmp1);
    Z(&acc) = mw_neg(&tmp1);

    //mw_printf("Acceleration[AX,AY,AZ] = [%.15f,%.15f,%.15f]\n",X(acc),Y(acc),Z(acc));
    
    /*mw_printf("disk x acc: %.15f\n", X(acc));
    mw_printf("disk y acc: %.15f\n", Y(acc));
    mw_printf("disk z acc: %.15f\n", Z(acc));*/
    return acc;
}

/*WARNING: This potential is currently not in use as it is non-physical and difficult to test.*/
static inline mwvector freemanDiskAccel(const Disk* disk, mwvector* pos)     /*From Freeman 1970*/
{
    mwvector acc;
    /*Reset acc vector*/
    X(&acc) = ZERO_REAL;
    Y(&acc) = ZERO_REAL;
    Z(&acc) = ZERO_REAL;

//    const mwvector r_hat = mw_mulvs(pos, inv(r));
//    const real r_proj = mw_hypot(X(pos), Y(pos));

//    mwvector theta_hat;
//    X(theta_hat) = X(pos)*Z(pos)/r/r_proj;
//    Y(theta_hat) = Y(pos)*Z(pos)/r/r_proj;
//    Z(theta_hat) = -r_proj/r;

//    const real scl = disk->scaleLength;
//    const real M = disk->mass;
//    const real costheta = Z(pos)/r;
//    const real sintheta = r_proj/r;

//    real r_f = 0.0;
//    real theta_f = 0.0;
//    real a_r = 0.0;
//    real b_r = 0.0;
//    real p0 = 0.0;
//    real pcos = 0.0;

//    for (int l = 0; l < 5; l++) /*Only even terms add, so counter is multiplied by 2*/
//    {
//        if (l>0)
//        {
//            a_r = GenExpIntegral(2*l,r/scl);
//        }
//        else
//        {
//            a_r = 0.0;
//        }
//        b_r = mw_pow(scl/r,2*l+2)*lower_gamma(2*l+2,r/scl);
//        p0 = leg_pol_zeros(2*l);
//        pcos = leg_pol(costheta,2*l);
//        r_f += p0*pcos*(2*l*a_r - (2*l+1)*b_r);

//        if (l > 0)
//        {
//            theta_f -= p0*(a_r+b_r)*sintheta*leg_pol_derv(costheta,2*l);
//        }
//    }

//    mwvector r_comp = mw_mulvs(r_hat, r_f*M/sqr(scl));
//    mwvector theta_comp = mw_mulvs(theta_hat, theta_f*M/sqr(scl));

//    X(acc) = X(r_comp) + X(theta_comp);
//    Y(acc) = Y(r_comp) + Y(theta_comp);
//    Z(acc) = Z(r_comp) + Z(theta_comp);

    //if (r_f > 0)
    //{
    //    mw_printf("ERROR: Repulsive acceleration!\n");
    //    mw_printf("r_magnitude = %.15f\n",r_f*M/sqr(scl));
    //    mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    //    mw_printf("Acceleration = %.15f \n",mw_absv(acc));
    //}

    return acc;
}

/*WARNING: This potential is currently failing the Poisson test. Use at your own risk!*/
static inline mwvector doubleExponentialDiskAccel(const Disk* disk, mwvector* pos)
{
    //mw_printf("Calculating Acceleration\n");
    //mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    mwvector acc;

    const real R = mw_hypot(&X(pos), &Y(pos));
    mwvector R_hat;
    R_hat.x = mw_div(&X(pos), &R);
    R_hat.y = mw_div(&Y(pos), &R);
    R_hat.z = ZERO_REAL;

    mwvector Z_hat;
    Z_hat.x = ZERO_REAL;
    Z_hat.y = ZERO_REAL;
    Z_hat.z = mw_real_const(1.0);
    //mw_printf("Z_hat = [ %.15f, %.15f, %.15f ]\n", showRealValue(&Z_hat.x), showRealValue(&Z_hat.y), showRealValue(&Z_hat.z));

    const real M    = mw_real_var(disk->mass, DISK_MASS_POS);
    const real Rd   = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);
    const real zd   = mw_real_var(disk->scaleHeight, DISK_HEIGHT_POS);
    const real z    = Z(pos);
    const real_0 h  = 0.001;

    const real a = inv(&Rd);
    const real b = inv(&zd);
    
    real R_piece = ZERO_REAL;
    real z_piece = ZERO_REAL;

    real_0 j0_zero;
    real_0 psi_in_0;
    real_0 psi_0;
    real_0 psi_prime_0;
    real_0 j0_x;
    real_0 j0_w;
    real fun_0;

    real_0 j1_zero;
    real_0 psi_in_1;
    real_0 psi_1;
    real_0 psi_prime_1;
    real_0 j1_x;
    real_0 j1_w;
    real fun_1;

    real tmp1, tmp2, tmp3, tmp4;

    for (int n = 0; n < 150; n+=1)    // Hidenori Ogata's Numerical Integration Formula Based on the Bessel Functions
    {
        j0_zero = besselJ0_zero(n)/M_PI;
        j1_zero = besselJ1_zero(n)/M_PI;
        psi_in_0 = h * j0_zero;
        psi_in_1 = h * j1_zero;
        psi_0 = psi_in_0 * mw_sinh_0(M_PI / 2.0 * mw_sinh_0(psi_in_0)) / mw_cosh_0(M_PI / 2.0 * mw_sinh_0(psi_in_0));
        psi_1 = psi_in_1 * mw_sinh_0(M_PI / 2.0 * mw_sinh_0(psi_in_1)) / mw_cosh_0(M_PI / 2.0 * mw_sinh_0(psi_in_1));
        psi_prime_0 = (mw_sinh_0(M_PI * mw_sinh_0(psi_in_0)) + M_PI * psi_in_0 * mw_cosh_0(psi_in_0)) / (mw_cosh_0(M_PI * mw_sinh_0(psi_in_0)) + 1.0);
        psi_prime_1 = (mw_sinh_0(M_PI * mw_sinh_0(psi_in_1)) + M_PI * psi_in_1 * mw_cosh_0(psi_in_1)) / (mw_cosh_0(M_PI * mw_sinh_0(psi_in_1)) + 1.0);
        j0_x = M_PI / h * psi_0;
        j1_x = M_PI / h * psi_1;

        j0_w = 2.0 / (M_PI * j0_zero * sqr_0(besselJ1(M_PI * j0_zero))) * besselJ0(j0_x) * psi_prime_0;
        j1_w = 2.0 / (M_PI * j1_zero * sqr_0(besselJ2(M_PI * j1_zero))) * besselJ1(j1_x) * psi_prime_1;

        tmp1 = inv(&R);
        real j0_x_R = mw_mul_s(&tmp1, j0_x);
        real j1_x_R = mw_mul_s(&tmp1, j1_x);

        tmp1 = mw_hypot(&a, &j0_x_R);
        tmp1 = cube(&tmp1);
        tmp1 = inv(&tmp1);
        tmp2 = mw_abs(&z);
        tmp2 = mw_mul(&j0_x_R, &tmp2);
        tmp2 = mw_exp(&tmp2);
        tmp2 = inv(&tmp2);
        tmp3 = mw_abs(&z);
        tmp3 = mw_mul(&b, &tmp3);
        tmp3 = mw_exp(&tmp3);
        tmp3 = inv(&tmp3);
        tmp2 = mw_sub(&tmp2, &tmp3);
        tmp3 = sqr(&b);
        tmp4 = sqr(&j0_x_R);
        tmp3 = mw_sub(&tmp3, &tmp4);
        tmp2 = mw_div(&tmp2, &tmp3);
        tmp1 = mw_mul(&tmp1, &tmp2);
        fun_0 = mw_mul(&tmp1, &j0_x_R);

        tmp1 = mw_hypot(&a, &j1_x_R);
        tmp1 = cube(&tmp1);
        tmp1 = inv(&tmp1);
        tmp2 = mw_abs(&z);
        tmp2 = mw_mul(&j1_x_R, &tmp2);
        tmp2 = mw_exp(&tmp2);
        tmp2 = inv(&tmp2);
        tmp2 = mw_mul(&b, &tmp2);
        tmp3 = mw_abs(&z);
        tmp3 = mw_mul(&b, &tmp3);
        tmp3 = mw_exp(&tmp3);
        tmp3 = inv(&tmp3);
        tmp3 = mw_mul(&j1_x_R, &tmp3);
        tmp2 = mw_sub(&tmp2, &tmp3);
        tmp3 = sqr(&b);
        tmp4 = sqr(&j1_x_R);
        tmp3 = mw_sub(&tmp3, &tmp4);
        tmp2 = mw_div(&tmp2, &tmp3);
        tmp1 = mw_mul(&tmp1, &tmp2);
        fun_1 = mw_mul_s(&tmp1, j1_x);

        tmp1 = mw_mul(&a, &b);
        tmp1 = mw_div(&tmp1, &R);
        tmp1 = mw_mul(&tmp1, &fun_0);
        tmp1 = mw_mul_s(&tmp1, j0_w);
        //mw_printf("fun_0 = %.15f\n", showRealValue(&fun_0));
        //mw_printf("j0_w  = %.15f\n", j0_w);
        //mw_printf("tmp1  = %.15f\n", showRealValue(&tmp1));
        real z_pieceAdd = mw_mul_s(&tmp1, 4.0 * M_PI);


        tmp1 = sqr(&R);
        tmp1 = mw_div(&a, &tmp1);
        tmp1 = mw_mul(&tmp1, &fun_1);
        tmp1 = mw_mul_s(&tmp1, j1_w);
        real R_pieceAdd = mw_mul_s(&tmp1, 4.0 * M_PI);
        
        z_piece = mw_add(&z_piece, &z_pieceAdd);
        R_piece = mw_add(&R_piece, &R_pieceAdd);
    }

    tmp1 = sqr(&a);
    tmp1 = mw_mul(&M, &tmp1);
    tmp1 = mw_mul(&tmp1, &b);
    tmp1 = mw_mul(&tmp1, &R_piece);
    tmp1 = mw_mul_s(&tmp1, -1.0 / 4.0 / M_PI);
    mwvector R_comp;
    R_comp.x = mw_mul(&R_hat.x, &tmp1);
    R_comp.y = mw_mul(&R_hat.y, &tmp1);
    R_comp.z = ZERO_REAL;
    //mw_printf("R_comp = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&R_comp)), showRealValue(&Y(&R_comp)), showRealValue(&X(&R_comp)));

    tmp1 = sqr(&a);
    tmp1 = mw_mul(&M, &tmp1);
    tmp1 = mw_mul(&tmp1, &b);
    tmp1 = mw_mul(&tmp1, &z_piece);
    tmp1 = mw_mul_s(&tmp1, -1.0 / 4.0 / M_PI);
    //mw_printf("tmp1  = %.15f\n", showRealValue(&tmp1));
    mwvector Z_comp;
    Z_comp.x = ZERO_REAL;
    Z_comp.y = ZERO_REAL;
    Z_comp.z = mw_mul(&Z_hat.z, &tmp1);
    //mw_printf("Z_comp = [%.15f,%.15f,%.15f]\n", showRealValue(&Z_comp.x), showRealValue(&Z_comp.y), showRealValue(&Z_comp.z));

    X(&acc) = mw_add(&X(&R_comp), &X(&Z_comp));
    Y(&acc) = mw_add(&Y(&R_comp), &Y(&Z_comp));
    Z(&acc) = mw_add(&Z(&R_comp), &Z(&Z_comp));

    //real magnitude = mw_hypot(mw_hypot(X(acc), Y(acc)), Z(acc));
    //mw_printf("acc = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acc)), showRealValue(&Y(&acc)), showRealValue(&X(&acc)));
    return acc;
}

/*WARNING: This potential can take a while to integrate if any part of the orbit extends past 100 times the scaleLength*/
static inline mwvector sech2ExponentialDiskAccel(const Disk* disk, mwvector* pos)
{
    //mw_printf("Calculating Acceleration\n");
    //mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    mwvector acc;

    const real R = mw_hypot(&X(pos), &Y(pos));
    mwvector R_hat;
    X(&R_hat) = mw_div(&X(pos), &R);
    Y(&R_hat) = mw_div(&Y(pos), &R);
    Z(&R_hat) = ZERO_REAL;

    mwvector Z_hat;
    X(&Z_hat) = ZERO_REAL;
    Y(&Z_hat) = ZERO_REAL;
    Z(&Z_hat) = mw_real_const(1.0);

    const real M    = mw_real_var(disk->mass, DISK_MASS_POS);
    const real Rd   = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);
    const real zd   = mw_real_var(disk->scaleHeight, DISK_HEIGHT_POS);
    const real z = Z(pos);
    real tmp;

    const int n = 15;
    const real_0 a = 0.0;
    const real_0 b = 60.0 / showRealValue(&R); // Should be infinity, but this should be the point where the integral becomes negligible
    real integralR = ZERO_REAL;
    real integralZ = ZERO_REAL;
    const real_0 weight[] = {0.236927,0.478629,0.568889,0.478629,0.236927};
    const real_0 point[] = {-0.90618,-0.538469,0.0,0.538469,0.90618};
    const real_0 h = (b-a)/(n*1.0);

    for (int k = 0; k < n; k++)     /*Five-point Gaussian Quadrature*/
    {
        real Rpiece = ZERO_REAL;
        real Zpiece = ZERO_REAL;
        real k_val  = ZERO_REAL;
        for (int j = 0; j < 5; j++)
        {
            k_val = mw_real_const(h*point[j]/2 + a + (k*1.0+0.5)*h);

            tmp = RSechIntegrand(&k_val,&R,&Rd,&z,&zd);
            tmp = mw_mul_s(&tmp, h*weight[j]/2.0);
            Rpiece = mw_add(&Rpiece, &tmp);

            tmp = ZSechIntegrand(&k_val,&R,&Rd,&z,&zd);
            tmp = mw_mul_s(&tmp, h*weight[j]/2.0);
            Zpiece = mw_add(&Zpiece, &tmp);
        }
        integralR  = mw_add(&integralR, &Rpiece);
        integralZ  = mw_add(&integralZ, &Zpiece);
    }

    tmp = mw_mul(&M, &integralR);
    tmp = mw_neg(&tmp);
    mwvector R_comp;
    R_comp.x = mw_mul(&R_hat.x, &tmp);
    R_comp.y = mw_mul(&R_hat.y, &tmp);
    R_comp.z = mw_mul(&R_hat.z, &tmp);

    tmp = mw_mul(&M, &integralZ);
    tmp = mw_neg(&tmp);
    mwvector Z_comp;
    Z_comp.x = mw_mul(&Z_hat.x, &tmp);
    Z_comp.y = mw_mul(&Z_hat.y, &tmp);
    Z_comp.z = mw_mul(&Z_hat.z, &tmp);

    X(&acc) = mw_add(&X(&R_comp), &X(&Z_comp));
    Y(&acc) = mw_add(&Y(&R_comp), &Y(&Z_comp));
    Z(&acc) = mw_add(&Z(&R_comp), &Z(&Z_comp));

    //real magnitude = mw_sqrt(sqr(X(acc))+sqr(Y(acc))+sqr(Z(acc)));

    //mw_printf("Acceleration[AX,AY,AZ] = [%.15f,%.15f,%.15f]   Magnitude = %.15f\n",X(acc),Y(acc),Z(acc),magnitude);

    return acc;
}

//Softened needle bar potential
static inline mwvector orbitingBarAccel(const Disk* disk, mwvector* pos, real_0 time)
{
    real amp = mw_real_var(disk->mass, DISK_MASS_POS);
    real a = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);
    real b = mw_real_const(1.4);                 //Triaxial softening length
    real c = mw_real_const(1.0);                 //Prolate softening length
    real tmp1, tmp2, tmp3;
    
    real_0 curAngle = (disk->patternSpeed * time * -1)+disk->startAngle;
    //first rotate pos curAngle * -1 radians to emulate the current angle of the bar
    real Radi = mw_hypot(&X(pos), &Y(pos));

    tmp1 = mw_div(&Y(pos), &X(pos));
    real Phi = mw_atan(&tmp1);
    Phi = mw_add_s(&Phi, -curAngle);
    if(showRealValue(&X(pos)) < 0){
        Radi = mw_neg(&Radi);
    }
    tmp1 = mw_cos(&Phi);
    real x = mw_mul(&Radi, &tmp1);
    tmp1 = mw_sin(&Phi);
    real y = mw_mul(&Radi, &tmp1); 
    real z = Z(pos);

    //calculate force in accordance with the galpy implementation
    tmp1 = mw_hypot(&z,&c);
    tmp1 = mw_add(&b, &tmp1);
    real secondpart = mw_hypot(&y, &tmp1);

    tmp1 = mw_add(&a,&x);
    real Tp = mw_hypot(&tmp1, &secondpart);
    tmp1 = mw_sub(&a,&x);
    real Tm = mw_hypot(&tmp1, &secondpart);

    tmp1 = mw_add(&Tp, &Tm);
    tmp2 = sqr(&x);
    tmp2 = mw_div(&tmp2, &tmp1);
    tmp2 = mw_mul_s(&tmp2,4.0);
    real thirdpart = mw_sub(&tmp1, &tmp2);

    mwvector force;

    tmp1 = mw_mul(&Tp, &Tm);
    tmp2 = mw_add(&Tp, &Tm);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp1 = mw_div(&x, &tmp1);
    force.x = mw_mul_s(&tmp1, -2.0);

    tmp1 = mw_div(&y, &Tp);
    tmp1 = mw_div(&tmp1, &Tm);
    tmp1 = mw_mul(&tmp1, &thirdpart);
    tmp2 = sqr(&secondpart);
    tmp1 = mw_div(&tmp1, &tmp2);
    force.y = mw_mul_s(&tmp1, -0.5);

    tmp1 = mw_mul(&force.y, &z);
    tmp1 = mw_div(&tmp1, &y);
    tmp2 = mw_hypot(&z,&c);
    tmp2 = mw_add(&b, &tmp2);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = mw_hypot(&z,&c);
    force.z = mw_div(&tmp1, &tmp2);
    
    //undo the pos rotation and calculate acceleration from the force vector we got
    mwvector acc;
    real_0 cp = mw_cos_0(curAngle);
    real_0 sp = mw_sin_0(curAngle);

    tmp1 = mw_mul_s(&force.x, cp);
    tmp2 = mw_mul_s(&force.y, sp);
    acc.x = mw_sub(&tmp1, &tmp2);

    tmp1 = mw_mul_s(&force.x, sp);
    tmp2 = mw_mul_s(&force.y, cp);
    acc.y = mw_add(&tmp1, &tmp2);

    acc.z = force.z;
    
    acc.x = mw_mul(&acc.x, &amp);
    acc.y = mw_mul(&acc.y, &amp);
    acc.z = mw_mul(&acc.z, &amp);

    return acc;
}


/*If you want to test the time dependency of the bar with the bar as a point mass, comment out 
the above function and uncomment the one below*/

/*
static inline mwvector orbitingBarAccel(const Disk* disk, mwvector* pos, real* r, real_0 time)
{
    //mw_printf("Calculating Acceleration\n");
    //mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",showRealValue(&X(pos)),showRealValue(&Y(pos)),showRealValue(&Z(pos)));
    //mw_printf("r = %.15f\n", showRealValue(&r));

    mwvector pointPos;
    pointPos.z = ZERO_REAL;
    real_0 curAngle = (disk->patternSpeed * time * -1)+disk->startAngle;
    curAngle = curAngle - M_PI;//this is because the sun is negative in our coordinate system
    real diskLength = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);
    real diskMass = mw_real_var(disk->mass, DISK_MASS_POS);
    pointPos.x = mw_mul_s(&diskLength, mw_cos_0(curAngle)); //this is assuming top-down
    pointPos.y = mw_mul_s(&diskLength, mw_sin_0(curAngle));

    real dist = mw_distv(pos, &pointPos);

    mwvector acc = mw_subv(pos, &pointPos);  //get direction from pos to pointPos
    real tmp = cube(&dist);
    real totalAcc = mw_div(&diskMass, &tmp);
    acc = mw_mulvs(acc, totalAcc);           //a = Gm/r^2

    //mw_printf("curAngle: %.15f\n", curAngle);
    //mw_printf("pointPos: [%.15f,%.15f,%.15f]\n", showRealValue(&X(&pointPos)),showRealValue(&Y(&pointPos)),showRealValue(&Z(&pointPos)));
    //mw_printf("Accel: [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acc)),showRealValue(&Y(&acc)),showRealValue(&Z(&acc)));
    return acc;
}
*/

static inline mwvector logHaloAccel(const Halo* halo, mwvector* pos)
{
    mwvector acc;
    real tmp1, tmp2;

    const real v0 = mw_real_var(halo->vhalo, HALO_MASS_POS);
    const real d  = mw_real_var(halo->scaleLength, HALO_RADIUS_POS);
    const real q  = mw_real_var(halo->flattenZ, HALO_ZFLATTEN_POS);

    tmp1 = mw_hypot(&d, &X(pos));
    tmp1 = mw_hypot(&tmp1, &Y(pos));
    tmp2 = mw_div(&Z(pos), &q);
    const real denom = mw_hypot(&tmp1, &tmp2);

    tmp1 = sqr(&v0);
    tmp1 = mw_div(&tmp1, &denom);
    const real k = mw_mul_s(&tmp1, -2.0);

    X(&acc) = mw_mul(&k, &X(pos));
    Y(&acc) = mw_mul(&k, &Y(pos));

    tmp1 = sqr(&q);
    tmp1 = mw_div(&Z(pos), &tmp1);
    Z(&acc) = mw_mul(&k, &tmp1);

    return acc;
}

static inline mwvector nfwHaloAccel(const Halo* h, mwvector* pos, real* r)
{
    mwvector acc;
    real tmp1, tmp2, tmp3;
    const real v0 = mw_real_var(h->vhalo, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);

    tmp1 = sqr(&v0);
    tmp1 = mw_mul(&tmp1, &a);
    const real M = mw_mul_s(&tmp1, inv_0(0.2162165954)); /*Maximum of [ln(1+x)/x - 1/(1+x)]*/
    const real ar = mw_add(&a, r);

    tmp1 = mw_neg(&M);
    tmp2 = sqr(r);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_div(&ar, &a);
    tmp2 = mw_log(&tmp2);
    tmp2 = mw_div(&tmp2, r);
    tmp3 = inv(&ar);
    tmp2 = mw_sub(&tmp2, &tmp3);
    const real c = mw_mul(&tmp1, &tmp2);

    acc.x = mw_mul(&pos->x, &c);
    acc.y = mw_mul(&pos->y, &c);
    acc.z = mw_mul(&pos->z, &c);

    return acc;
}

/* CHECKME: Seems to have precision related issues for a small number of cases for very small qy */
static inline mwvector triaxialHaloAccel(const Halo* h, mwvector* pos)  /** Triaxial Logarithmic **/
{
    mwvector acc;
    real tmp1, tmp2, tmp3;

    /* TODO: More things here can be cached */
    tmp1 = mw_real_var(h->flattenZ, HALO_ZFLATTEN_POS);
    const real qzs      = sqr(&tmp1);

    tmp1 = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real rhalosqr = sqr(&tmp1);

    tmp1 = mw_real_var(h->vhalo, HALO_MASS_POS);
    tmp1 = sqr(&tmp1);
    const real mvsqr = mw_neg(&tmp1);

    const real xsqr = sqr(&X(pos));
    const real ysqr = sqr(&Y(pos));
    const real zsqr = sqr(&Z(pos));

    tmp1 = mw_mul_s(&xsqr, h->c1);
    tmp2 = mw_mul(&X(pos), &Y(pos));
    tmp2 = mw_mul_s(&tmp2, h->c3);
    tmp3 = mw_mul_s(&ysqr, h->c2);
    tmp2 = mw_add(&tmp2, &tmp3);
    tmp1 = mw_add(&tmp1, &tmp2);
    const real arst  = mw_add(&rhalosqr, &tmp1);

    tmp1 = mw_div(&zsqr, &qzs);
    const real arst2 = mw_add(&tmp1, &arst);

    tmp1 = mw_mul_s(&X(pos), (2.0 * h->c1));
    tmp2 = mw_mul_s(&Y(pos), h->c3);
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp1 = mw_div(&tmp1, &arst2);
    X(&acc) = mw_mul(&mvsqr, &tmp1);

    tmp1 = mw_mul_s(&Y(pos), (2.0 * h->c2));
    tmp2 = mw_mul_s(&X(pos), h->c3);
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp1 = mw_div(&tmp1, &arst2);
    Y(&acc) = mw_mul(&mvsqr, &tmp1);

    tmp1 = mw_mul(&mvsqr, &Z(pos));
    tmp1 = mw_mul_s(&tmp1, 2.0);
    tmp2 = mw_mul(&qzs, &arst);
    tmp2 = mw_add(&tmp2, &zsqr);
    Z(&acc) = mw_div(&tmp1, &tmp2);

    return acc;
}

static inline mwvector ASHaloAccel(const Halo* h, mwvector* pos, real* r)
{
    mwvector acc;
    real tmp1, tmp2;
    const real_0 gam = h->gamma;
    const real_0 lam = h->lambda;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real scaleR = mw_div(r, &a);
    tmp1 = inv(&a);
    const real scaleL = mw_mul_s(&tmp1, lam);
    real factor;
    real c;

    if (showRealValue(&r) < lam)
    {
        tmp1 = mw_mul(&a, r);
        tmp1 = mw_div(&M, &tmp1);
        factor = mw_neg(&tmp1);

        tmp1 =  mw_real_const(gam-1.0);
        tmp1 = mw_pow(&scaleR, &tmp1);
        tmp2 = mw_real_const(gam-1.0);
        tmp2 = mw_pow(&scaleR, &tmp2);
        tmp2 = mw_add_s(&tmp2, 1.0);
        tmp1 = mw_div(&tmp1, &tmp2);
        c = mw_mul(&factor, &tmp1);
        //c = -(M/(a*r))*mw_pow(scaleR,gam-1.0)/(1.0+mw_pow(scaleR,gam-1.0));
    }
    else
    {
        tmp1 = sqr(r);
        tmp1 = mw_div(&M, &tmp1);
        factor = mw_neg(&tmp1);

        tmp1 = mw_real_const(gam);
        tmp1 = mw_pow(&scaleL, &tmp1);
        tmp2 = mw_real_const(gam-1.0);
        tmp2 = mw_pow(&scaleL, &tmp2);
        tmp2 = mw_add_s(&tmp2, 1.0);
        tmp1 = mw_div(&tmp1, &tmp2);
        c = mw_mul(&factor, &tmp1);
        //c = -(M/sqr(r))*mw_pow(scaleL,gam)/(1.0+mw_pow(scaleL,gam-1.0));
    }

    tmp1 = mw_div(&c, r);

    acc.x = mw_mul(&pos->x, &tmp1);
    acc.y = mw_mul(&pos->y, &tmp1);
    acc.z = mw_mul(&pos->z, &tmp1);

    return acc;
}

static inline mwvector WEHaloAccel(const Halo* h, mwvector* pos, real* r)
{
    mwvector acc;
    real tmp1, tmp2, tmp3, tmp4;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real hypot = mw_hypot(&a, r);

    tmp1 = mw_neg(&M);
    tmp1 = mw_div(&tmp1, r);
    tmp2 = mw_add(&a, &hypot);
    tmp3 = sqr(&hypot);
    tmp4 = mw_mul(&a, &hypot);
    tmp3 = mw_add(&tmp3, &tmp4);
    tmp2 = mw_div(&tmp2, &tmp3);
    const real c = mw_mul(&tmp1, &tmp2);

    tmp1 = mw_div(&c, r);

    acc.x = mw_mul(&pos->x, &tmp1);
    acc.y = mw_mul(&pos->y, &tmp1);
    acc.z = mw_mul(&pos->z, &tmp1);

    return acc;

}

static inline mwvector NFWMHaloAccel(const Halo* h, mwvector* pos, real* r)
{
    mwvector acc;
    real tmp1, tmp2, tmp3;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real ar = mw_add(&a, r);

    tmp1 = mw_neg(&M);
    tmp2 = sqr(r);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_div(&ar, &a);
    tmp2 = mw_log(&tmp2);
    tmp2 = mw_div(&tmp2, r);
    tmp3 = inv(&ar);
    tmp2 = mw_sub(&tmp2, &tmp3);
    const real c = mw_mul(&tmp1, &tmp2);

    acc.x = mw_mul(&pos->x, &c);
    acc.y = mw_mul(&pos->y, &c);
    acc.z = mw_mul(&pos->z, &c);

    return acc;

}

static inline mwvector plummerHaloAccel(const Halo* h, mwvector* pos, real* r)
{
    mwvector acc;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    real tmp = mw_hypot(&a,r);

    tmp = cube(&tmp);
    tmp = mw_div(&M, &tmp);
    tmp = mw_neg(&tmp);

    acc.x = mw_mul(&pos->x, &tmp);
    acc.y = mw_mul(&pos->y, &tmp);
    acc.z = mw_mul(&pos->z, &tmp);

    return acc;
}

static inline mwvector hernquistHaloAccel(const Halo* h, mwvector* pos, real* r)
{
    mwvector acc;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    real tmp = mw_add(&a, r);

    tmp = sqr(&tmp);
    tmp = mw_mul(r, &tmp);
    tmp = mw_div(&M, &tmp);
    tmp = mw_neg(&tmp);

    acc.x = mw_mul(&pos->x, &tmp);
    acc.y = mw_mul(&pos->y, &tmp);
    acc.z = mw_mul(&pos->z, &tmp);

    return acc;
}

static inline mwvector ninkovicHaloAccel(const Halo* h, mwvector* pos, real* r)      /*Special case of Ninkovic Halo (l1=0,l2=3,l3=2) (Ninkovic 2017)*/
{
    mwvector acc;
    real tmp1, tmp2, tmp3;
    const real rho0 = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real_0 lambda = h->lambda;

    const real z  = mw_div(r, &a);

    tmp1 = inv(&a);
    const real zl = mw_mul_s(&tmp1, lambda);

    tmp1 = cube(&a);
    tmp1 = mw_mul(&rho0, &tmp1);
    const real f  = mw_mul_s(&tmp1, 4.0*M_PI/3.0);

    real mass_enc;

    if (showRealValue(&r) > lambda)
    {
        tmp1 = cube(&zl);
        tmp1 = mw_log1p(&tmp1);
        tmp2 = cube(&zl);
        tmp3 = mw_add_s(&tmp2, 1.0);
        tmp2 = mw_div(&tmp2, &tmp3);
        tmp1 = mw_sub(&tmp1, &tmp2);
        mass_enc = mw_mul(&f, &tmp1);
    }
    else
    {
        tmp1 = cube(&z);
        tmp1 = mw_log1p(&tmp1);
        tmp2 = cube(&z);
        tmp3 = cube(&zl);
        tmp3 = mw_add_s(&tmp3, 1.0);
        tmp2 = mw_div(&tmp2, &tmp3);
        tmp1 = mw_sub(&tmp1, &tmp2);
        mass_enc = mw_mul(&f, &tmp1);
    }

    tmp1 = mw_neg(&mass_enc);
    tmp2 = cube(r);
    tmp1 = mw_div(&tmp1, &tmp2);

    acc.x = mw_mul(&pos->x, &tmp1);
    acc.y = mw_mul(&pos->y, &tmp1);
    acc.z = mw_mul(&pos->z, &tmp1);

    return acc;
}

mwvector nbExtAcceleration(const Potential* pot, mwvector* pos, real_0 time)
{
    //mw_printf("POS = [%.15f,%.15f,%.15f]\n", showRealValue(&X(pos)), showRealValue(&Y(pos)), showRealValue(&Z(pos)));
    mwvector acc, acctmp;
    real_0 limit_val = mw_pow_0(2.0,-8.0);

    /* Change r if less than limit. Done this way to pipeline this step*/
    real r = mw_absv(pos);
    real limit = r;
    setRealValue(&limit, limit_val);

    real check1 = mw_mul_s(&r, (showRealValue(&r)>limit_val));
    real check2 = mw_mul_s(&limit, (showRealValue(&r)<=limit_val));
    r = mw_add(&check1, &check2);

    /*Calculate the Disk Accelerations*/
    switch (pot->disk.type)
    {
        case FreemanDisk:
            acc = freemanDiskAccel(&pot->disk, pos);
            break;
        case MiyamotoNagaiDisk:
            acc = miyamotoNagaiDiskAccel(&pot->disk, pos);
            break;
        case DoubleExponentialDisk:
            acc = doubleExponentialDiskAccel(&pot->disk, pos);
            break;
        case Sech2ExponentialDisk:
            acc = sech2ExponentialDiskAccel(&pot->disk, pos);
            break;
        case NoDisk:
            X(&acc) = ZERO_REAL;
            Y(&acc) = ZERO_REAL;
            Z(&acc) = ZERO_REAL;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid primary disk type in external acceleration\n");
    }
    if( isnan(showRealValue(&X(&acc))) || isnan(showRealValue(&Y(&acc))) || isnan(showRealValue(&Y(&acc))) )
    {
        mw_printf("BAD DISK: %s\n", showDiskT(pot->disk.type));
    }
    //mw_printf("Disk  Acceleration = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acc)), showRealValue(&Y(&acc)), showRealValue(&Z(&acc)));

    /*Calculate Second Disk Accelerations*/
    switch (pot->disk2.type)
    {
        case FreemanDisk:
            acctmp = freemanDiskAccel(&pot->disk2, pos);
            break;
        case MiyamotoNagaiDisk:
            acctmp = miyamotoNagaiDiskAccel(&pot->disk2, pos);
            break;
        case DoubleExponentialDisk:
            acctmp = doubleExponentialDiskAccel(&pot->disk2, pos);
            break;
        case Sech2ExponentialDisk:
            acctmp = sech2ExponentialDiskAccel(&pot->disk2, pos);
            break;
        case OrbitingBar:
            acctmp = orbitingBarAccel(&pot->disk2, pos, time);
            break;
        case NoDisk:
            X(&acctmp) = ZERO_REAL;
            Y(&acctmp) = ZERO_REAL;
            Z(&acctmp) = ZERO_REAL;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid secondary disk type in external acceleration\n");
    }
    if( isnan(showRealValue(&X(&acctmp))) || isnan(showRealValue(&Y(&acctmp))) || isnan(showRealValue(&Y(&acctmp))) )
    {
        mw_printf("BAD DISK2: %s\n", showDiskT(pot->disk2.type));
    }
    acc = mw_addv(&acc, &acctmp);
    //mw_printf("Disk2 Acceleration = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acctmp)), showRealValue(&Y(&acctmp)), showRealValue(&Z(&acctmp)));
    //mw_printf("Total Acceleration = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acc)), showRealValue(&Y(&acc)), showRealValue(&Z(&acc)));

    /*Calculate the Halo Accelerations*/
    switch (pot->halo.type)
    {
        case LogarithmicHalo:
            acctmp = logHaloAccel(&pot->halo, pos);
            break;
        case NFWHalo:
            acctmp = nfwHaloAccel(&pot->halo, pos, &r);
            break;
        case TriaxialHalo:
            acctmp = triaxialHaloAccel(&pot->halo, pos);
            break;
        case CausticHalo:
            acctmp = causticHaloAccel(&pot->halo, pos, &r);
            break;
        case AllenSantillanHalo:
            acctmp = ASHaloAccel(&pot->halo, pos, &r);
            break;
        case WilkinsonEvansHalo:
            acctmp = WEHaloAccel(&pot->halo, pos, &r);
            break;
        case NFWMassHalo:
            acctmp = NFWMHaloAccel(&pot->halo, pos, &r);
            break;
        case PlummerHalo:
            acctmp = plummerHaloAccel(&pot->halo, pos, &r);
            break;
        case HernquistHalo:
            acctmp = hernquistHaloAccel(&pot->halo, pos, &r);
            break;
        case NinkovicHalo:
            acctmp = ninkovicHaloAccel(&pot->halo, pos, &r);
            break;
        case NoHalo:
            X(&acctmp) = ZERO_REAL;
            Y(&acctmp) = ZERO_REAL;
            Z(&acctmp) = ZERO_REAL;
            break;
        case InvalidHalo:
        default:
            mw_fail("Invalid halo type in external acceleration\n");
    }
    if( isnan(showRealValue(&X(&acctmp))) || isnan(showRealValue(&Y(&acctmp))) || isnan(showRealValue(&Y(&acctmp))) )
    {
        mw_printf("BAD HALO: %s\n", showHaloT(pot->halo.type));
    }
    acc = mw_addv(&acc, &acctmp);
    //mw_printf("Halo  Acceleration = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acctmp)), showRealValue(&Y(&acctmp)), showRealValue(&Z(&acctmp)));
    //mw_printf("Total Acceleration = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acc)), showRealValue(&Y(&acc)), showRealValue(&Z(&acc)));

    /*Calculate the Bulge Accelerations*/
    switch (pot->sphere[0].type)
    {
        case HernquistSpherical:
            acctmp = hernquistSphericalAccel(&pot->sphere[0], pos, &r);
            break;
        case PlummerSpherical:
            acctmp = plummerSphericalAccel(&pot->sphere[0], pos, &r);
            break;
        case NoSpherical:
            X(&acctmp) = ZERO_REAL;
            Y(&acctmp) = ZERO_REAL;
            Z(&acctmp) = ZERO_REAL;
            break;
        case InvalidSpherical:
        default:
            mw_fail("Invalid bulge type in external acceleration\n");
    }
    if( isnan(showRealValue(&X(&acctmp))) || isnan(showRealValue(&Y(&acctmp))) || isnan(showRealValue(&Y(&acctmp))) )
    {
        mw_printf("BAD BULGE: %s\n", showSphericalT(pot->sphere[0].type));
    }
    acc = mw_addv(&acc, &acctmp);
    //mw_printf("Bulge Acceleration = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acctmp)), showRealValue(&Y(&acctmp)), showRealValue(&Z(&acctmp)));
    //mw_printf("Total Acceleration = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acc)), showRealValue(&Y(&acc)), showRealValue(&Z(&acc)));

    /*For debugging acceleration values*/
//    if (!isfinite(mw_absv(acc)))
//    {
//        mw_printf("ERROR: UNNATURAL ACCELERATION CALCULATED!\n");
//        mw_printf("POS = [%.15f,%.15f,%.15f]\n", showRealValue(&X(pos)), showRealValue(&Y(pos)), showRealValue(&Z(pos)));
//        mw_printf("ACC = [%.15f,%.15f,%.15f]\n", showRealValue(&X(&acc)), showRealValue(&Y(&acc)), showRealValue(&Z(&acc)));
//    }

    return acc;
}







