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

static inline real leg_pol(real x, int l)
{
    real sum = ZERO_REAL;
    for (int m = 0; m < mw_floor_0(l/2)+1; m++)
    {
        sum = mw_add(sum, mw_mul_s(mw_pow(x, mw_real_const(l-2*m)), mw_pow_0(-1, m)/mw_pow_0(2,l)*binom(l,m)*binom(2*(l-m),l)));
    }
    /*mw_printf("P(%.15f, %u) = %.15f \n",x,l,sum);*/
    return sum;
}

static inline real_0 leg_pol_zeros(int l)
{
    if (l % 2 == 0)
    {
        return mw_pow_0(-1,l/2)*mw_exp_0(lnfact(l)-2*lnfact(l/2) - l*mw_log(2));
    }
    else
    {
        return 0.0;
    }
}

static inline real leg_pol_derv(real x, int l)
{
    real sum = ZERO_REAL;
    for (int m = 0; m < mw_floor((l-1)/2)+1; m++)
    {
        sum = mw_add(sum, mw_mul_s(mw_pow(x, l - 2*m - 1), mw_pow_0(-1, m)/mw_pow_0(2,l)*binom(l,m)*binom(2*(l-m),l)*(l - 2*m)));
    }
    /*mw_printf("P'(%.15f, %u) = %.15f \n",x,l,sum);*/
    return sum;
}

static inline real lower_gamma(int n, real x)
{
    real sum = ZERO_REAL;
    for (int k = 0; k < n; k++)
    {
       sum = mw_add(sum, mw_mul_s(mw_pow(x,k), inv_0(mw_exp_0(lnfact(k)))));
    }
    /*mw_printf("g(%u, %.15f) = %.15f \n",n,x,mw_exp(lnfact(n-1))*(1-mw_exp(-x)*sum));*/
    return mw_exp(lnfact(n-1))*(1-mw_exp(-x)*sum);
}

static inline real G_cont_frac(real x, int n, int k)
{
    if (k > 15)
    {
        return mw_add(x, mw_real_const(n + k));
    }
    else
    {
        return mw_add(x, mw_mul_s(inv(mw_add(mw_real_const(1.0), mw_mul_s(inv(G_cont_frac(x,n,k+1)), (1+k)))), (n+k)));
    }
}

static inline real GenExpIntegral(int n, real x)
{
    return inv(mw_mul(G_cont_frac(x,n,0), mw_exp(x)));
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

static inline real RSechIntegrand (real k, real R, real Rd, real z, real zd)
{
    real ExpStuff = mw_sub(mw_mul(aExp(k,R,Rd),besselK1(mw_mul(k,R))), mw_mul(bExp(k,R,Rd),besselI1(mw_mul(k,R))));
    real TrigStuff = mw_div(mw_cos(mw_mul(k,z)), mw_sinh(mw_mul_s(mw_mul(k,zd),M_PI/2.0)));
    real val = mw_mul(mw_mul(mw_mul(sqr(k),zd),TrigStuff),ExpStuff);
    //mw_printf("RSech(%.15f,%.15f,%.15f,%.15f,%.15f) = %.15f\n",k,R,Rd,z,zd,val);
    return val;
}

static inline real ZSechIntegrand (real k, real R, real Rd, real z, real zd)
{
    real ExpStuff = mw_sub(mw_mul(aExp(k,R,Rd),besselK0(mw_mul(k,R))), mw_mul(bExp(k,R,Rd),besselI0(mw_mul(k,R))));
    real TrigStuff = mw_div(mw_sin(mw_mul(k,z)), mw_sinh(mw_mul_s(mw_mul(k,zd),M_PI/2.0)));
    real val = mw_mul(mw_mul(mw_mul(sqr(k),zd),TrigStuff),ExpStuff);
    //mw_printf("ZSech = %.15f\n",val);
    return val;
}

/**********************************************************************************************************************************************************************************************/

mwvector pointAccel(const mwvector pos, const mwvector pos1, const real mass)
{
    mwvector v = mw_subv(pos1, pos);
    real dist = mw_distv(pos, pos1);
    real tmp = mw_div(mass, cube(dist));
    mw_incmulvs(v, tmp);
    return v;
}

mwvector plummerAccel(const mwvector pos, const mwvector pos1, const real mass, const real scale)
{
    mwvector v = mw_subv(pos1, pos);
    real dist = mw_distv(pos, pos1);
    real tmp = mw_hypot(scale, dist);
    mw_incmulvs(v, mw_div(mass, cube(tmp)));
    return v;
}


/*spherical bulge potentials*/

static inline mwvector hernquistSphericalAccel(const Spherical* sph, mwvector pos, real r)
{
    const real M = mw_real_var(sph->scale, 11);
    const real a = mw_real_var(sph->scale, 12);
    const real tmp = mw_add(a, r);

    const real factor = mw_neg(mw_div(M, mw_mul(r, sqr(tmp))));

    return mw_mulvs(pos, factor);
}

static inline mwvector plummerSphericalAccel(const Spherical* sph, mwvector pos, real r)
{
    const real M = mw_real_var(sph->scale, 11);
    const real a = mw_real_var(sph->scale, 12);
    const real tmp = mw_hypot(a,r);

    const real factor = mw_neg(mw_div(M, cube(tmp)));

    return mw_mulvs(pos, factor);
}

/* Disk Potentials */

static inline mwvector miyamotoNagaiDiskAccel(const Disk* disk, mwvector pos, real r)
{
    mwvector acc;
    const real M   = mw_real_var(disk->mass, 13);
    const real a   = mw_real_var(disk->scaleLength, 14);
    const real b   = mw_real_var(disk->scaleHeight, 15);
    const real zp  = mw_hypot(Z(pos), b);
    const real azp = mw_add(a, zp);

    const real rp  = mw_hypot(mw_hypot(X(pos), Y(pos)), azp);
    const real rth = cube(rp);  /* rp ^ (3) */

    X(acc) = mw_mul_s(mw_div(mw_mul(M, X(pos)), rth), -1.0);
    Y(acc) = mw_mul_s(mw_div(mw_mul(M, Y(pos)), rth), -1.0);
    Z(acc) = mw_mul_s(mw_mul(M, mw_mul(Z(pos), mw_div(azp, mw_mul(zp, rth)))), -1.0);

    //mw_printf("Acceleration[AX,AY,AZ] = [%.15f,%.15f,%.15f]\n",X(acc),Y(acc),Z(acc));
    
    /*mw_printf("disk x acc: %.15f\n", X(acc));
    mw_printf("disk y acc: %.15f\n", Y(acc));
    mw_printf("disk z acc: %.15f\n", Z(acc));*/
    return acc;
}

/*WARNING: This potential is currently not in use as it is non-physical and difficult to test.*/
static inline mwvector freemanDiskAccel(const Disk* disk, mwvector pos, real r)     /*From Freeman 1970*/
{
    mwvector acc;
    /*Reset acc vector*/
    X(acc) = ZERO_REAL;
    Y(acc) = ZERO_REAL;
    Z(acc) = ZERO_REAL;

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
static inline mwvector doubleExponentialDiskAccel(const Disk* disk, mwvector pos, real r)
{
    //mw_printf("Calculating Acceleration\n");
    //mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    mwvector acc;

    const real R = mw_hypot(X(pos), Y(pos));
    mwvector R_hat;
    X(R_hat) = mw_div(X(pos), R);
    Y(R_hat) = mw_div(Y(pos), R);
    Z(R_hat) = ZERO_REAL;

    mwvector Z_hat;
    X(Z_hat) = ZERO_REAL;
    Y(Z_hat) = ZERO_REAL;
    Z(Z_hat) = mw_div(Z(pos), mw_abs(Z(pos)));

    const real M    = mw_real_var(disk->mass, 13);
    const real Rd   = mw_real_var(disk->scaleLength, 14);
    const real zd   = mw_real_var(disk->scaleHeight, 15);
    const real z    = Z(pos);
    const real_0 h  = 0.001;

    const real a = inv(Rd);
    const real b = inv(zd);
    
    real R_piece = ZERO_REAL;
    real z_piece = ZERO_REAL;

    real_0 j0_zero;
    real_0 psi_in_0;
    real_0 psi_0;
    real_0 psi_prime_0;
    real_0 j0_x;
    real j0_w;
    real fun_0;

    real_0 j1_zero;
    real_0 psi_in_1;
    real_0 psi_1;
    real_0 psi_prime_1;
    real_0 j1_x;
    real j1_w;
    real fun_1;

    for (int n = 0; n < 150; n+=1)    // Hidenori Ogata's Numerical Integration Formula Based on the Bessel Functions
    {
        j0_zero = besselJ0_zero(n)/M_PI;
        j1_zero = besselJ1_zero(n)/M_PI;
        psi_in_0 = h * j0_zero;
        psi_in_1 = h * j1_zero;
        psi_0 = psi_in_0 * mw_sinh_0(M_PI / 2.0 * mw_sinh_0(psi_in_0)) / mw_cosh_0(M_PI / 2.0 * mw_sinh_0(psi_in_0));
        psi_1 = psi_in_1 * mw_sinh_0(M_PI / 2.0 * mw_sinh_0(psi_in_1)) / mw_cosh_0(M_PI / 2.0 * mw_sinh_0(psi_in_1));
        psi_prime_0 = (mw_sinh_0(M_PI * mw_sinh_0(psi_in_0)) + M_PI * psi_in_0 * mw_cosh_0(psi_in_0)) / (mw_cosh_0(M_PI * mw_sinh(psi_in_0)) + 1.0);
        psi_prime_1 = (mw_sinh_0(M_PI * mw_sinh_0(psi_in_1)) + M_PI * psi_in_1 * mw_cosh_0(psi_in_1)) / (mw_cosh_0(M_PI * mw_sinh(psi_in_1)) + 1.0);
        j0_x = M_PI / h * psi_0;
        j1_x = M_PI / h * psi_1;

        j0_w = mw_mul_s(mw_div(besselJ0(mw_real_const(j0_x)), sqr(besselJ1(mw_real_const(M_PI * j0_zero)))), 2.0 / (M_PI * j0_zero) * psi_prime_0);
        j1_w = mw_mul_s(mw_div(besselJ1(mw_real_const(j1_x)), sqr(besselJ2(mw_real_const(M_PI * j1_zero)))), 2.0 / (M_PI * j1_zero) * psi_prime_1);

        real j0_x_R = mw_mul_s(inv(R), j0_x);
        real j1_x_R = mw_mul_s(inv(R), j1_x);

        fun_0 = mw_mul(mw_mul(inv(cube(mw_hypot(a,j0_x_R))), mw_div(mw_sub(inv(mw_exp(mw_mul(j0_x_R, mw_abs(z)))), inv(mw_exp(mw_mul(b, mw_abs(z))))), mw_sub(sqr(b), sqr(j0_x_R)))), j0_x_R);
        fun_1 = mw_mul_s(mw_mul(inv(cube(mw_hypot(a,j1_x_R))), mw_div(mw_sub(mw_mul(b, inv(mw_exp(j1_x_R * mw_abs(z)))), mw_mul(j1_x_R, inv(mw_exp(b * mw_abs(z))))), mw_sub(sqr(b), sqr(j1_x_R)))), j1_x);

        real z_pieceAdd = mw_mul_s(mw_mul(mw_mul(mw_div(mw_mul(a, b), R), fun_0), j0_w), 4.0 * M_PI);
        real R_pieceAdd = mw_mul_s(mw_mul(mw_mul(mw_div(a, sqr(R)), fun_1), j1_w), 4.0 * M_PI);
        
        z_piece = mw_add(z_piece, z_pieceAdd);
        R_piece = mw_add(R_piece, R_pieceAdd);
    }

    mwvector R_comp = mw_mulvs(R_hat, mw_mul_s(mw_mul(mw_mul(mw_mul(M, sqr(a)), b), R_piece), -1.0 / 4.0 / M_PI));
    mwvector Z_comp = mw_mulvs(Z_hat, mw_mul_s(mw_mul(mw_mul(mw_mul(M, sqr(a)), b), z_piece), -1.0 / 4.0 / M_PI));

    X(acc) = mw_add(X(R_comp), X(Z_comp));
    Y(acc) = mw_add(Y(R_comp), Y(Z_comp));
    Z(acc) = mw_add(Z(R_comp), Z(Z_comp));

    //real magnitude = mw_hypot(mw_hypot(X(acc), Y(acc)), Z(acc));
    //mw_printf("Acceleration[AX,AY,AZ] = [%.15f,%.15f,%.15f]   Magnitude = %.15f\n",X(acc),Y(acc),Z(acc),magnitude);
    return acc;
}

/*WARNING: This potential can take a while to integrate if any part of the orbit extends past 100 times the scaleLength*/
static inline mwvector sech2ExponentialDiskAccel(const Disk* disk, mwvector pos, real r)
{
    //mw_printf("Calculating Acceleration\n");
    //mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    mwvector acc;

    const real R = mw_hypot(X(pos), Y(pos));
    mwvector R_hat;
    X(R_hat) = mw_div(X(pos), R);
    Y(R_hat) = mw_div(Y(pos), R);
    Z(R_hat) = ZERO_REAL;

    mwvector Z_hat;
    X(Z_hat) = ZERO_REAL;
    Y(Z_hat) = ZERO_REAL;
    Z(Z_hat) = mw_real_const(1.0);

    const real M    = mw_real_var(disk->mass, 13);
    const real Rd   = mw_real_var(disk->scaleLength, 14);
    const real zd   = mw_real_var(disk->scaleHeight, 15);
    const real z = Z(pos);

    const int n = 15;
    const real_0 a = 0.0;
    const real_0 b = 60.0 / showRealValue(R); // Should be infinity, but this should be the point where the integral becomes negligible
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
            Rpiece = mw_add(Rpiece, mw_mul_s(RSechIntegrand(k_val,R,Rd,z,zd), h*weight[j]/2.0));
            Zpiece = mw_add(Zpiece, mw_mul_s(ZSechIntegrand(k_val,R,Rd,z,zd), h*weight[j]/2.0));
        }
        integralR  = mw_add(integralR, Rpiece);
        integralZ  = mw_add(integralZ, Zpiece);
    }

    mwvector R_comp = mw_mulvs(R_hat, mw_mul_s(mw_mul(M,integralR), -1.0));
    mwvector Z_comp = mw_mulvs(Z_hat, mw_mul_s(mw_mul(M,integralZ), -1.0));

    X(acc) = mw_add(X(R_comp), X(Z_comp));
    Y(acc) = mw_add(Y(R_comp), Y(Z_comp));
    Z(acc) = mw_add(Z(R_comp), Z(Z_comp));

    //real magnitude = mw_sqrt(sqr(X(acc))+sqr(Y(acc))+sqr(Z(acc)));

    //mw_printf("Acceleration[AX,AY,AZ] = [%.15f,%.15f,%.15f]   Magnitude = %.15f\n",X(acc),Y(acc),Z(acc),magnitude);

    return acc;
}

//Softened needle bar potential
static inline mwvector orbitingBarAccel(const Disk* disk, mwvector pos, real r, real_0 time)
{
    real amp = mw_real_var(disk->mass, 13);
    real a = mw_real_var(disk->scaleLength, 14);
    real b = mw_real_const(1.4);                 //Triaxial softening length
    real c = mw_real_const(1.0);                 //Prolate softening length
    
    real_0 curAngle = (disk->patternSpeed * time * -1)+disk->startAngle;
    //first rotate pos curAngle * -1 radians to emulate the current angle of the bar
    real Radi = mw_hypot(X(pos), Y(pos));
    real Phi = mw_atan(mw_div(Y(pos), X(pos)));
    Phi = mw_sub(Phi, mw_real_const(curAngle));
    if(showRealValue(X(pos)) < 0){
        Radi = mw_mul_s(Radi, -1.0);
    }
    real x = mw_mul(Radi, mw_cos(Phi));
    real y = mw_mul(Radi, mw_sin(Phi)); 
    real z = Z(pos);

    //calculate force in accordance with the galpy implementation
    real secondpart = mw_hypot(y, mw_add(b, mw_hypot(z,c)));
    real Tp = mw_hypot(mw_add(a,x), secondpart);
    real Tm = mw_hypot(mw_sub(a,x), secondpart);

    real thirdpart = mw_sub(mw_add(Tp, Tm), mw_mul_s(mw_div(sqr(x), mw_add(Tp, Tm)),4.0));

    mwvector force;
    force.x = mw_mul_s(mw_div(x, mw_mul(mw_mul(Tp, Tm), mw_add(Tp, Tm))), -2.0);
    force.y = mw_mul_s(mw_div(mw_mul(mw_div(mw_div(y, Tp), Tm),thirdpart),sqr(secondpart)), -0.5);
    force.z = mw_div(mw_mul(mw_div(mw_mul(force.y, z), y), mw_add(b, mw_hypot(z,c))), mw_hypot(z,c));
    
    //undo the pos rotation and calculate acceleration from the force vector we got
    mwvector acc;
    real_0 cp = mw_cos_0(curAngle);
    real_0 sp = mw_sin_0(curAngle);
    acc.x = mw_sub(mw_mul_s(force.x, cp), mw_mul_s(force.y, sp));
    acc.y = mw_add(mw_mul_s(force.x, sp), mw_mul_s(force.y, cp));
    acc.z = force.z;
    
    acc.x = mw_mul(acc.x, amp);
    acc.y = mw_mul(acc.y, amp);
    acc.z = mw_mul(acc.z, amp);

    return acc;
}


/*If you want to test the time dependency of the bar with the bar as a point mass, comment out 
the above function and uncomment the one below*/

/*
static inline mwvector orbitingBarAccel(const Disk* disk, mwvector pos, real r, real time)
{
    //mw_printf("Calculating Acceleration\n");
    //mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    //mw_printf("r = %.15f\n", r);

    mwvector pointPos;
    pointPos.z = 0;
    real curAngle = (disk->patternSpeed * time * -1)+disk->startAngle;
    curAngle = curAngle - M_PI;//this is because the sun is negative in our coordinate system
    pointPos.x = cos (curAngle) * disk->scaleLength; //this is assuming top-down
    pointPos.y = sin (curAngle) * disk->scaleLength;

    real dist = mw_distv(pos, pointPos);

    mwvector acc = mw_divvs(mw_subv(pos, pointPos), dist);//get direction from pos to pointPos
    real totalAcc = disk->mass/(dist*dist);//a = Gm/r^2
    acc = mw_mulvs(acc, totalAcc);

    //mw_printf("curAngle: %.15f\n", curAngle);
    //mw_printf("pointPos: [%.15f,%.15f,%.15f]\n", X(pointPos), Y(pointPos), Z(pointPos));
    //mw_printf("Accel: [%.15f,%.15f,%.15f]\n", X(acc), Y(acc), Z(acc));
    //mw_printf("point x acc: %.15f\n", X(acc));
    //mw_printf("point y acc: %.15f\n", Y(acc));
    //mw_printf("point z acc: %.15f\n", Z(acc));
    return acc;
}
*/

static inline mwvector logHaloAccel(const Halo* halo, mwvector pos)
{
    mwvector acc;

    const real v0 = mw_real_var(halo->vhalo, 16);
    const real d  = mw_real_var(halo->scaleLength, 17);
    const real q  = mw_real_var(halo->flattenZ, 18);

    const real denom = mw_add(mw_add(mw_add(sqr(d), sqr(X(pos))), sqr(Y(pos))), sqr(mw_div(Z(pos), q)));
    const real k = mw_mul_s(mw_div(sqr(v0), denom), -2.0);

    X(acc) = mw_mul(k, X(pos));
    Y(acc) = mw_mul(k, Y(pos));
    Z(acc) = mw_mul(k, mw_div(Z(pos), sqr(q)));

    return acc;
}

static inline mwvector nfwHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real v0 = mw_real_var(h->vhalo, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real M = mw_mul_s(mw_mul(sqr(v0), a), inv_0(0.2162165954)); /*Maximum of [ln(1+x)/x - 1/(1+x)]*/
    const real ar = mw_add(a, r);

    const real c = mw_mul(mw_div(mw_mul_s(M,-1.0), sqr(r)), mw_sub(mw_div(mw_log(mw_div(ar, a)), r), inv(ar)));

    return mw_mulvs(pos, c);
}

/* CHECKME: Seems to have precision related issues for a small number of cases for very small qy */
static inline mwvector triaxialHaloAccel(const Halo* h, mwvector pos, real r)  /** Triaxial Logarithmic **/
{
    mwvector acc;

    /* TODO: More things here can be cached */
    const real qzs      = sqr(mw_real_var(h->flattenZ, 18));
    const real rhalosqr = sqr(mw_real_var(h->scaleLength, 17));
    const real mvsqr    = mw_mul_s(sqr(mw_real_var(h->vhalo, 16)), -1.0);

    const real xsqr = sqr(X(pos));
    const real ysqr = sqr(Y(pos));
    const real zsqr = sqr(Z(pos));

    const real arst  = rhalosqr + mw_mul_s(xsqr, h->c1) + mw_mul_s(mw_mul(X(pos), Y(pos)), h->c3) + mw_mul_s(ysqr, h->c2);
    const real arst2 = mw_add(mw_div(zsqr, qzs), arst);

    X(acc) = mw_mul(mvsqr, mw_div(mw_add(mw_mul_s(X(pos), (2.0 * h->c1)), mw_mul_s(Y(pos), h->c3)), arst2));

    Y(acc) = mw_mul(mvsqr, mw_div(mw_add(mw_mul_s(Y(pos), (2.0 * h->c2)), mw_mul_s(X(pos), h->c3)), arst2));

    Z(acc) = mw_div(mw_mul_s(mw_mul(mvsqr, Z(pos)), 2.0), mw_add(mw_mul(qzs, arst), zsqr));

    return acc;
}

static inline mwvector ASHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real_0 gam = h->gamma;
    const real_0 lam = h->lambda;
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real scaleR = mw_div(r, a);
    const real scaleL = mw_mul_s(inv(a), lam);
    real factor;
    real c;

    if (r<lam)
    {
        factor = mw_neg(mw_div(M, mw_mul(a, r)));
        c = mw_mul(factor, mw_div(mw_pow(scaleR, mw_real_const(gam-1.0)), mw_add(mw_real_const(1.0), mw_pow(scaleR, mw_real_const(gam-1.0)))));
        //c = -(M/(a*r))*mw_pow(scaleR,gam-1.0)/(1.0+mw_pow(scaleR,gam-1.0));
    }
    else
    {
        factor = mw_neg(mw_div(M, sqr(r)));
        c = mw_mul(factor, mw_div(mw_pow(scaleL, mw_real_const(gam)), mw_add(mw_real_const(1.0), mw_pow(scaleL, mw_real_const(gam-1.0)))));
        //c = -(M/sqr(r))*mw_pow(scaleL,gam)/(1.0+mw_pow(scaleL,gam-1.0));
    }

    return mw_mulvs(pos, mw_div(c, r));
}

static inline mwvector WEHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real sum2 = mw_add(sqr(a), sqr(r));

    const real c = mw_mul(mw_div(mw_neg(M), r), mw_div(mw_add(a, mw_sqrt(sum2)), mw_add(sum2, mw_mul(a, mw_sqrt(sum2)))));

    return mw_mulvs(pos, mw_div(c, r));

}

static inline mwvector NFWMHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real ar = mw_add(a, r);

    const real c = mw_mul(mw_div(mw_mul_s(M,-1.0), sqr(r)), mw_sub(mw_div(mw_log(mw_div(ar, a)), r), inv(ar)));

    return mw_mulvs(pos, c);

}

static inline mwvector plummerHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real tmp = mw_hypot(a,r);

    const real factor = mw_neg(mw_div(M, cube(tmp)));

    return mw_mulvs(pos, factor);
}

static inline mwvector hernquistHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real tmp = mw_add(a, r);

    const real factor = mw_neg(mw_div(M, mw_mul(r, sqr(tmp))));

    return mw_mulvs(pos, factor);
}

static inline mwvector ninkovicHaloAccel(const Halo* h, mwvector pos, real r)      /*Special case of Ninkovic Halo (l1=0,l2=3,l3=2) (Ninkovic 2017)*/
{
    const real rho0 = mw_real_var(h->rho0, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real_0 lambda = h->lambda;

    const real z  = mw_div(r, a);
    const real zl = mw_mul_s(inv(a), lambda);
    const real f  = mw_mul_s(mw_mul(rho0, cube(a)), 4.0*M_PI/3.0);

    real mass_enc;

    if (r > lambda)
    {
        mass_enc = mw_mul(f, mw_sub(mw_log1p(cube(zl)), mw_div(cube(zl), mw_add(mw_real_const(1.0), cube(zl)))));
    }
    else
    {
        mass_enc = mw_mul(f, mw_sub(mw_log1p(cube(z)), mw_div(cube(z), mw_add(mw_real_const(1.0), cube(zl)))));
    }

    mwvector acc = mw_mulvs(pos, mw_div(mw_neg(mass_enc), cube(r)));

    return acc;
}

mwvector nbExtAcceleration(const Potential* pot, mwvector pos, real_0 time)
{
    mwvector acc, acctmp;
    real_0 limit_val = mw_pow(2.0,-8.0);

    /* Change r if less than limit. Done this way to pipeline this step*/
    real r = mw_absv(pos);
    real limit = r;
    setRealValue(&limit, limit_val);

    r = mw_add(mw_mul_s(limit, (showRealValue(r) <= limit_val)), mw_mul(r, (showRealValue(r) > limit)));

    /*Calculate the Disk Accelerations*/
    switch (pot->disk.type)
    {
        case FreemanDisk:
            acc = freemanDiskAccel(&pot->disk, pos, r);
            break;
        case MiyamotoNagaiDisk:
            acc = miyamotoNagaiDiskAccel(&pot->disk, pos, r);
            break;
        case DoubleExponentialDisk:
            acc = doubleExponentialDiskAccel(&pot->disk, pos, r);
            break;
        case Sech2ExponentialDisk:
            acc = sech2ExponentialDiskAccel(&pot->disk, pos, r);
            break;
        case NoDisk:
            X(acc) = ZERO_REAL;
            Y(acc) = ZERO_REAL;
            Z(acc) = ZERO_REAL;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid primary disk type in external acceleration\n");
    }

    /*Calculate Second Disk Accelerations*/
    switch (pot->disk2.type)
    {
        case FreemanDisk:
            acctmp = freemanDiskAccel(&pot->disk2, pos, r);
            break;
        case MiyamotoNagaiDisk:
            acctmp = miyamotoNagaiDiskAccel(&pot->disk2, pos, r);
            break;
        case DoubleExponentialDisk:
            acctmp = doubleExponentialDiskAccel(&pot->disk2, pos, r);
            break;
        case Sech2ExponentialDisk:
            acctmp = sech2ExponentialDiskAccel(&pot->disk2, pos, r);
            break;
        case OrbitingBar:
            acctmp = orbitingBarAccel(&pot->disk2, pos, r, time);
            break;
        case NoDisk:
            X(acctmp) = ZERO_REAL;
            Y(acctmp) = ZERO_REAL;
            Z(acctmp) = ZERO_REAL;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid secondary disk type in external acceleration\n");
    }
    mw_incaddv(acc, acctmp);

    /*Calculate the Halo Accelerations*/
    switch (pot->halo.type)
    {
        case LogarithmicHalo:
            acctmp = logHaloAccel(&pot->halo, pos);
            break;
        case NFWHalo:
            acctmp = nfwHaloAccel(&pot->halo, pos, r);
            break;
        case TriaxialHalo:
            acctmp = triaxialHaloAccel(&pot->halo, pos, r);
            break;
        case CausticHalo:
            acctmp = causticHaloAccel(&pot->halo, pos, r);
            break;
        case AllenSantillanHalo:
            acctmp = ASHaloAccel(&pot->halo, pos, r);
            break;
        case WilkinsonEvansHalo:
            acctmp = WEHaloAccel(&pot->halo, pos, r);
            break;
        case NFWMassHalo:
            acctmp = NFWMHaloAccel(&pot->halo, pos, r);
            break;
        case PlummerHalo:
            acctmp = plummerHaloAccel(&pot->halo, pos, r);
            break;
        case HernquistHalo:
            acctmp = hernquistHaloAccel(&pot->halo, pos, r);
            break;
        case NinkovicHalo:
            acctmp = ninkovicHaloAccel(&pot->halo, pos, r);
            break;
        case NoHalo:
            X(acctmp) = ZERO_REAL;
            Y(acctmp) = ZERO_REAL;
            Z(acctmp) = ZERO_REAL;
            break;
        case InvalidHalo:
        default:
            mw_fail("Invalid halo type in external acceleration\n");
    }
    mw_incaddv(acc, acctmp);

    /*Calculate the Bulge Accelerations*/
    switch (pot->sphere[0].type)
    {
        case HernquistSpherical:
            acctmp = hernquistSphericalAccel(&pot->sphere[0], pos, r);
            break;
        case PlummerSpherical:
            acctmp = plummerSphericalAccel(&pot->sphere[0], pos, r);
            break;
        case NoSpherical:
            X(acctmp) = ZERO_REAL;
            Y(acctmp) = ZERO_REAL;
            Z(acctmp) = ZERO_REAL;
            break;
        case InvalidSpherical:
        default:
            mw_fail("Invalid bulge type in external acceleration\n");
    }
    mw_incaddv(acc, acctmp);

    /*For debugging acceleration values*/
//    if (!isfinite(mw_absv(acc)))
//    {
//        mw_printf("ERROR: UNNATURAL ACCELERATION CALCULATED!\n");
//        mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
//        mw_printf("Acceleration = %.15f \n",mw_absv(acc));
//    }

    return acc;
}







