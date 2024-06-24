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

static const real pi = 3.1415926535;

static inline real lnfact(int n)
{
     int counter;
     real result = 0.0;
     if (n > 0)
     {
          for (counter = n; counter >= 1; counter--)
          {
               result += mw_log((real) counter);
          }
     }
     /*mw_printf("ln(%u!) = %.15f \n",n,result);*/
     return result;
}

static inline real binom(int n, int k)
{
    return mw_exp(lnfact(n)-lnfact(k)-lnfact(n-k));
}

static inline real leg_pol(real x, int l)
{
    real sum = 0.0;
    for (int m = 0; m < mw_floor(l/2)+1; m++)
    {
        sum += mw_pow(-1, m)*mw_pow(x, l-2*m)/mw_pow(2,l)*binom(l,m)*binom(2*(l-m),l);
    }
    if (x == 0.0)
    {
        if (l*1.0/2.0 == mw_ceil(l*1.0/2.0))
        {
            sum = mw_pow(-1,l/2)*mw_exp(lnfact(l)-2*lnfact(l/2) - l*mw_log(2));
        }
        else
        {
            sum = 0.0;
        }
    }
    /*mw_printf("P(%.15f, %u) = %.15f \n",x,l,sum);*/
    return sum;
}

static inline real leg_pol_derv(real x, int l)
{
    real sum = 0.0;
    for (int m = 0; m < mw_floor((l-1)/2)+1; m++)
    {
        sum += mw_pow(-1, m)*mw_pow(x, l - 2*m - 1)/mw_pow(2,l)*binom(l,m)*binom(2*(l-m),l)*(l - 2*m);
    }
    /*mw_printf("P'(%.15f, %u) = %.15f \n",x,l,sum);*/
    return sum;
}

static inline real lower_gamma(int n, real x)
{
    real sum = 0;
    for (int k = 0; k < n; k++)
    {
       sum += mw_pow(x,k)/mw_exp(lnfact(k));
    }
    /*mw_printf("g(%u, %.15f) = %.15f \n",n,x,mw_exp(lnfact(n-1))*(1-mw_exp(-x)*sum));*/
    return mw_exp(lnfact(n-1))*(1-mw_exp(-x)*sum);
}

static inline real GenExpIntegral(int n, real x) /*Optimized for convergence at n=1,x=0.5*/
{
    return mw_exp(-x)/(x + n/(1+1/(x + (n+1)/(1 + 2/(x + (n+2)/(1 + 3/(x + (n+3)/(1 + 4/(x + (n+4)/(1 + 5/(x + (n+5)/(1 + 6/(x + (n+6)/(1 + 7/(x + (n+7)/(1 + 8/(x + (n+8)/(1 + 9/(x + (n+9)/11)))))))))))))))))));
}

static inline real RExpIntegrand (real k, real R, real Rd, real z, real zd)
{
    real val = k*mw_cos(k*z)*(aExp(k,R,Rd)*besselK1(k*R) - bExp(k,R,Rd)*besselI1(k*R))/(sqr(zd*k) + 1);
    //mw_printf("RExp(%.15f,%.15f,%.15f,%.15f,%.15f) = %.15f\n",k,R,Rd,z,zd,val);
    return val;
}

static inline real ZExpIntegrand (real k, real R, real Rd, real z, real zd)
{
    real val = k*mw_sin(k*z)*(aExp(k,R,Rd)*besselK0(k*R) + bExp(k,R,Rd)*besselI0(k*R))/(sqr(zd*k) + 1);
    //mw_printf("ZExp = %.15f\n",val);
    return val;
}

static inline real RSechIntegrand (real k, real R, real Rd, real z, real zd)
{
    real val = sqr(k)*zd/mw_sinh(pi*k*zd/2.0)*mw_cos(k*z)*(aExp(k,R,Rd)*besselK1(k*R) - bExp(k,R,Rd)*besselI1(k*R));
    //mw_printf("RSech(%.15f,%.15f,%.15f,%.15f,%.15f) = %.15f\n",k,R,Rd,z,zd,val);
    return val;
}

static inline real ZSechIntegrand (real k, real R, real Rd, real z, real zd)
{
    real val = sqr(k)*zd/mw_sinh(pi*k*zd/2.0)*mw_sin(k*z)*(aExp(k,R,Rd)*besselK0(k*R) + bExp(k,R,Rd)*besselI0(k*R));
    //mw_printf("ZSech = %.15f\n",val);
    return val;
}

/**********************************************************************************************************************************************************************************************/

mwvector pointAccel(const mwvector pos, const mwvector pos1, const real mass)
{
    mwvector v = mw_subv(pos1, pos);
    real dist = mw_distv(pos, pos1);
    real tmp = mass/mw_pow(dist,3.0);
    mw_incmulvs(v, tmp);
    return v;
}

//Plummer accel
/*
mwvector plummerAccel(const mwvector pos, const mwvector pos1, const real mass, const real scale, const real scale2)
{
    //mw_printf("-----------------------------------------\n");
    //mw_printf("body position = {%.15f,%.15f,%.15f}\n", pos.x, pos.y, pos.z);
    //mw_printf("LMC position = {%.15f,%.15f,%.15f}\n", pos1.x, pos1.y, pos1.z);
    mwvector v = mw_subv(pos1, pos);
    real dist = mw_distv(pos, pos1);
    real tmp = mw_sqrt(mw_pow(scale,2.0) + mw_pow(dist,2.0));
    mw_incmulvs(v, mass/mw_pow(tmp,3.0));
    return v;
}
*/

//NFW accel from Erkal
/*
mwvector plummerAccel(const mwvector pos, const mwvector pos1, const real mass, const real scale, const real scale2)
{
    //mw_printf("-----------------------------------------\n");
    //mw_printf("body position = {%.15f,%.15f,%.15f}\n", pos.x, pos.y, pos.z);
    //mw_printf("LMC position = {%.15f,%.15f,%.15f}\n", pos1.x, pos1.y, pos1.z);    
    mwvector v = mw_subv(pos1, pos);
    //mw_printf("v = {%.15f,%.15f,%.15f}\n", v.x, v.y, v.z);
    real dist = mw_distv(pos, pos1);
    //mw_printf("dist = %.15f\n", dist);

    if (dist <= scale2)
    {
	real tmp1 = mw_log(1.0 + (dist/scale)) - (dist/(dist + scale));
	real tmp2 = mw_log(1.0 + (scale2/scale)) - (scale2/(scale2 + scale));
	real mass_LMC = mass*tmp1/tmp2;
	mw_incmulvs(v, mass_LMC/mw_pow(dist,3.0));
    }
    else
    {
        mw_incmulvs(v, mass/mw_pow(dist,3.0));
    }    

    //mw_printf("accel = {%.15f,%.15f,%.15f}\n", v.x, v.y, v.z);
    return v;
}
*/

//Hernquist accel

mwvector plummerAccel(const mwvector pos, const mwvector pos1, const real mass, const real scale, const real scale2)
{
    mwvector v = mw_subv(pos1, pos);
    real dist = mw_distv(pos, pos1);
    real tmp = mw_pow((scale + dist), 2.0);
    mw_incmulvs(v, mass/(tmp*dist));
    return v;
}


/*spherical bulge potentials*/

static inline mwvector hernquistSphericalAccel(const Spherical* sph, mwvector pos, real r)
{
    const real tmp = sph->scale + r;

    return mw_mulvs(pos, -sph->mass / (r * sqr(tmp)));
}

static inline mwvector plummerSphericalAccel(const Spherical* sph, mwvector pos, real r)
{
    const real tmp = mw_sqrt(sqr(sph->scale) + sqr(r));

    return mw_mulvs(pos, -sph->mass / mw_pow(tmp,3.0));
}

/* Disk Potentials */

static inline mwvector miyamotoNagaiDiskAccel(const Disk* disk, mwvector pos, real r)
{
    mwvector acc;
    const real a   = disk->scaleLength;
    const real b   = disk->scaleHeight;
    const real zp  = mw_pow(mw_pow(Z(pos),2.0) + mw_pow(b,2.0), 0.5);
    const real azp = a + zp;

    const real rp  = mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0) + mw_pow(azp,2.0);
    const real rth = mw_pow(rp,1.5);  /* rp ^ (3/2) */

    X(acc) = -disk->mass * X(pos) / rth;
    Y(acc) = -disk->mass * Y(pos) / rth;
    Z(acc) = -disk->mass * Z(pos) * azp / (zp * rth);

    //mw_printf("Acceleration[AX,AY,AZ] = [%.15f,%.15f,%.15f]\n",X(acc),Y(acc),Z(acc));
    
    /*mw_printf("disk x acc: %.15f\n", X(acc));
    mw_printf("disk y acc: %.15f\n", Y(acc));
    mw_printf("disk z acc: %.15f\n", Z(acc));*/
    return acc;
}

/*WARNING: This potential uses incomplete gamma functions and the generalized exponential integral function. This potential will take longer than other potentials to run.*/
static inline mwvector freemanDiskAccel(const Disk* disk, mwvector pos, real r)     /*From Freeman 1970*/
{
    mwvector acc;
    /*Reset acc vector*/
    X(acc) = 0.0;
    Y(acc) = 0.0;
    Z(acc) = 0.0;

    const mwvector r_hat = mw_mulvs(pos, 1.0/r);
    const real r_proj = mw_sqrt(sqr(X(pos)) + sqr(Y(pos)));

    mwvector theta_hat;
    X(theta_hat) = X(pos)*Z(pos)/r/r_proj;
    Y(theta_hat) = Y(pos)*Z(pos)/r/r_proj;
    Z(theta_hat) = -r_proj/r;

    const real scl = disk->scaleLength;
    const real M = disk->mass;
    const real costheta = Z(pos)/r;
    const real sintheta = r_proj/r;

    real r_f = 0.0;
    real theta_f = 0.0;
    real a_r = 0.0;
    real b_r = 0.0;
    real p0 = 0.0;
    real pcos = 0.0;

    for (int l = 0; l < 5; l++) /*Only even terms add, so counter is multiplied by 2*/
    {
        if (l>0)
        {
            a_r = GenExpIntegral(2*l,r/scl);
        }
        else
        {
            a_r = 0.0;
        }
        b_r = mw_pow(scl/r,2*l+2)*lower_gamma(2*l+2,r/scl);
        p0 = leg_pol(0,2*l);
        pcos = leg_pol(costheta,2*l);
        r_f += p0*pcos*(2*l*a_r - (2*l+1)*b_r);

        if (l > 0)
        {
            theta_f -= p0*(a_r+b_r)*sintheta*leg_pol_derv(costheta,2*l);
        }
    }

    mwvector r_comp = mw_mulvs(r_hat, r_f*M/sqr(scl));
    mwvector theta_comp = mw_mulvs(theta_hat, theta_f*M/sqr(scl));

    X(acc) = X(r_comp) + X(theta_comp);
    Y(acc) = Y(r_comp) + Y(theta_comp);
    Z(acc) = Z(r_comp) + Z(theta_comp);

    //if (r_f > 0)
    //{
    //    mw_printf("ERROR: Repulsive acceleration!\n");
    //    mw_printf("r_magnitude = %.15f\n",r_f*M/sqr(scl));
    //    mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    //    mw_printf("Acceleration = %.15f \n",mw_absv(acc));
    //}

    return acc;
}

/*WARNING: This potential can take a while to integrate if any part of the orbit extends past 100 times the scaleLength*/
static inline mwvector doubleExponentialDiskAccel(const Disk* disk, mwvector pos, real r)
{
    //mw_printf("Calculating Acceleration\n");
    //mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    mwvector acc;

    const real R = mw_sqrt(sqr(X(pos)) + sqr(Y(pos)));
    mwvector R_hat;
    X(R_hat) = X(pos) / R;
    Y(R_hat) = Y(pos) / R;
    Z(R_hat) = 0.0;

    mwvector Z_hat;
    X(Z_hat) = 0.0;
    Y(Z_hat) = 0.0;
    Z(Z_hat) = Z(pos) / mw_abs(Z(pos));

    const real Rd = disk->scaleLength;
    const real zd = disk->scaleHeight;
    const real M = disk->mass;
    const real z = Z(pos);
    const real h = 0.001;

    const real a = 1.0 / Rd;
    const real b = 1.0 / zd;
    
    real R_piece = 0.0;
    real z_piece = 0.0;

    real j0_zero;
    real psi_in_0;
    real psi_0;
    real psi_prime_0;
    real j0_x;
    real j0_w;
    real fun_0;

    real j1_zero;
    real psi_in_1;
    real psi_1;
    real psi_prime_1;
    real j1_x;
    real j1_w;
    real fun_1;

    for (int n = 0; n < 150; n+=1)    // Hidenori Ogata's Numerical Integration Formula Based on the Bessel Functions
    {
        j0_zero = besselJ0_zero(n)/pi;
        j1_zero = besselJ1_zero(n)/pi;
        psi_in_0 = h * j0_zero;
        psi_in_1 = h * j1_zero;
        psi_0 = psi_in_0 * mw_sinh(pi / 2.0 * mw_sinh(psi_in_0)) / mw_cosh(pi / 2.0 * mw_sinh(psi_in_0));
        psi_1 = psi_in_1 * mw_sinh(pi / 2.0 * mw_sinh(psi_in_1)) / mw_cosh(pi / 2.0 * mw_sinh(psi_in_1));
        psi_prime_0 = (mw_sinh(pi * mw_sinh(psi_in_0)) + pi * psi_in_0 * mw_cosh(psi_in_0)) / (mw_cosh(pi * mw_sinh(psi_in_0)) + 1.0);
        psi_prime_1 = (mw_sinh(pi * mw_sinh(psi_in_1)) + pi * psi_in_1 * mw_cosh(psi_in_1)) / (mw_cosh(pi * mw_sinh(psi_in_1)) + 1.0);
        j0_x = pi / h * psi_0;
        j1_x = pi / h * psi_1;

        j0_w = 2.0 / (pi * j0_zero * mw_pow(besselJ1(pi * j0_zero), 2)) * besselJ0(j0_x) * psi_prime_0;
        j1_w = 2.0 / (pi * j1_zero * mw_pow(besselJ2(pi * j1_zero),2)) * besselJ1(j1_x) * psi_prime_1;

        fun_1 = j1_x * mw_pow(mw_pow(a,2.0) + mw_pow(j1_x / R, 2.0), -1.5) * (b * mw_exp(-j1_x / R * mw_abs(z)) - j1_x / R * mw_exp(-b * mw_abs(z))) / (mw_pow(b, 2.0) - mw_pow(j1_x / R, 2.0));
        fun_0 = mw_pow(mw_pow(a, 2.0) + mw_pow(j0_x / R, 2.0), -1.5) * j0_x / R * (mw_exp(-j0_x / R * mw_abs(z)) - mw_exp(-b * mw_abs(z))) / (mw_pow(b, 2.0) - mw_pow(j0_x / R, 2.0));

        real z_pieceAdd = (4.0 * pi * a * b / R) * fun_0 * j0_w;
        real R_pieceAdd = (4.0 * pi * a / mw_pow(R, 2.0)) * fun_1 * j1_w;
        
        z_piece += z_pieceAdd;
        R_piece += R_pieceAdd;
    }

    mwvector R_comp = mw_mulvs(R_hat, -M * mw_pow(a,2) * b / 4 / pi * R_piece);
    mwvector Z_comp = mw_mulvs(Z_hat, -M * mw_pow(a, 2) * b / 4 / pi * z_piece);

    X(acc) = X(R_comp) + X(Z_comp);
    Y(acc) = Y(R_comp) + Y(Z_comp);
    Z(acc) = Z(R_comp) + Z(Z_comp);

    real magnitude = mw_sqrt(sqr(X(acc))+sqr(Y(acc))+sqr(Z(acc)));
    //mw_printf("Acceleration[AX,AY,AZ] = [%.15f,%.15f,%.15f]   Magnitude = %.15f\n",X(acc),Y(acc),Z(acc),magnitude);
    return acc;
}

/*WARNING: This potential can take a while to integrate if any part of the orbit extends past 100 times the scaleLength*/
static inline mwvector sech2ExponentialDiskAccel(const Disk* disk, mwvector pos, real r)
{
    //mw_printf("Calculating Acceleration\n");
    //mw_printf("[X,Y,Z] = [%.15f,%.15f,%.15f]\n",X(pos),Y(pos),Z(pos));
    mwvector acc;

    const real R = mw_sqrt(sqr(X(pos)) + sqr(Y(pos)));
    mwvector R_hat;
    X(R_hat) = X(pos)/R;
    Y(R_hat) = Y(pos)/R;
    Z(R_hat) = 0.0;

    mwvector Z_hat;
    X(Z_hat) = 0.0;
    Y(Z_hat) = 0.0;
    Z(Z_hat) = 1.0;

    const real Rd = disk->scaleLength;
    const real zd = disk->scaleHeight;
    const real M = disk->mass;
    const real z = Z(pos);

    const int n = 15;
    const real a = 0.0;
    const real b = 60.0/R;
    real integralR = 0.0;
    real integralZ = 0.0;
    const real weight[] = {0.236927,0.478629,0.568889,0.478629,0.236927};
    const real point[] = {-0.90618,-0.538469,0.0,0.538469,0.90618};
    const real h = (b-a)/(n*1.0);

    for (int k = 0; k < n; k++)     /*Five-point Gaussian Quadrature*/
    {
        real Rpiece = 0.0;
        real Zpiece = 0.0;
        real k_val = 0.0;
        for (int j = 0; j < 5; j++)
        {
            k_val = h*point[j]/2 + a+(k*1.0+0.5)*h;
            Rpiece = Rpiece + h*weight[j]*RSechIntegrand(k_val,R,Rd,z,zd)/2.0;
            Zpiece = Zpiece + h*weight[j]*ZSechIntegrand(k_val,R,Rd,z,zd)/2.0;
        }
        integralR  = integralR + Rpiece;
        integralZ  = integralZ + Zpiece;
    }

    mwvector R_comp = mw_mulvs(R_hat, -1.0*M*integralR);
    mwvector Z_comp = mw_mulvs(Z_hat, -1.0*M*integralZ);

    X(acc) = X(R_comp) + X(Z_comp);
    Y(acc) = Y(R_comp) + Y(Z_comp);
    Z(acc) = Z(R_comp) + Z(Z_comp);

    //real magnitude = mw_sqrt(sqr(X(acc))+sqr(Y(acc))+sqr(Z(acc)));

    //mw_printf("Acceleration[AX,AY,AZ] = [%.15f,%.15f,%.15f]   Magnitude = %.15f\n",X(acc),Y(acc),Z(acc),magnitude);

    return acc;
}

//Softened needle bar potential
static inline mwvector orbitingBarAccel(const Disk* disk, mwvector pos, real r, real time)
{
    real amp = disk->mass;
    real a = disk->scaleLength;
    real b = 1.4;//Triaxial softening length
    real c = 1;//Prolate softening length
    
    real curAngle = (disk->patternSpeed * time * -1)+disk->startAngle;
    //first rotate pos curAngle * -1 radians to emulate the current angle of the bar
    real Radi = mw_sqrt(pos.x*pos.x+pos.y*pos.y);
    real Phi = mw_atan(pos.y/pos.x);
    Phi -= curAngle;
    if(pos.x < 0){
        Radi = Radi * -1;
    }
    real x = Radi*cos(Phi);
    real y = Radi*sin(Phi); 
    real z = pos.z;

    //calculate force in accordance with the galpy implementation
    real secondpart = mw_pow(y,2.) + mw_pow(b+mw_sqrt(mw_pow(c,2) + mw_pow(z,2)),2);
    real Tp = mw_sqrt(mw_pow(a + x,2) + secondpart);
    real Tm = mw_sqrt(mw_pow(a - x,2) + secondpart);

    mwvector force;
    force.x = -2.*x/Tp/Tm/(Tp+Tm);
    force.y = -y/2./Tp/Tm*(Tp+Tm-4.*mw_pow(x,2.)/(Tp+Tm))/(mw_pow(y,2.)+mw_pow(b+mw_sqrt(mw_pow(z,2.)+mw_pow(c,2)),2.));
    force.z = force.y*z/y*(b+mw_sqrt(mw_pow(z,2)+mw_pow(c,2)))/mw_sqrt(mw_pow(z,2)+mw_pow(c,2));
    
    //undo the pos rotation and calculate acceleration from the force vector we got
    mwvector acc;
    real cp = cos (curAngle);
    real sp = sin (curAngle);
    acc.x=cp*force.x-sp*force.y;
    acc.y=sp*force.x+cp*force.y;
    acc.z=force.z;
    
    acc.x = acc.x*amp;
    acc.y = acc.y*amp;
    acc.z = acc.z*amp;

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

    const real v0 = halo->vhalo;
    const real q  = halo->flattenZ;
    const real d  = halo->scaleLength;

    const real denom = mw_pow(d,2.0) + mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0) + mw_pow(Z(pos)/q,2.0);
    const real k = -2.0*v0*v0/denom;

    X(acc) = k * X(pos);
    Y(acc) = k * Y(pos);
    Z(acc) = k * Z(pos)/(q*q);

    return acc;
}

static inline mwvector NFWerkalHaloAccel(const Halo* halo, mwvector pos)
{
    mwvector acc;

    const real a = halo->scaleLength;
    const real q  = halo->flattenZ;
    const real M  = halo->mass;

    real c = M/(mw_log(1.0 + 15.3) - (15.3/(1.0 + 15.3)));  //c = 15.3
    real d = mw_sqrt(mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0) + (mw_pow(Z(pos),2.0)/mw_pow(q,2.0)));
    real tmp1 = (a*mw_pow(d,2.0)) + mw_pow(d,3.0);
    real tmp2 = mw_log(1.0 + (d/a));

    X(acc) = c * ((X(pos)/tmp1) - (X(pos)*tmp2/mw_pow(d,3.0)));
    Y(acc) = c * ((Y(pos)/tmp1) - (Y(pos)*tmp2/mw_pow(d,3.0)));
    Z(acc) = c * ((Z(pos)/(q*q)/tmp1) - (Z(pos)*tmp2/(q*q)/mw_pow(d,3.0)));

    return acc;
}


static inline mwvector SphericalNFWerkalHaloAccel(const Halo* halo, mwvector pos, real r)
{
    const real a = halo->scaleLength;
    const real M  = halo->mass;

    const real ar = a + r;
    const real tmp = M/(mw_log(1.0 + 15.3) - (15.3/(1.0 + 15.3)));  //c = 15.3
    const real c = tmp*(r - (ar*mw_log(1 + (r/a))))/(r*r*r*ar);							

    return mw_mulvs(pos, c);
}

static inline mwvector nfwHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real a = h->scaleLength;
    const real M = sqr(h->vhalo)*a/0.2162165954; /*Maximum of [ln(1+x)/x - 1/(1+x)]*/
    const real ar = a + r;

    const real c = (-M/sqr(r))*(mw_log(ar/a)/r - 1/ar);

    return mw_mulvs(pos, c);
}

/* CHECKME: Seems to have precision related issues for a small number of cases for very small qy */
static inline mwvector triaxialHaloAccel(const Halo* h, mwvector pos, real r)  /** Triaxial Logarithmic **/
{
    mwvector acc;

    /* TODO: More things here can be cached */
    const real qzs      = sqr(h->flattenZ);
    const real rhalosqr = sqr(h->scaleLength);
    const real mvsqr    = -sqr(h->vhalo);

    const real xsqr = sqr(X(pos));
    const real ysqr = sqr(Y(pos));
    const real zsqr = sqr(Z(pos));

    const real arst  = rhalosqr + (h->c1 * xsqr) + (h->c3 * X(pos) * Y(pos)) + (h->c2 * ysqr);
    const real arst2 = (zsqr / qzs) + arst;

    X(acc) = mvsqr * (((2.0 * h->c1) * X(pos)) + (h->c3 * Y(pos)) ) / arst2;

    Y(acc) = mvsqr * (((2.0 * h->c2) * Y(pos)) + (h->c3 * X(pos)) ) / arst2;

    Z(acc) = (2.0 * mvsqr * Z(pos)) / ((qzs * arst) + zsqr);

    return acc;
}

static inline mwvector ASHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real gam = h->gamma;
    const real lam = h->lambda;
    const real M = h->mass;
    const real a = h->scaleLength;
    const real scaleR = r/a;
    const real scaleL = lam/a;
    real c;

    if (r<lam)
    {    c = -(M/(a*r))*mw_pow(scaleR,gam-1.0)/(1.0+mw_pow(scaleR,gam-1.0));
    }
    else
    {    c = -(M/mw_pow(r,2.0))*mw_pow(scaleL,gam)/(1.0+mw_pow(scaleL,gam-1.0));
    }

    return mw_mulvs(pos, c/r);
}

static inline mwvector WEHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;
    const real sum2 = sqr(a) + sqr(r);

    const real c = (-M/r)*(a + mw_sqrt(sum2))/(sum2 + a*mw_sqrt(sum2));

    return mw_mulvs(pos, c/r);

}

static inline mwvector NFWMHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;
    const real ar = a + r;

    const real c = (-M/sqr(r))*(mw_log(ar/a)/r - 1/ar);

    return mw_mulvs(pos, c);

}

static inline mwvector plummerHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real tmp = mw_sqrt(sqr(h->scaleLength) + sqr(r));

    return mw_mulvs(pos, -h->mass / cube(tmp));
}

static inline mwvector hernquistHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real tmp = h->scaleLength + r;

    return mw_mulvs(pos, -h->mass / (r * sqr(tmp)));
}

static inline mwvector ninkovicHaloAccel(const Halo* h, mwvector pos, real r)      /*Special case of Ninkovic Halo (l1=0,l2=3,l3=2) (Ninkovic 2017)*/
{
    const real rho0 = h->rho0;
    const real a = h->scaleLength;
    const real lambda = h->lambda;

    const real z = r/a;
    const real zl = lambda/a;
    const real f = 4.0*pi/3.0*rho0*cube(a);

    real mass_enc;

    if (r > lambda)
    {
        mass_enc = f*(mw_log(1.0+cube(zl)) - cube(zl)/(1+cube(zl)));
    }
    else
    {
        mass_enc = f*(mw_log(1.0+cube(z)) - cube(z)/(1+cube(zl)));
    }

    mwvector acc = mw_mulvs(pos, -mass_enc/cube(r));

    return acc;
}

mwvector nbExtAcceleration(const Potential* pot, mwvector pos, real time)
{
    mwvector acc, acctmp;
    real limit = mw_pow(2.0,-8.0);

    /* Change r if less than limit. Done this way to pipeline this step*/
    real r = (mw_absv(pos) <= limit)*limit + (mw_absv(pos) > limit)*mw_absv(pos);

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
            X(acc) = 0.0;
            Y(acc) = 0.0;
            Z(acc) = 0.0;
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
            X(acctmp) = 0.0;
            Y(acctmp) = 0.0;
            Z(acctmp) = 0.0;
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
	case NFWerkalHalo:
	    acctmp = NFWerkalHaloAccel(&pot->halo, pos);
	    break;
	case SphericalNFWerkalHalo:
	    acctmp = SphericalNFWerkalHaloAccel(&pot->halo, pos, r);
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
            X(acctmp) = 0.0;
            Y(acctmp) = 0.0;
            Z(acctmp) = 0.0;
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
            X(acctmp) = 0.0;
            Y(acctmp) = 0.0;
            Z(acctmp) = 0.0;
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







