/*********************** DISCLAIMER *************************
* This code is included for its potential usefulness but    *
* has not been fully tested. It is not currently used in    *
* any other milkyway@home files. Use at your own risk.      *
************************************************************/

#include "nbody_priv.h"
#include "nbody_potential.h"
#include "nbody_density.h"
#include "nbody_potential_types.h"
#include "milkyway_util.h"
#include "nbody_caustic.h"
#include "nbody_bessel.h"
#include "nbody_show.h"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

/*Spherical Buldge Densities*/
static inline real hernquistSphericalDensity(const Spherical* sph, real* r)
{
    const real a = mw_real_var(sph->scale, BULGE_RADIUS_POS);
    const real M = mw_real_var(sph->mass, BULGE_MASS_POS);
    real tmp1, tmp2;

    /*return 0 rather than get a divide by 0 error*/
    if(showRealValue(r) == 0) {
        return ZERO_REAL;
    }

    tmp1 = mw_mul(&M,&a);
    tmp2 = mw_add(r,&a);
    tmp2 = cube(&tmp2);
    tmp2 = mw_mul(r, &tmp2);
    tmp1 = mw_div(&tmp1, &tmp2);
    return mw_mul_s(&tmp1, inv_0(2*M_PI));

  /*return ((sph->mass)/(2*pi*mw_pow(a, 3)))*mw_pow(a, 4)/(r*mw_pow(r+a, 3));*/
}

static inline real plummerSphericalDensity(const Spherical* sph, real* r)
{
    const real a = mw_real_var(sph->scale, BULGE_RADIUS_POS);
    const real M = mw_real_var(sph->mass, BULGE_MASS_POS);
    real tmp1, tmp2;
    
    if(showRealValue(&a) == 0)
    {
        return ZERO_REAL;
    }
    tmp1 = mw_div(r,&a);
    tmp1 = sqr(&tmp1);
    tmp1 = mw_add_s(&tmp1, 1.0);
    tmp1 = minusfivehalves(&tmp1);

    tmp2 = cube(&a);
    tmp2 = mw_div(&M, &tmp2);
    tmp2 = mw_mul_s(&tmp2, 3.0/(4.0*M_PI));

    return mw_mul(&tmp2, &tmp1);

  /*real r_a = r/a;
    real rho_peak = 3*M/(4.0*pi*mw_pow(a, 3.0));
    return rho_peak*mw_pow((1.0+mw_pow(r_a,2.0)), -2.5);*/
}

/*Disk Densities*/
static inline real miyamotoNagaiDiskDensity(const Disk* disk, mwvector* pos)
{
    const real M   = mw_real_var(disk->mass, DISK_MASS_POS);
    const real a   = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);
    const real b   = mw_real_var(disk->scaleHeight, DISK_HEIGHT_POS);

    real R   = mw_hypot(&X(pos), &Y(pos));
    real zp  = mw_hypot(&Z(pos), &b);
    real azp = mw_add(&a, &zp);

    real tmp1, tmp2, tmp3;

    tmp1 = sqr(&R);
    tmp1 = mw_mul(&a, &tmp1);
    tmp2 = mw_mul_s(&zp, 3.0);
    tmp2 = mw_add(&tmp2, &a);
    tmp3 = sqr(&azp);
    tmp2 = mw_mul(&tmp2, &tmp3);
    real numer = mw_add(&tmp1, &tmp2);

    tmp1 = sqr(&R);
    tmp2 = sqr(&azp);
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp1 = fivehalves(&tmp1);
    tmp2 = cube(&zp);
    real denom = mw_mul(&tmp1, &tmp2);

    tmp1 = sqr(&b);
    tmp1 = mw_mul(&M, &tmp1);
    real factor = mw_mul_s(&tmp1, inv_0(4.0*M_PI));

    tmp1 = mw_div(&numer, &denom);
    real den = mw_mul(&factor, &tmp1);

  /*const real R   = mw_pow(mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0), 0.5);
    const real zp  = mw_pow(mw_pow(Z(pos),2.0) + mw_pow(b,2.0), 0.5);
    const real azp = a + zp;
    real numer = M*b*b*(a*R*R + (a + 3.0*zp) * mw_pow(azp, 2));
    real denom = 4.0*pi*mw_pow(R*R + azp*azp, 2.5)*mw_pow(zp, 3.0);
    return numer/denom;*/

    return den;
}

static inline real doubleExponentialDiskDensity(const Disk* disk, mwvector* pos)
{
    const real M   = mw_real_var(disk->mass, DISK_MASS_POS);
    const real d_r = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);
    const real d_z = mw_real_var(disk->scaleHeight, DISK_HEIGHT_POS);
    const real R   = mw_hypot(&X(pos), &Y(pos));
    real tmp1, tmp2, tmp3;

    tmp1 = sqr(&d_r);
    tmp1 = mw_mul(&d_z,&tmp1);
    tmp1 = mw_div(&M, &tmp1);
    tmp2 = mw_div(&R,&d_r);
    tmp3 = mw_abs(&Z(pos));
    tmp3 = mw_div(&tmp3,&d_z);
    tmp2 = mw_add(&tmp2, &tmp3);
    tmp2 = mw_exp(&tmp2);
    tmp1 = mw_div(&tmp1, &tmp2);

    return mw_mul_s(&tmp1, inv_0(4.0*M_PI));

  /*const real R   = mw_sqrt(mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0));
    real den = M/(4.0*pi*d_z*mw_pow(d_r,2.0))*mw_exp(-R/d_r)*mw_exp(-mw_abs(Z(pos))/d_z);
    return den;*/
}

static inline real sech2ExponentialDiskDensity(const Disk* disk, mwvector* pos)
{
    const real M   = mw_real_var(disk->mass, DISK_MASS_POS);
    const real d_r = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);
    const real d_z = mw_real_var(disk->scaleHeight, DISK_HEIGHT_POS);
    const real R   = mw_hypot(&X(pos), &Y(pos));
    real tmp1, tmp2;

    tmp1 = sqr(&d_r);
    tmp1 = mw_mul(&d_z, &tmp1);
    tmp1 = mw_div(&M, &tmp1);
    tmp2 = mw_div(&R, &d_r);
    tmp2 = mw_exp(&tmp2);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_div(&Z(pos), &d_z);
    tmp2 = mw_cosh(&tmp2);
    tmp2 = sqr(&tmp2);
    tmp1 = mw_div(&tmp1, &tmp2);

    return mw_mul_s(&tmp1, inv_0(4.0*M_PI));

  /*const real R   = mw_sqrt(mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0));
    real den = M/(4.0*pi*d_z*d_r*d_r)*mw_exp(-R/d_r)/mw_pow(mw_cosh(Z(pos)/d_z),2.0);
    return den;*/
}

static inline real orbitingBarDensity(const Disk* disk, mwvector* pos, real_0 time)
{
    real tmp1, tmp2;
    real M = mw_real_var(disk->mass, DISK_MASS_POS);
    real a = mw_real_var(disk->scaleLength, DISK_LENGTH_POS);  //Bar half-length
    real b = mw_real_const(1.4);                               //Triaxial softening length
    real c = mw_real_const(1.0);                               //Prolate softening length

    real pSpeed = mw_real_var(disk->patternSpeed, 0); //Several of these parameters are not assigned in AUTODIFF, so we use 0 as a placeholder//
    real sAngle = mw_real_var(disk->startAngle, 0);

    tmp1 = mw_mul_s(&pSpeed, time);
    real curAngle =  mw_sub(&sAngle, &tmp1);
    //first rotate pos curAngle * -1 radians to emulate the current angle of the bar
    real Radi = mw_hypot(&X(pos), &Y(pos));
    real Phi = mw_atan2(&Y(pos), &X(pos));
    Phi = mw_sub(&Phi, &curAngle);

    tmp1 = mw_cos(&Phi);
    real x = mw_mul(&Radi, &tmp1);
    tmp1 = mw_sin(&Phi);
    real y = mw_mul(&Radi, &tmp1); 
    real z = Z(pos);

    real zc = mw_hypot(&z,&c);
    tmp1 = mw_add(&b,&zc);
    real bzc2 = sqr(&tmp1);

    tmp1 = sqr(&y);
    tmp1 = mw_mul(&b,&tmp1);
    tmp2 = mw_mul_s(&zc,3.0);
    tmp2 = mw_add(&b, &tmp2);
    tmp2 = mw_mul(&tmp2, &bzc2);
    real bigA = mw_add(&tmp1, &tmp2);
    tmp1 = sqr(&y);
    real bigC = mw_add(&tmp1, &bzc2);

    tmp1 = mw_mul(&bigA,&bigC);
    real AC3 = mw_mul_s(&tmp1, 3.0);
    tmp1 = mw_mul_s(&bigA,2.0);
    tmp2 = mw_mul(&b,&bigC);
    real A2bC = mw_add(&tmp1, &tmp2);

    tmp1 = sqr(&c);
    tmp1 = mw_div(&tmp1,&a);
    tmp2 = sqr(&bigC);
    tmp1 = mw_div(&tmp1,&tmp2);
    tmp2 = cube(&zc);
    tmp1 = mw_div(&tmp1,&tmp2);
    real factor = mw_mul_s(&tmp1,inv_0(24*M_PI));

    tmp1 = mw_add(&x,&a);
    tmp2 = mw_add(&x,&a);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_mul(&A2bC, &tmp2);
    tmp2 = mw_add(&AC3, &tmp2);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = mw_add(&x,&a);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_add(&bigC,&tmp2);
    tmp2 = threehalves(&tmp2);
    real part1 = mw_div(&tmp1,&tmp2);


    tmp1 = mw_sub(&x,&a);
    tmp2 = mw_sub(&x,&a);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_mul(&A2bC, &tmp2);
    tmp2 = mw_add(&AC3, &tmp2);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = mw_sub(&x,&a);
    tmp2 = sqr(&tmp2);
    tmp2 = mw_add(&bigC,&tmp2);
    tmp2 = threehalves(&tmp2);
    real part2 = mw_div(&tmp1,&tmp2);

    part1 = mw_sub(&part1, &part2);

    real unscaledDens = mw_mul(&factor, &part1);

    return mw_mul(&unscaledDens, &M);

  /*real curAngle = (disk->patternSpeed * time * -1)+disk->startAngle;
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

    real zc = mw_sqrt(mw_pow(z,2)+mw_pow(c,2));
    real bzc2 = mw_pow(b+zc,2);
    real bigA = b*mw_pow(y,2) + (b+3*zc)*bzc2;
    real bigC = mw_pow(y,2)+bzc2;
    real unscaledDens = mw_pow(c,2)/24/pi/a/mw_pow(bigC,2)/mw_pow(zc,3)*
    ((x+a)*(3*bigA*bigC+(2*bigA+b*bigC)*mw_pow(x+a,2))/
    mw_pow(bigC+mw_pow(x+a,2),1.5)-(x-a)*(3*bigA*bigC+(2*bigA+b*bigC)*
    mw_pow(x-a,2))/mw_pow(bigC+mw_pow(x-a,2),1.5));

    return unscaledDens * disk->mass;*/
}

/*Halo Densities*/
static inline real logarithmicHaloDensity(const Halo* h, mwvector* pos) /** flattenZ should be greater than 1/sqrt(2) to keep positive definite **/
{
    real tmp1, tmp2, tmp3;

    const real v  = mw_real_var(h->vhalo, HALO_MASS_POS);
    const real q  = mw_real_var(h->flattenZ, HALO_ZFLATTEN_POS);
    const real a  = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real R = mw_hypot(&X(pos), &Y(pos));

    tmp1 = sqr(&q);
    tmp1 = mw_mul_s(&tmp1, 2.0);
    tmp1 = mw_add_s(&tmp1, 1.0);
    tmp2 = sqr(&a);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp2 = sqr(&R);
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp2 = sqr(&q);
    tmp2 = inv(&tmp2);
    tmp2 = mw_neg(&tmp2);
    tmp2 = mw_add_s(&tmp2, 2.0);
    tmp3 = sqr(&Z(pos));
    tmp2 = mw_mul(&tmp3, &tmp2);
    real numer = mw_add(&tmp1, &tmp2);

    tmp1 = mw_div(&Z(pos), &q);
    tmp1 = sqr(&tmp1);
    tmp2 = sqr(&a);
    tmp1 = mw_add(&tmp2, &tmp1);
    tmp2 = sqr(&R);
    tmp1 = mw_add(&tmp2, &tmp1);
    tmp2 = sqr(&q);
    tmp1 = sqr(&tmp1);
    real denom = mw_mul(&tmp2, &tmp1);

    real density = sqr(&v);
    density = mw_mul(&density, &numer);
    density = mw_div(&density, &denom);
    density = mw_mul_s(&density, 0.5/M_PI);
    return density;

  /*const real R2 = mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0);

    real numer = (2.0*q*q + 1.0)*a*a + R2 + mw_pow(Z(pos),2.0)*(2.0-1.0/q/q);
    real denom = q*q*mw_pow(R2 + a*a + mw_pow(Z(pos)/q,2.0),2.0);

    return v*v*numer/2.0/pi/denom;*/
}

static inline real NFWHaloDensity(const Halo* h,  real* r)
{
    if(showRealValue(r) == 0) return ZERO_REAL;

    real tmp1, tmp2;
    const real v  = mw_real_var(h->vhalo, HALO_MASS_POS);
    const real a  = mw_real_var(h->scaleLength, HALO_RADIUS_POS);

    tmp1 = mw_div(&v,&a);
    tmp1 = sqr(&tmp1);
    tmp1 = mw_mul_s(&tmp1, inv_0(4.0*M_PI*0.2162165954));
    tmp2 = mw_div(r,&a);
    tmp1 = mw_div(&tmp1, &tmp2);
    tmp2 = mw_div(r,&a);
    tmp2 = mw_add_s(&tmp2, 1.0);
    tmp2 = sqr(&tmp2);
    
    return mw_div(&tmp1, &tmp2);

  /*real rho = v*v/4.0/pi/a/a/0.2162165954;
    return rho / (r/a) / mw_pow(1.0+(r/a),2.0);*/
}

static inline real triaxialHaloDensity(const Halo* h, mwvector* pos)
{
    real tmp1, tmp2;
    const real v  = mw_real_var(h->vhalo, HALO_MASS_POS);
    const real a  = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real q  = mw_real_var(h->flattenZ, HALO_ZFLATTEN_POS);

    tmp1 = sqr(&a);
    tmp2 = sqr(&X(pos));
    tmp2 = mw_mul_s(&tmp2, (h->c1));
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp2 = sqr(&Y(pos));
    tmp2 = mw_mul_s(&tmp2, (h->c2));
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp2 = mw_mul(&X(pos),&Y(pos));
    tmp2 = mw_mul_s(&tmp2, (h->c3));
    tmp1 = mw_add(&tmp1, &tmp2);
    tmp2 = mw_div(&Z(pos),&q);
    tmp2 = sqr(&tmp2);
    const real D = mw_add(&tmp1, &tmp2);

    tmp1 = mw_mul_s(&D, (2.0*(h->c1) + 2.0*(h->c2)));
    tmp2 = sqr(&q);
    tmp2 = mw_div(&D, &tmp2);
    tmp2 = mw_mul_s(&tmp2, 2.0);
    const real num1 = mw_add(&tmp1, &tmp2);

    tmp1 = mw_mul_s(&X(pos), 2.0*(h->c1));
    tmp2 = mw_mul_s(&Y(pos), (h->c3));
    tmp1 = mw_add(&tmp1, &tmp2);
    const real num2 = sqr(&tmp1);

    tmp1 = mw_mul_s(&Y(pos), 2.0*(h->c2));
    tmp2 = mw_mul_s(&X(pos), (h->c3));
    tmp1 = mw_add(&tmp1, &tmp2);
    const real num3 = sqr(&tmp1);

    tmp1 = sqr(&q);
    tmp1 = mw_div(&Z(pos), &tmp1);
    tmp1 = mw_mul_s(&tmp1, 2.0);
    const real num4 = sqr(&tmp1);

    tmp1 = mw_add(&num2, &num3);
    tmp1 = mw_add(&tmp1, &num4);
    const real num = mw_sub(&num1, &tmp1);

    tmp1 = sqr(&v);
    tmp1 = mw_mul(&tmp1,&num);
    tmp2 = sqr(&D);
    tmp1 = mw_div(&tmp1,&tmp2);

    return mw_mul_s(&tmp1, inv_0(4.0*M_PI));

  /*const real D   = a*a + (h->c1)*mw_pow(X(pos),2.0) + (h->c2)*mw_pow(Y(pos),2.0) + (h->c3)*X(pos)*Y(pos) + mw_pow(Z(pos),2.0)/q/q;
    const real num = 2.0*(h->c1)*D + 2.0*(h->c2)*D + 2.0*D/q/q - mw_pow(2.0*(h->c1)*X(pos) + (h->c3)*Y(pos),2.0) - mw_pow(2.0*(h->c2)*Y(pos) + (h->c3)*X(pos),2.0) - mw_pow(2.0*Z(pos)/q/q,2.0);

    return v*v*num/4.0/pi/D/D;*/
}

static inline real hernquistHaloDensity(const Halo* h,  real* r)
{
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    real tmp1, tmp2;

    /*return 0 rather than get a divide by 0 error*/
    if(showRealValue(r) == 0) {
        return ZERO_REAL;
    }

    tmp1 = mw_mul(&M,&a);
    tmp2 = mw_add(r,&a);
    tmp2 = cube(&tmp2);
    tmp2 = mw_mul(r, &tmp2);
    tmp1 = mw_div(&tmp1, &tmp2);
    return mw_mul_s(&tmp1, inv_0(2*M_PI));

  /*return ((h->mass)/(2*pi*mw_pow(a, 3)))*mw_pow(a, 4)/(r*mw_pow(r+a, 3));*/
}

static inline real plummerHaloDensity(const Halo* h, real* r)
{
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    real tmp1, tmp2;
    
    if(showRealValue(&a) == 0)
    {
        return ZERO_REAL;
    }
    tmp1 = mw_div(r,&a);
    tmp1 = sqr(&tmp1);
    tmp1 = mw_add_s(&tmp1, 1.0);
    tmp1 = minusfivehalves(&tmp1);

    tmp2 = cube(&a);
    tmp2 = mw_div(&M, &tmp2);
    tmp2 = mw_mul_s(&tmp2, 3.0/(4.0*M_PI));

    return mw_mul(&tmp2, &tmp1);

  /*real r_a = r/a;
    real rho_peak = 3*M/(4.0*pi*mw_pow(a, 3.0));
    return rho_peak*mw_pow((1.0+mw_pow(r_a,2.0)), -2.5);*/
}

static inline real NFWMHaloDensity(const Halo* h,  real* r)
{
    if(showRealValue(r) == 0) return ZERO_REAL;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    real tmp1, tmp2;

    tmp1 = mw_div(&M, r);
    tmp2 = mw_add(&a,r);
    tmp2 = sqr(&tmp2);
    tmp1 = mw_div(&tmp1, &tmp2);
    
    return mw_mul_s(&tmp1, inv_0(4.0*M_PI));

  /*return M / (r*(a+r)*(a+r)*4.0*pi);*/

}

static inline real allenSantillanHaloDensity(const Halo* h, real* r)
{
    if(showRealValue(r)==0.0) return ZERO_REAL;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real_0 lam = h->lambda;
    const real gam = mw_real_const(h->gamma);
    if(showRealValue(r) > lam) return ZERO_REAL;
    real tmp1, tmp2;

    real b = mw_add_s(&gam, -1.0);

    tmp1 = mw_div(r,&a);
    tmp1 = mw_pow(&tmp1,&b);
    tmp2 = mw_div(r,&a);
    tmp2 = mw_pow(&tmp2,&b);
    tmp2 = mw_add(&tmp2, &gam);
    real numer = mw_mul(&tmp1, &tmp2);

    tmp1 = sqr(r);
    tmp2 = mw_div(r,&a);
    tmp2 = mw_pow(&tmp2,&b);
    tmp2 = mw_add_s(&tmp2, 1.0);
    tmp2 = sqr(&tmp2);
    real denom = mw_mul(&tmp1,&tmp2);

    tmp1 = mw_mul(&M,&numer);
    tmp1 = mw_div(&tmp1,&a);
    tmp1 = mw_div(&tmp1,&denom);

    return mw_mul_s(&tmp1, inv_0(4.0*M_PI));

  /*real b = gam - 1.0;
    real numer = mw_pow(r/a,b)*(mw_pow(r/a,b) + b + 1.0);
    real denom = r*r*mw_pow(1.0 + mw_pow(r/a,b),2.0);

    return M*numer/a/denom/4.0/pi;*/
}

static inline real wilkinsonEvansHaloDensity(const Halo* h, real* r)
{
    if(showRealValue(r)==0.0) return ZERO_REAL;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    real tmp1, tmp2;

    tmp1 = mw_div(r,&a);
    tmp1 = sqr(&tmp1);
    tmp2 = mw_hypot(r,&a);
    tmp2 = cube(&tmp2);
    tmp1 = mw_mul(&tmp1, &tmp2);
    tmp1 = mw_div(&M, &tmp1);

    return mw_mul_s(&tmp1, inv_0(4*M_PI));
  /*return (1/(4*pi)) * M*a*a/(r*r*mw_pow(r*r + a*a,1.5));*/
}

static inline real ninkovicHaloDensity(const Halo* h, real* r)
{
    const real rho = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    const real lam = mw_real_const(h->lambda);
    if(showRealValue(r) > showRealValue(&lam))  return ZERO_REAL; 
    real tmp1, tmp2;

    tmp1 = mw_div(r,&a);
    tmp1 = cube(&tmp1);
    tmp1 = mw_add_s(&tmp1, 1.0);
    tmp1 = inv(&tmp1);
    tmp2 = mw_div(&lam,&a);
    tmp2 = cube(&tmp2);
    tmp2 = mw_add_s(&tmp2, 1.0);
    tmp2 = inv(&tmp2);
    tmp1 = mw_sub(&tmp1, &tmp2);

    return mw_mul(&rho,&tmp1);
  /*return rho*(1.0/(1.0 + mw_pow(r/a,3.0)) - 1.0/(1.0 + mw_pow(lam/a,3.0)));*/
}

static inline real KVHalo(const Halo* h, real* r) /*What is this one?*/
{
    if(showRealValue(r) == 0) return ZERO_REAL;
    const real M = mw_real_var(h->mass, HALO_MASS_POS);
    const real a = mw_real_var(h->scaleLength, HALO_RADIUS_POS);
    real tmp1, tmp2;

    tmp1 = mw_add(r,&a);
    tmp1 = sqr(&tmp1);
    tmp1 = mw_mul(r, &tmp1);
    tmp1 = inv(&tmp1);
    tmp2 = mw_add(r,&a);
    tmp2 = cube(&tmp2);
    tmp2 = inv(&tmp2);
    tmp2 = mw_mul_s(&tmp2, 2.0);
    tmp1 = mw_sub(&tmp1, &tmp2);
    tmp1 = mw_mul(&M, &tmp1);

    return mw_mul_s(&tmp1, inv_0(4*M_PI));
  /* return (1/(4*pi)) * (M/(r*mw_pow(r+a, 2))) - ((2*M)/(mw_pow(r+a, 3)));*/
}

real nbExtDensity(const Potential* pot, mwvector* pos, real_0 time)
{
    //mw_printf("POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(pos)), showRealValue(&Y(pos)), showRealValue(&Z(pos)) );
    real density = ZERO_REAL;
    const real_0 limit_val = mw_pow_0(2.0,-8.0);
    real r = mw_absv(pos);
    real limit = r;
    setRealValue(&limit, limit_val);
    real bulge_den;
    real disk_den;
    real disk2_den;
    real halo_den;

    /* Change r if less than limit. Done this way to pipeline this step*/
    real check1 = mw_mul_s(&r, (showRealValue(&r)>limit_val));
    real check2 = mw_mul_s(&limit, (showRealValue(&r)<=limit_val));
    r = mw_add(&check1, &check2);
    //mw_printf("r = %.15f\n", showRealValue(&r));

    switch (pot->sphere[0].type)
    {
        case HernquistSpherical:
            bulge_den = hernquistSphericalDensity(&(pot->sphere[0]), &r);
            break;
        case PlummerSpherical:
            bulge_den = plummerSphericalDensity(&(pot->sphere[0]), &r);
            break;
        case NoSpherical:
            bulge_den = ZERO_REAL;
            break;
        case InvalidSpherical:
        default:
            mw_fail("Invalid bulge type in density\n");
    }
    //mw_printf("Bulge Density = %.15f\n", showRealValue(&bulge_den));
    density = mw_add(&density, &bulge_den);

    switch (pot->disk.type)
    {
        case FreemanDisk: /*Density negligible since infinitely thin*/
            break;
        case MiyamotoNagaiDisk:
            disk_den = miyamotoNagaiDiskDensity(&(pot->disk), pos);
            break;
        case DoubleExponentialDisk:
            disk_den = doubleExponentialDiskDensity(&(pot->disk), pos);
            break;
        case Sech2ExponentialDisk:
            disk_den = sech2ExponentialDiskDensity(&(pot->disk), pos);
            break;
        case OrbitingBar:
            disk_den = orbitingBarDensity(&(pot->disk), pos, time);
            break;
        case NoDisk:
            disk_den = ZERO_REAL;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid primary disk type in density\n");
    }
    //mw_printf("Disk Density  = %.15f\n", showRealValue(&disk_den));
    density = mw_add(&density, &disk_den);

    switch (pot->disk2.type)
    {
        case FreemanDisk: /*Density negligible since infinitely thin*/
            break;
        case MiyamotoNagaiDisk:
            disk2_den = miyamotoNagaiDiskDensity(&(pot->disk2), pos);
            break;
        case DoubleExponentialDisk:
            disk2_den = doubleExponentialDiskDensity(&(pot->disk2), pos);
            break;
        case Sech2ExponentialDisk:
            disk2_den = sech2ExponentialDiskDensity(&(pot->disk2), pos);
            break;
        case OrbitingBar:
            disk2_den = orbitingBarDensity(&(pot->disk2), pos, time);
            break;
        case NoDisk:
            disk2_den = ZERO_REAL;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid primary disk type in density\n");
    }
    //mw_printf("Disk2 Density = %.15f\n", showRealValue(&disk2_den));
    density = mw_add(&density, &disk2_den);

    switch (pot->halo.type)
    {
        case LogarithmicHalo:
            halo_den = logarithmicHaloDensity(&(pot->halo), pos);
            break;
        case NFWHalo:
            halo_den = NFWHaloDensity(&(pot->halo), &r);
            break;
        case TriaxialHalo:
            halo_den = triaxialHaloDensity(&(pot->halo), pos);
            break;
        case CausticHalo: /*FIXME: Add density profile for caustic halo when we actually plan on making this work*/
            halo_den = ZERO_REAL;
            break;
        case AllenSantillanHalo:
            halo_den = allenSantillanHaloDensity(&(pot->halo), &r);
            break;
        case WilkinsonEvansHalo:
            halo_den = wilkinsonEvansHaloDensity(&(pot->halo), &r);
	    break;
        case NFWMassHalo:
            halo_den = NFWMHaloDensity(&(pot->halo), &r);
            break;
        case PlummerHalo:
            halo_den = plummerHaloDensity(&(pot->halo), &r);
            break;
        case HernquistHalo:
            halo_den = hernquistHaloDensity(&(pot->halo), &r);
            break;
        case NinkovicHalo:
            halo_den = ninkovicHaloDensity(&(pot->halo), &r);
            break;
        case NoHalo:
            halo_den = ZERO_REAL;
            break;
        case InvalidHalo:
        default:
            mw_fail("Invalid halo type in density\n");
    }
    //mw_printf("Halo Density  = %.15f\n", showRealValue(&halo_den));
    density = mw_add(&density, &halo_den);

    if (showRealValue(&density) < 0.0)
    {
        mw_fail("Negative density calculated!\n    Faulty Potential = %s\n", showPotential(pot));
    }

    return density;

}
