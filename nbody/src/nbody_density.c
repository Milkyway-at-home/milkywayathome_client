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

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

const real pi = 3.1415926535;

/*Methods to be called by densities*/

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
    real val = sqr(k)*zd/mw_sinh(3.1415926535*k*zd/2.0)*mw_cos(k*z)*(aExp(k,R,Rd)*besselK1(k*R) - bExp(k,R,Rd)*besselI1(k*R));
    //mw_printf("RSech(%.15f,%.15f,%.15f,%.15f,%.15f) = %.15f\n",k,R,Rd,z,zd,val);
    return val;
}

static inline real ZSechIntegrand (real k, real R, real Rd, real z, real zd)
{
    real val = sqr(k)*zd/mw_sinh(3.1415926535*k*zd/2.0)*mw_sin(k*z)*(aExp(k,R,Rd)*besselK0(k*R) + bExp(k,R,Rd)*besselI0(k*R));
    //mw_printf("ZSech = %.15f\n",val);
    return val;
}


/*************************************************************************************************************************************************************************************************************************************/
/*Spherical Buldge Densities*/
static inline real hernquistSphericalDensity(const Spherical* sph, real r)
{
    const real a = sph->scale;

    /*return 0 rather than get a divide by 0 error*/
    if(r*mw_pow(r+a, 3) == 0) {
        return 0;
    }

    return ((sph->mass)/(2*pi*mw_pow(a, 3)))*mw_pow(a, 4)/(r*mw_pow(r+a, 3));
}

static inline real plummerSphericalDensity(const Spherical* sph, real r)
{
    const real a = sph->scale;
    if((4*pi*mw_pow(a, 3))*mw_pow((1+(r/a)*(r/a)), -5/2) == 0) return 0;

    return ((3*sph->mass)/(4*pi*mw_pow(a, 3))*mw_pow((1+(r/a)*(r/a)), -5/2));
}

/*Disk Densities*/
static inline real miyamotoNagaiDiskDensity(const Disk* disk, mwvector pos)
{
    const real a   = disk->scaleLength;
    const real b   = disk->scaleHeight;
    const real zp  = mw_sqrt(sqr(Z(pos)) + sqr(b));
    const real R   = mw_sqrt(mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0));
    

    real numer = disk->mass*b*b*(a*R*R + (a + 3.0*mw_pow(zp*zp+b*b, 0.5) * mw_pow(a+mw_pow(zp*zp+b*b, 0.5), 2)));
    real denom = 4.0*pi*mw_pow(R*R + mw_pow(a+mw_pow(zp*zp+b*b, 0.5), 2.0), 2.5)*mw_pow(zp*zp+b*b, 1.5);

    if(denom == 0) return 0;

    return numer/denom;
}

static inline real doubleExponentialDiskDensity(const Disk* disk, mwvector pos)
{
    const real M   = disk->mass;
    const real d_r = disk->scaleLength;
    const real d_z = disk->scaleHeight;
    const real R   = mw_sqrt(mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0));

    real den = M/(4.0*pi*d_z*mw_pow(d_r,2.0))*mw_exp(-R/d_r)*mw_exp(-mw_abs(Z(pos))/d_z);

    return den;

}

static inline real sech2ExponentialDiskDensity(const Disk* disk, mwvector pos)
{
    const real M   = disk->mass;
    const real d_r = disk->scaleLength;
    const real d_z = disk->scaleHeight;
    const real R   = mw_sqrt(mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0));

    real den = M/(4.0*pi*d_z*d_r*d_r)*mw_exp(-R/d_r)/mw_pow(mw_cosh(Z(pos)/d_z),2.0);

    return den;

}

/*Halo Densities*/
static inline real logarithmicHaloDensity(const Halo* h, mwvector pos)
{
    const real v  = h->vhalo;
    const real q  = h->flattenZ;
    const real a  = h->scaleLength;
    const real R2 = mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0);

    real numer = (2.0*q*q + 1.0)*a*a + R2 + mw_pow(Z(pos),2.0)*(2.0-1.0/q/q);
    real denom = mw_pow(R2+a*a+mw_pow(Z(pos),2.0)/q/q,2.0);

    return v*v*numer/2.0/pi/denom;
}

static inline real triaxialHaloDensity(const Halo* h, mwvector pos)
{
    const real v   = h->vhalo;
    const real q   = h->flattenZ;
    const real a   = h->scaleLength;

    const real D   = a*a + (h->c1)*mw_pow(X(pos),2.0) + (h->c2)*mw_pow(Y(pos),2.0) + (h->c3)*X(pos)*Y(pos) + mw_pow(Z(pos),2.0)/q/q;
    const real num = 2.0*(h->c1)*D + 2.0*(h->c2)*D + 2.0*Z(pos)/q/q - mw_pow(2.0*(h->c1)*X(pos) + (h->c3)*Y(pos),2.0) - mw_pow(2.0*(h->c2)*Y(pos) + (h->c3)*X(pos),2.0) - mw_pow(2.0*Z(pos)/q/q,2.0);

    return v*v*num/4.0/pi/D/D;
}

static inline real hernquistHaloDensity(const Halo* h,  real r)
{
    const real a = h->scaleLength;

    if(r*mw_pow(r+a, 3) == 0 || 2*pi*mw_pow(a, 3) == 0) return 0;

    return ((h->mass)/(2*pi*mw_pow(a, 3)))*mw_pow(a, 4)/(r*mw_pow(r+a, 3));
}

static inline real plummerHaloDensity(const Halo* h, real r)
{
    const real a = h->scaleLength;
    if(a == 0 || 4*pi*mw_pow(a, 3)*mw_pow((1+(r/a)*(r/a)), -5/2) == 0) return 0;

    return ((3*h->mass)/((4*pi*mw_pow(a, 3))*mw_pow((1+(r/a)*(r/a)), -5/2)));
}

static inline real NFWMHaloDensity(const Halo* h,  real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;

    if(r == 0) return 0;
    
    return M / (r*(a+r)*(a+r)*4.0*pi);

}

static inline real allenSantillanHaloDensity(const Halo* h, real r)
{
    const real a   = h->scaleLength;
    const real M   = h->mass;
    const real lam = h->lambda;
    const real gam = h->gamma;

    real b = gam - 1.0;
    real numer = mw_pow(r/a,b)*(mw_pow(r/a,b) + b + 1.0);
    real denom = r*r*mw_pow(1.0 + mw_pow(r/a,b),2.0);

    if(r==0.0) return 0.0;

    if(r > lam) return 0.0;

    return M*numer/a/denom;
}

static inline real wilkinsonEvansHaloDensity(const Halo* h, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;

    if((r*r*mw_pow(r*r + a*a, 1/2)) == 0)  return 0; 

    return (1/(4*pi)) * ((M/(r*r*mw_pow(r*r + a*a, 1/2))) - ((M * (2*r*r + a*a))/(r*r*mw_pow(r*r + a*a, 3/2))));

}

static inline real ninkovicHaloDensity(const Halo* h, real r)
{
    const real a = h->scaleLength;
    const real rho = h->rho0;
    const real lam = h->lambda;

    const real density = 1.0/(1.0 + mw_pow(r/a,3.0)) - 1.0/(1.0 + mw_pow(lambda/a,3.0))

    if(density < 0)  return 0.0; 

    return density;
}

static inline real KVHalo(const Halo* h, real r) /*What is this one?*/
{
    const real a = h->scaleLength;
    const real M = h->mass;
    if(r + a == 0) return 0;

    return (1/(4*pi)) * (M/(r*mw_pow(r+a, 2))) - ((2*M)/(mw_pow(r+a, 3)));
}
 
