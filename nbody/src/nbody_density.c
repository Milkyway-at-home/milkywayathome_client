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

static inline real plummerSpherical(const Spherical* sph, real r)
{
    const real a = sph->scale;
    if((4*pi*mw_pow(a, 3))*mw_pow((1+(r/a)*(r/a)), -5/2) == 0) return 0;

    return ((3*sph->mass)/(4*pi*mw_pow(a, 3))*mw_pow((1+(r/a)*(r/a)), -5/2));
}

/*Disk Densities*/
static inline real miyamotoNagaiDiskAccel(const Disk* disk, mwvector pos, real r)
{
    const real a   = disk->scaleLength;
    const real b   = disk->scaleHeight;
    const real zp  = mw_sqrt(sqr(Z(pos)) + sqr(b));
    

    real numer = disk->mass*b*b*(a*r*r + (a + 3*mw_pow(zp*zp+b*b, 1/2) * mw_pow(a+mw_pow(zp*zp+b*b, 1/2), 2)));
    real denom = 4*pi*mw_pow(r*r + mw_pow(a+mw_pow(zp*zp+b*b, 1/2), 2), 5/2)*mw_pow(zp*zp+b*b, 3/2);

    if(denom == 0) return 0;

    return numer/denom;
}

/*Halo Densities*/
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

    if(a == 0 || ((r/a)*(1 + r/a)*(1 + r/a)) == 0) return 0;
    
    return M / ((r/a)*(1 + r/a)*(1 + r/a)); /*the M in the numerator should be a scale density*/

}

static inline real wilkinsonHalo(const Halo* h, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;

    if((r*r*mw_pow(r*r + a*a, 1/2)) == 0)  return 0; 

    return (1/(4*pi)) * ((M/(r*r*mw_pow(r*r + a*a, 1/2))) - ((M * (2*r*r + a*a))/(r*r*mw_pow(r*r + a*a, 3/2))));

}

static inline real KVHalo(const Halo* h, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;
    if(r + a == 0) return 0;

    return (1/(4*pi)) * (M/(r*mw_pow(r+a, 2))) - ((2*M)/(mw_pow(r+a, 3)));
}
 
