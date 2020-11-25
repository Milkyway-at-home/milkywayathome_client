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

static const real pi = 3.1415926535;

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

static inline real NFWHaloDensity(const Halo* h,  real r)
{
    const real a = h->scaleLength;
    const real v = h->vhalo;

    real rho = v*v/4.0/pi/a/a/0.2162165954;

    if(r == 0) return 0;
    
    return rho / (r/a) / mw_pow(1.0+(r/a),2.0);

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

    const real density = rho*(1.0/(1.0 + mw_pow(r/a,3.0)) - 1.0/(1.0 + mw_pow(lam/a,3.0)));

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

real nbExtDensity(const Potential* pot, mwvector pos)
{
    real density = 0.0;
    const real limit = mw_pow(2.0,-8.0);

    /* Change r if less than limit. Done this way to pipeline this step*/
    real r = (mw_absv(pos) <= limit)*limit + (mw_absv(pos) > limit)*mw_absv(pos);

    switch (pot->sphere[0].type)
    {
        case HernquistSpherical:
            density += hernquistSphericalDensity(&(pot->sphere[0]), r);
            break;
        case PlummerSpherical:
            density += plummerSphericalDensity(&(pot->sphere[0]), r);
            break;
        case NoSpherical:
            density += 0.0;
            break;
        case InvalidSpherical:
        default:
            mw_fail("Invalid bulge type in density\n");
    }

    switch (pot->disk.type)
    {
        case FreemanDisk:
            density += 0.0; /*Density negligible since infinitely thin*/
            break;
        case MiyamotoNagaiDisk:
            density += miyamotoNagaiDiskDensity(&(pot->disk), pos);
            break;
        case DoubleExponentialDisk:
            density += doubleExponentialDiskDensity(&(pot->disk), pos);
            break;
        case Sech2ExponentialDisk:
            density += sech2ExponentialDiskDensity(&(pot->disk), pos);
            break;
        case NoDisk:
            density += 0.0;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid primary disk type in density\n");
    }

    switch (pot->disk2.type)
    {
        case FreemanDisk:
            density += 0.0; /*Density negligible since infinitely thin*/
            break;
        case MiyamotoNagaiDisk:
            density += miyamotoNagaiDiskDensity(&(pot->disk2), pos);
            break;
        case DoubleExponentialDisk:
            density += doubleExponentialDiskDensity(&(pot->disk2), pos);
            break;
        case Sech2ExponentialDisk:
            density += sech2ExponentialDiskDensity(&(pot->disk2), pos);
            break;
        case NoDisk:
            density += 0.0;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid secondary disk type in density\n");
    }

    switch (pot->halo.type)
    {
        case LogarithmicHalo:
            density += logarithmicHaloDensity(&(pot->halo), pos);
            break;
        case NFWHalo:
            density += NFWHaloDensity(&(pot->halo), r);
            break;
        case TriaxialHalo:
            density += triaxialHaloDensity(&(pot->halo), pos);
            break;
        case CausticHalo:
            density += 0.0; /*FIXME: Add density profile for caustic halo when we actually plan on making this work*/
            break;
        case AllenSantillanHalo:
            density += allenSantillanHaloDensity(&(pot->halo), r);
            break;
        case WilkinsonEvansHalo:
            density += wilkinsonEvansHaloDensity(&(pot->halo), r);
	    break;
        case NFWMassHalo:
            density += NFWMHaloDensity(&(pot->halo), r);
            break;
        case PlummerHalo:
            density += plummerHaloDensity(&(pot->halo), r);
            break;
        case HernquistHalo:
            density += hernquistHaloDensity(&(pot->halo), r);
            break;
        case NinkovicHalo:
            density += ninkovicHaloDensity(&(pot->halo), r);
            break;
        case NoHalo:
            density += 0.0;
            break;
        case InvalidHalo:
        default:
            mw_fail("Invalid halo type in density\n");
    }

    return density;

}

 
