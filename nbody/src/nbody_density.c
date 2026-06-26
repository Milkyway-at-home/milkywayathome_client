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

static const real pi = 3.1415926535;

/*Spherical Buldge Densities*/
static inline real hernquistSphericalDensity(const Spherical* sph, real r)
{
    const real a = sph->scale;

    /*return 0 rather than get a divide by 0 error*/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(r*mw_pow(r+a, 3) == 0) {
#pragma GCC diagnostic pop
        return 0;
    }

    return ((sph->mass)/(2*pi*mw_pow(a, 3)))*mw_pow(a, 4)/(r*mw_pow(r+a, 3));
}

static inline real plummerSphericalDensity(const Spherical* sph, real r)
{
    const real a = sph->scale;
    const real M = sph->mass;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(a == 0) {
#pragma GCC diagnostic pop
        return 0;
    }
    real r_a = r/a;

    real rho_peak = 3*M/(4.0*pi*mw_pow(a, 3.0));

    return rho_peak*mw_pow((1.0+mw_pow(r_a,2.0)), -2.5);
}

/*Disk Densities*/
static inline real miyamotoNagaiDiskDensity(const Disk* disk, mwvector pos)
{
    const real a   = disk->scaleLength;
    const real b   = disk->scaleHeight;
    const real M   = disk->mass;
    const real R   = mw_pow(mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0), 0.5);
    const real zp  = mw_pow(mw_pow(Z(pos),2.0) + mw_pow(b,2.0), 0.5);
    const real azp = a + zp;

    real numer = M*b*b*(a*R*R + (a + 3.0*zp) * mw_pow(azp, 2));
    real denom = 4.0*pi*mw_pow(R*R + azp*azp, 2.5)*mw_pow(zp, 3.0);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(denom == 0) return 0;
#pragma GCC diagnostic pop

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
static inline real logarithmicHaloDensity(const Halo* h, mwvector pos) /** flattenZ should be greater than 1/sqrt(2) to keep positive definite **/
{
    const real v  = h->vhalo;
    const real q  = h->flattenZ;
    const real a  = h->scaleLength;
    const real R2 = mw_pow(X(pos),2.0) + mw_pow(Y(pos),2.0);

    real numer = (2.0*q*q + 1.0)*a*a + R2 + mw_pow(Z(pos),2.0)*(2.0-1.0/q/q);
    real denom = q*q*mw_pow(R2 + a*a + mw_pow(Z(pos)/q,2.0),2.0);

    return v*v*numer/2.0/pi/denom;
}

static inline real SphericalNFWerkalHaloDensity(const Halo* h, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;

    const real c = mw_log(1.0 + 15.3) - (15.3/(1.0 + 15.3));

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(r == 0) return 0;
#pragma GCC diagnostic pop

    return M / (c*r*(a+r)*(a+r)*4.0*pi);

}

static inline real NFWHaloDensity(const Halo* h,  real r)
{
    const real a = h->scaleLength;
    const real v = h->vhalo;

    real rho = v*v/4.0/pi/a/a/0.2162165954;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(r == 0) return 0;
#pragma GCC diagnostic pop
    
    return rho / (r/a) / mw_pow(1.0+(r/a),2.0);

}

static inline real triaxialHaloDensity(const Halo* h, mwvector pos)
{
    const real v   = h->vhalo;
    const real q   = h->flattenZ;
    const real a   = h->scaleLength;

    const real D   = a*a + (h->c1)*mw_pow(X(pos),2.0) + (h->c2)*mw_pow(Y(pos),2.0) + (h->c3)*X(pos)*Y(pos) + mw_pow(Z(pos),2.0)/q/q;
    const real num = 2.0*(h->c1)*D + 2.0*(h->c2)*D + 2.0*D/q/q - mw_pow(2.0*(h->c1)*X(pos) + (h->c3)*Y(pos),2.0) - mw_pow(2.0*(h->c2)*Y(pos) + (h->c3)*X(pos),2.0) - mw_pow(2.0*Z(pos)/q/q,2.0);

    return v*v*num/4.0/pi/D/D;
}

static inline real orbitingBarDensity(const Disk* disk, mwvector pos, real time)
{
    // Bar parameters
    const real a = disk->scaleLength; // Bar half-length
    const real b = 1.4;               // Triaxial softening length
    const real c = 1.0;               // Prolate softening length

    // Calculate current bar angle
    const real curAngle = disk->startAngle - disk->patternSpeed * time;

    // Calculate cylindrical radius and angle using atan2 for robustness
    const real px = pos.x;
    const real py = pos.y;
    const real pz = pos.z;
    const real Radi = mw_sqrt(px * px + py * py);
    const real Phi = mw_atan2(py, px) - curAngle;

    // Rotate position into bar frame
    const real x = Radi * mw_cos(Phi);
    const real y = Radi * mw_sin(Phi);
    const real z = pz;

    // Intermediate calculations for density
    const real z2 = z * z;
    const real c2 = c * c;
    const real zc = mw_sqrt(z2 + c2);
    const real bzc = b + zc;
    const real bzc2 = bzc * bzc;
    const real y2 = y * y;
    const real bigA = b * y2 + (b + 3.0 * zc) * bzc2;
    const real bigC = y2 + bzc2;

    // Precompute powers and denominators
    const real bigC2 = bigC * bigC;
    const real zc3 = zc * zc * zc;
    const real denom = 24.0 * pi * a * bigC2 * zc3;
    const real c2_over_denom = c2 / denom;

    // Compute terms for (x+a) and (x-a)
    const real xp = x + a;
    const real xm = x - a;
    const real xp2 = xp * xp;
    const real xm2 = xm * xm;
    const real bigCp = bigC + xp2;
    const real bigCm = bigC + xm2;
    const real bigCp15 = mw_pow(bigCp, 1.5);
    const real bigCm15 = mw_pow(bigCm, 1.5);

    const real termp = (xp * (3.0 * bigA * bigC + (2.0 * bigA + b * bigC) * xp2)) / bigCp15;
    const real termm = (xm * (3.0 * bigA * bigC + (2.0 * bigA + b * bigC) * xm2)) / bigCm15;

    const real unscaledDens = c2_over_denom * (termp - termm);

    return unscaledDens * disk->mass;
}

/*Halo Densities*/
static inline real hernquistHaloDensity(const Halo* h,  real r)
{
    const real a = h->scaleLength;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(r*mw_pow(r+a, 3) == 0 || 2*pi*mw_pow(a, 3) == 0) return 0;
#pragma GCC diagnostic pop

    return ((h->mass)/(2*pi*mw_pow(a, 3)))*mw_pow(a, 4)/(r*mw_pow(r+a, 3));
}

static inline real plummerHaloDensity(const Halo* h, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(a == 0) {
#pragma GCC diagnostic pop
        return 0;
    }
    real r_a = r/a;

    real rho_peak = 3*M/(4.0*pi*mw_pow(a, 3.0));

    return rho_peak*mw_pow((1.0+mw_pow(r_a,2.0)), -2.5);
}

static inline real NFWMHaloDensity(const Halo* h,  real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(r == 0) return 0;
#pragma GCC diagnostic pop
    
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(r==0.0) return 0.0;
#pragma GCC diagnostic pop

    if(r > lam) return 0.0;

    return M*numer/a/denom/4.0/pi;
}

static inline real wilkinsonEvansHaloDensity(const Halo* h, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(r==0.0) return 0.0;
#pragma GCC diagnostic pop

    return (1/(4*pi)) * M*a*a/(r*r*mw_pow(r*r + a*a,1.5));

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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    if(r + a == 0) return 0;
#pragma GCC diagnostic pop

    return (1/(4*pi)) * (M/(r*mw_pow(r+a, 2))) - ((2*M)/(mw_pow(r+a, 3)));
}

real nbExtDensity(const Potential* pot, mwvector pos, real time)
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
        case OrbitingBar:
            density += orbitingBarDensity(&(pot->disk), pos, time);
            break;
        case NoDisk:
            density += 0.0;
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid primary disk type in density\n");
            break;
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
        case OrbitingBar:
            density += orbitingBarDensity(&(pot->disk2), pos, time);
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
        case SphericalNFWerkalHalo:
	        density += SphericalNFWerkalHaloDensity(&(pot->halo), r);
	        break;
        case NoHalo:
            density += 0.0;
            break;
        case InvalidHalo:
        default:
            mw_fail("Invalid halo type in density\n");
    }

    if (density < 0.0)
    {
        mw_fail("Negative density calculated!\n    Faulty Potential = %s\n", showPotential(pot));
    }

    return density;

}

 
