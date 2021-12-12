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

    /*When initializing a real parameter that we want to differentiate over,
      we must specify which rows and columns of the gradient and hessian attributes
      carry the derivatives with respect to that parameter. These are the parameters
      we currently differentiate over and their assigned column number n.

      +----PARAMETER---------------------------n----------+
           Backwards Evolution Time            0
           Time Ratio                          1
           Baryonic Plummer Radius             2
           Radius Ratio                        3
           Baryonic Mass                       4
           Mass Ratio                          5
           Orbital b coord                     6
           Orbital r coord                     7
           Orbital vx coord                    8
           Orbital vy coord                    9
           Orbital vz coord                   10
           MW Bulge Mass                      11
           MW Bulge Radius                    12
           MW Disk Mass                       13
           MW Disk Scale Length               14
           MW Disk Scale Height               15
           MW Halo Mass (or Scale Velocity)   16
           MW Halo Scale Radius               17
           MW Halo Flattening Parameter (Z)   18
           LMC Mass                           19
           LMC Scale Radius                   20

      FIXME: Most of our Halo models use more parameters that we are not listed here.
      As such, only a Logarithmic Halo and a single Miyamoto-Nagai disk can be used
      with AUTODIFF. To use other models, you will need to either reassign some numbers
      or add more parameters to the list. The number of model parameters to differentiate
      over can be found in milkyway_math_autodiff.h*/

/*Spherical Buldge Densities*/
static inline real hernquistSphericalDensity(const Spherical* sph, real r)
{
    const real a = mw_real_var(sph->scale, 12);
    const real M = mw_real_var(sph->mass, 11);

    /*return 0 rather than get a divide by 0 error*/
    if(showRealValue(r) == 0) {
        return ZERO_REAL;
    }

    return mw_mul_s(mw_div(mw_mul(M,a), mw_mul(r,cube(mw_add(r,a)))), inv_0(2*M_PI));
}

static inline real plummerSphericalDensity(const Spherical* sph, real r)
{
    const real a = mw_real_var(sph->scale, 12);
    const real M = mw_real_var(sph->mass, 11);
    if(showRealValue(a) == 0)
    {
        return ZERO_REAL;
    }
    real r_a = mw_div(r,a);

    real rho_peak = mw_mul_s(mw_div(M,cube(a)), 3.0/(4.0*M_PI));

    return mw_mul(rho_peak, minusfivehalves(mw_add(mw_real_const(1.0), sqr(r_a))));
}

/*Disk Densities*/
static inline real miyamotoNagaiDiskDensity(const Disk* disk, mwvector pos)
{
    const real M   = mw_real_var(disk->mass, 13);
    const real a   = mw_real_var(disk->scaleLength, 14);
    const real b   = mw_real_var(disk->scaleHeight, 15);

    const real R   = mw_hypot(X(pos), Y(pos));
    const real zp  = mw_hypot(Z(pos), b);
    const real azp = mw_add(a, zp);

    real numer = mw_mul(M, mw_mul(sqr(b),mw_add(mw_mul(a,sqr(R)), mw_mul(mw_add(a, mw_mul_s(zp,3.0)),sqr(azp)))));
    real denom = mw_mul_s(mw_mul(fifth(mw_hypot(R, azp)),cube(zp)), 4.0*M_PI);

    if(showRealValue(denom) == 0) return ZERO_REAL;

    return mw_div(numer, denom);
}

static inline real doubleExponentialDiskDensity(const Disk* disk, mwvector pos)
{
    const real M   = mw_real_var(disk->mass, 13);
    const real d_r = mw_real_var(disk->scaleLength, 14);
    const real d_z = mw_real_var(disk->scaleHeight, 15);
    const real R   = mw_hypot(X(pos), Y(pos));

    real den = mw_mul_s(mw_div(mw_div(M, mw_mul(d_z,sqr(d_r))), mw_exp(mw_add(mw_div(R,d_r), mw_div(mw_abs(Z(pos)),d_z)))), inv_0(4.0*M_PI));

    return den;

}

static inline real sech2ExponentialDiskDensity(const Disk* disk, mwvector pos)
{
    const real M   = mw_real_var(disk->mass, 13);
    const real d_r = mw_real_var(disk->scaleLength, 14);
    const real d_z = mw_real_var(disk->scaleHeight, 15);
    const real R   = mw_hypot(X(pos), Y(pos));

    real den = mw_mul_s(mw_div(mw_div(mw_div(M, mw_mul(d_z, sqr(d_r))), mw_exp(mw_div(R, d_r))), sqr(mw_cosh(mw_div(Z(pos), d_z)))), inv_0(4.0*M_PI));

    return den;

}

static inline real orbitingBarDensity(const Disk* disk, mwvector pos, real_0 time)
{
    real M = mw_real_var(disk->mass, 13);
    real a = mw_real_var(disk->scaleLength, 14);  // Bar half-length
    real b = mw_real_const(1.4);                  //Triaxial softening length
    real c = mw_real_const(1.0);                  //Prolate softening length

    real pSpeed = mw_real_var(disk->patternSpeed, 1); //Several of these parameters are not assigned in AUTODIFF, so we use 1 as a placeholder//
    real sAngle = mw_real_var(disk->startAngle, 1);

    real curAngle =  mw_sub(sAngle, mw_mul_s(pSpeed, time));
    //first rotate pos curAngle * -1 radians to emulate the current angle of the bar
    real Radi = mw_hypot(X(pos), Y(pos));
    real Phi = mw_atan2(Y(pos), X(pos));
    Phi = mw_sub(Phi, curAngle);

    real x = mw_mul(Radi, mw_cos(Phi));
    real y = mw_mul(Radi, mw_sin(Phi)); 
    real z = Z(pos);

    real zc = mw_hypot(z,c);
    real bzc2 = sqr(mw_add(b,zc));
    real bigA = mw_add(mw_mul(b,sqr(y)), mw_mul(mw_add(b, mw_mul_s(zc,3.0)), bzc2));
    real bigC = mw_add(sqr(y), bzc2);

    real AC3 = mw_mul_s(mw_mul(bigA,bigC), 3.0);
    real A2bC = mw_add(mw_mul_s(bigA,2.0), mw_mul(b,bigC));

    real factor = mw_mul_s(mw_div(mw_div(mw_div(sqr(c),a),sqr(bigC)),cube(zc)),inv_0(24*M_PI));
    real part1 = mw_div(mw_mul(mw_add(x,a), mw_add(AC3, mw_mul(A2bC, sqr(mw_add(x,a))))),threehalves(mw_add(bigC,sqr(mw_add(x,a)))));
    real part2 = mw_div(mw_mul(mw_sub(x,a), mw_add(AC3, mw_mul(A2bC, sqr(mw_sub(x,a))))),threehalves(mw_add(bigC,sqr(mw_sub(x,a)))));

    real unscaledDens = mw_mul(factor, mw_sub(part1, part2));

    return mw_mul(unscaledDens, M);
}

/*Halo Densities*/
static inline real logarithmicHaloDensity(const Halo* h, mwvector pos) /** flattenZ should be greater than 1/sqrt(2) to keep positive definite **/
{
    const real v  = mw_real_var(h->vhalo, 16);
    const real a  = mw_real_var(h->scaleLength, 17);
    const real q  = mw_real_var(h->flattenZ, 18);

    const real R2 = mw_add(sqr(X(pos)), sqr(Y(pos)));

    real numer = mw_add(mw_add(mw_mul(mw_add(mw_mul_s(sqr(q),2.0), mw_real_const(1.0)), sqr(a)), R2), mw_mul(sqr(Z(pos)),mw_sub(mw_real_const(2.0), inv(sqr(q)))));
    real denom = mw_mul(sqr(q),sqr(mw_add(mw_add(R2, sqr(a)), sqr(mw_div(Z(pos), q)))));

    return mw_mul_s(mw_div(mw_mul(sqr(v),numer),denom), inv_0(2.0*M_PI));
}

static inline real NFWHaloDensity(const Halo* h,  real r)
{
    const real v  = mw_real_var(h->vhalo, 16);
    const real a  = mw_real_var(h->scaleLength, 17);

    real rho = mw_mul_s(sqr(mw_div(v,a)), inv_0(4.0*M_PI*0.2162165954));

    if(showRealValue(r) == 0) return ZERO_REAL;
    
    return mw_div(mw_div(rho, mw_div(v,a)), sqr(mw_add(mw_real_const(1.0), mw_div(v,a))));

}

static inline real triaxialHaloDensity(const Halo* h, mwvector pos)
{
    const real v  = mw_real_var(h->vhalo, 16);
    const real a  = mw_real_var(h->scaleLength, 17);
    const real q  = mw_real_var(h->flattenZ, 18);

    const real D   = mw_add(mw_add(mw_add(mw_add(sqr(a), mw_mul_s(sqr(X(pos)), (h->c1))), mw_mul_s(sqr(Y(pos)), (h->c2))), mw_mul_s(mw_mul(X(pos),Y(pos)), (h->c3))), sqr(mw_div(Z(pos),q)));

    const real num1 = mw_add(mw_mul_s(D, (2.0*(h->c1) + 2.0*(h->c2))), mw_mul_s(mw_div(D, sqr(q)), 2.0));
    const real num2 = sqr(mw_add(mw_mul_s(X(pos), 2.0*(h->c1)), mw_mul_s(Y(pos), (h->c3))));
    const real num3 = sqr(mw_add(mw_mul_s(Y(pos), 2.0*(h->c2)), mw_mul_s(X(pos), (h->c3))));
    const real num4 = sqr(mw_mul_s(mw_div(Z(pos), sqr(q)), 2.0));

    const real num = mw_sub(num1, mw_add(mw_add(num2, num3), num4));

    return mw_mul_s(mw_div(mw_mul(sqr(v),num),sqr(D)), inv_0(4.0*M_PI));
}

static inline real hernquistHaloDensity(const Halo* h,  real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);

    if(showRealValue(r) == 0 || showRealValue(a) == 0) return ZERO_REAL;

    return mw_mul_s(mw_div(mw_mul(M,a), mw_mul(r,cube(mw_add(r,a)))), inv_0(2*M_PI));
}

static inline real plummerHaloDensity(const Halo* h, real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);

    if(showRealValue(a) == 0)
    {
        return ZERO_REAL;
    }
    real r_a = mw_div(r,a);

    real rho_peak = mw_mul_s(mw_div(M,cube(a)), 3.0/(4.0*M_PI));

    return mw_mul(rho_peak, minusfivehalves(mw_add(mw_real_const(1.0), sqr(r_a))));
}

static inline real NFWMHaloDensity(const Halo* h,  real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);

    if(showRealValue(r) == 0) return ZERO_REAL;
    
    return mw_mul_s(mw_div(mw_div(M, r), sqr(mw_add(a,r))), inv_0(4.0*M_PI));

}

static inline real allenSantillanHaloDensity(const Halo* h, real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real_0 lam = h->lambda;
    const real gam = mw_real_const(h->gamma);

    real b = mw_sub(gam, mw_real_const(1.0));
    real r_a = mw_div(r,a);
    real numer = mw_mul(mw_pow(r_a,b), mw_add(mw_pow(r_a,b), gam));
    real denom = mw_mul(sqr(r),sqr(mw_add(mw_real_const(1.0), mw_pow(r_a,b))));

    if(showRealValue(r)==0.0) return ZERO_REAL;

    if(showRealValue(r) > lam) return ZERO_REAL;

    return mw_mul_s(mw_div(mw_div(mw_mul(M,numer),a),denom), inv_0(4.0*M_PI));
}

static inline real wilkinsonEvansHaloDensity(const Halo* h, real r)
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);

    real r_a = mw_div(r,a);
    real tmp = mw_hypot(r,a);

    if(showRealValue(r)==0.0) return ZERO_REAL;

    //(1/(4*pi)) * M/(sqr(r_a)*cube(tmp))
    return mw_mul_s(mw_div(M, mw_mul(sqr(r_a), cube(tmp))), inv_0(4*M_PI));
}

static inline real ninkovicHaloDensity(const Halo* h, real r)
{
    const real rho = mw_real_var(h->rho0, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    const real lam = mw_real_const(h->lambda);

    real r_a = mw_div(r,a);
    real lam_a  = mw_div(lam,a);

    const real density = mw_mul(rho,mw_sub(inv(mw_add(mw_real_const(1.0), cube(r_a))), inv(mw_add(mw_real_const(1.0), cube(lam_a)))));

    if(showRealValue(r) > showRealValue(lam))  return ZERO_REAL; 

    return density;
}

static inline real KVHalo(const Halo* h, real r) /*What is this one?*/
{
    const real M = mw_real_var(h->mass, 16);
    const real a = mw_real_var(h->scaleLength, 17);
    if(showRealValue(r) == 0) return ZERO_REAL;

    real ra = mw_add(r,a);

    return mw_mul_s(mw_mul(M, mw_sub(inv(mw_mul(r, sqr(ra))), mw_mul_s(inv(cube(ra)), 2.0))), inv_0(4*M_PI));
}

real nbExtDensity(const Potential* pot, mwvector pos, real_0 time)
{
    real density = ZERO_REAL;
    const real_0 limit = mw_pow_0(2.0,-8.0);

    /* Change r if less than limit. Done this way to pipeline this step*/
    real r = mw_add(mw_mul_s(mw_real_const(limit), (real_0)(showRealValue(mw_absv(pos)) <= limit)), mw_mul_s(mw_absv(pos), (real_0)(showRealValue(mw_absv(pos)) > limit)));

    switch (pot->sphere[0].type)
    {
        case HernquistSpherical:
            density = mw_add(density, hernquistSphericalDensity(&(pot->sphere[0]), r));
            break;
        case PlummerSpherical:
            density = mw_add(density, plummerSphericalDensity(&(pot->sphere[0]), r));
            break;
        case NoSpherical:
            break;
        case InvalidSpherical:
        default:
            mw_fail("Invalid bulge type in density\n");
    }

    switch (pot->disk.type)
    {
        case FreemanDisk: /*Density negligible since infinitely thin*/
            break;
        case MiyamotoNagaiDisk:
            density = mw_add(density, miyamotoNagaiDiskDensity(&(pot->disk), pos));
            break;
        case DoubleExponentialDisk:
            density = mw_add(density, doubleExponentialDiskDensity(&(pot->disk), pos));
            break;
        case Sech2ExponentialDisk:
            density = mw_add(density, sech2ExponentialDiskDensity(&(pot->disk), pos));
            break;
        case OrbitingBar:
            density = mw_add(density, orbitingBarDensity(&(pot->disk), pos, time));
        case NoDisk:
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid primary disk type in density\n");
    }

    switch (pot->disk2.type)
    {
        case FreemanDisk: /*Density negligible since infinitely thin*/
            break;
        case MiyamotoNagaiDisk:
            density = mw_add(density, miyamotoNagaiDiskDensity(&(pot->disk2), pos));
            break;
        case DoubleExponentialDisk:
            density = mw_add(density, doubleExponentialDiskDensity(&(pot->disk2), pos));
            break;
        case Sech2ExponentialDisk:
            density = mw_add(density, sech2ExponentialDiskDensity(&(pot->disk2), pos));
            break;
        case OrbitingBar:
            density = mw_add(density, orbitingBarDensity(&(pot->disk2), pos, time));
        case NoDisk:
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid primary disk type in density\n");
    }

    switch (pot->halo.type)
    {
        case LogarithmicHalo:
            density = mw_add(density, logarithmicHaloDensity(&(pot->halo), pos));
            break;
        case NFWHalo:
            density = mw_add(density, NFWHaloDensity(&(pot->halo), r));
            break;
        case TriaxialHalo:
            density = mw_add(density, triaxialHaloDensity(&(pot->halo), pos));
            break;
        case CausticHalo: /*FIXME: Add density profile for caustic halo when we actually plan on making this work*/
            break;
        case AllenSantillanHalo:
            density = mw_add(density, allenSantillanHaloDensity(&(pot->halo), r));
            break;
        case WilkinsonEvansHalo:
            density = mw_add(density, wilkinsonEvansHaloDensity(&(pot->halo), r));
	    break;
        case NFWMassHalo:
            density = mw_add(density, NFWMHaloDensity(&(pot->halo), r));
            break;
        case PlummerHalo:
            density = mw_add(density, plummerHaloDensity(&(pot->halo), r));
            break;
        case HernquistHalo:
            density = mw_add(density, hernquistHaloDensity(&(pot->halo), r));
            break;
        case NinkovicHalo:
            density = mw_add(density, ninkovicHaloDensity(&(pot->halo), r));
            break;
        case NoHalo:
            break;
        case InvalidHalo:
        default:
            mw_fail("Invalid halo type in density\n");
    }

    if (showRealValue(density) < 0.0)
    {
        mw_fail("Negative density calculated!\n    Faulty Potential = %s\n", showPotential(pot));
    }

    return density;

}

 
