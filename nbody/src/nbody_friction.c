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

#include "milkway_math.h"
#include "nbody_friction.h"

/** Helper function
	Calculates the integral of Maxwell's distrubution DF formula**/

/*FIXME: Find v_dispersion formulas for each type of halo*/
real isotropic_v_dispersion(real a, real v0, real r){
    return (v0*v0) / 2 * (4 * r * (a+r) * (a+r) * mw_log(1 + r/a) - a * r * (5*a + 4*r)) / (a*a * (2*a + r));
}

/** Adds the change in velocity that is the result of DF from a massive body. 
	Uses the Maxwell's distribution dynamical friction formula.**/
mwvector dynamicalFriction(mwvector pos, mwvector vel, real mass, const Potential* pot){
    mwvector result;        //Vector with change in velocity due to DF
    real density = 0.0;     //Galactic density placeholder

    real const G_CONST = 1; //(Time: Gyrs, Distance: kpc, Mass: SMU = 222288.47 solar masses)
    real k = 1.428;         /*FIXME: Need better way of calculating this quantity*/
    real r = mw_pow(mw_dotv(pos,pos), 0.5); 
    real a;
    real v0 = 0.0;
    real sigma = 0.0;
    real M = 0.0;
    real thresh = mw_pow(2,-8);

    //object velocity where center of gravity is initially at rest
    real objectVel = mw_pow(mw_dotv(vel,vel), 0.5); 
    //for Coulomb logarithm
    real lambda = r / 1.6 / k; 

    //Calculate densities from each individual component

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
            a = (pot->halo)->scaleLength;
            v0 = (pot->halo)->vhalo;
            sigma = mw_sqrt((v0*v0) / 2 * (4 * r * (a+r) * (a+r) * mw_log(1 + r/a) - a * r * (5*a + 4*r)) / (a*a * (2*a + r)));
            break;
        case NFWHalo:
            density += NFWHaloDensity(&(pot->halo), r);
            a = (pot->halo)->scaleLength;
            v0 = (pot->halo)->vhalo;
            sigma = 0.0; /*FIXME: Find velocity dispersion of NFW Profile*/
            break;
        case TriaxialHalo:
            density += triaxialHaloDensity(&(pot->halo), pos);
            a = (pot->halo)->scaleLength;
            v0 = (pot->halo)->vhalo;
            sigma = 0.0; /*FIXME: Find velocity dispersion of Triaxial Profile*/
            break;
        case CausticHalo:
            density += 0.0; /*FIXME: Add density profile for caustic halo when we actually plan on making this work*/
            a = (pot->halo)->scaleLength;
            v0 = 0.0;
            break;
        case AllenSantillanHalo:
            density += allenSantillanHaloDensity(&(pot->halo), r);
            a = (pot->halo)->scaleLength;
            M = 0.0;
            sigma = 0.0; /*FIXME: Find velocity dispersion of Allen-Santillan Profile*/
            break;
        case WilkinsonEvansHalo:
            density += wilkinsonEvansHaloDensity(&(pot->halo), r);
            a = (pot->halo)->scaleLength;
            M = 0.0;
            sigma = 0.0; /*FIXME: Find velocity dispersion of Wilkinson-Evans Profile*/
	    break;
        case NFWMassHalo:
            density += NFWMHaloDensity(&(pot->halo), r);
            a = (pot->halo)->scaleLength;
            M = 0.0;
            sigma = 0.0; /*FIXME: Find velocity dispersion of NFW Profile*/
            break;
        case PlummerHalo:
            density += plummerHaloDensity(&(pot->halo), r);
            a = (pot->halo)->scaleLength;
            M = (pot->halo)->mass;
            sigma = mw_sqrt(mass/mw_sqrt(r*r + a*a)/6.0);
            break;
        case HernquistHalo:
            density += hernquistHaloDensity(&(pot->halo), r);
            a = (pot->halo)->scaleLength;
            M = (pot->halo)->mass;
            real R = (!(r<thresh))*r/a + (r<thresh)*thresh;
            sigma = mw_sqrt((M/a)*(R*mw_pow(1+R,3)*mw_log((1+R)/R) - R*(25.0 + 52.0*R + 42.0*R*R + 12.0*R*R*R)/12.0/(1+R)));
            break;
        case NinkovicHalo:
            density += ninkovicHaloDensity(&(pot->halo), r);
            a = (pot->halo)->scaleLength;
            v0 = 0.0;
            sigma = 0.0; /*FIXME: Find velocity dispersion of Ninkovic Profile*/
            break;
        case NoHalo:
            density += 0.0;
            a = (pot->halo)->scaleLength;
            v0 = 0.0;
            break;
        case InvalidHalo:
        default:
            mw_fail("Invalid halo type in density\n");
    }
 
    //ratio of the velocity of the object to the modal velocity
    real X = objectVel / (mw_pow(2, 0.5) * sigma);

    //acceleration from DF
    real acc = (-4*pi*mw_pow(G_CONST, 2)*mass* mw_log(lambda)*density / mw_pow(objectVel, 2)) * (erf(X) - 2*X/mw_pow(pi, 0.5)*exp(-1.0*mw_pow(X, 2)));
    
    result.x = (acc * vel.x / objectVel);
    result.y = (acc * vel.y / objectVel);
    result.z = (acc * vel.z / objectVel);
    return result;
}
