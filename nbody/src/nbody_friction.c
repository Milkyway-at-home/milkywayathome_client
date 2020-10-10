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
    real a = (pot->halo)->scaleLength; 
    real v0 = 74.61;        /*FIXME: Calculate v0 for each halo type*/

    //object velocity where center of gravity is initially at rest
    real objectVel = mw_pow(mw_dotv(vel,vel), 0.5); 
    //for Coulomb logarithm
	real lambda = r / 1.6 / k; 
    //dispersion
    real sigma = isotropic_v_dispersion(a, v0, r); 
    //ratio of the velocity of the object to the modal velocity
    real X = objectVel / (mw_pow(2, 0.5) * sigma);

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

    //force from DF
    real F = (-4*pi*mw_pow(G_CONST, 2)*mw_pow(mass, 2)* mw_log(lambda)*density / mw_pow(objectVel, 2)) * (erf(X) - 2*X/mw_pow(pi, 0.5)*exp(-1.0*mw_pow(X, 2)));
    
    result.x = (F * vel.x / objectVel);
    result.y = (F * vel.y / objectVel);
    result.z = (F * vel.z / objectVel);
    return result;
}
