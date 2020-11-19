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
#include "nbody_potential.h"

#include "milkway_math.h"
#include "nbody_friction.h"

/** Isotropic Velocity Dispersion Formulas **/
real pi = 3.14159265;

static inline real dispIntegrand(const Potential* pot, const mwvector pos){
    mwvector acc = nbExtAcceleration(&pot, pos);

    real mag = mw_pow(mw_dotv(pos,pos), 0.5);
    mwvector unit_v = mw_mulvs(pos, 1.0/mag);

    real a_r = mw_abs(mw_dotv(acc,unit_v));

    real rho = nbExtDensity(&pot, pos);
    return rho*a_r;
}

static inline real velDispersion(const Potential* pot, const mwvector pos, real upperlimit){
    /** Use 5-point Gaussian Quadrature to calculate the RADIAL (1-dimensional) velocity dispersion integral **/
    int i,j;
    real a,b,integral,r;
    mwvector input_vec;
    const real weights[5];
    weights[0] = 0.2369268850561891;
    weights[1] = 0.4786286704993665;
    weights[2] = 0.5688888888888889;
    weights[3] = 0.4786286704993665;
    weights[4] = 0.2369268850561891;

    const real points[5];
    points[0] = -0.9061798459386640;
    points[1] = -0.5384693101056831;
    points[2] = 0.0;
    points[3] = 0.5384693101056831;
    points[4] = 0.9061798459386640;

    real dist = mw_pow(mw_dotv(pos,pos), 0.5);
    int nDivs = 10;
    real width = (upperlimit - dist)/(nDivs*1.0);
    integral = 0.0;
    for (i = 0; i < nDivs; i++) {
        a = dist + i*width;
        b = dist + (i+1)*width;
        for (j = 0; j < 5; j++) {
            r = width*points[j]/2.0 + (a+b)/2.0;
            input_vec = mw_mulvs(pos,r/dist);
            integral += weights[j]*dispIntegrand(&pot,input_vec);
        }
    }
    integral *= width/2.0;

    return integral/nbExtDensity(&pot,pos); /** Reutrns a velocity squared **/
}

/** Coulomb Logarithm for plummer spheres derived from Esquivel and Fuchs (2018) **/
static inline real CoulombLogPlummer(real scale_plummer, real scale_mwhalo){
    real u, c_log;
    if (scale_mwhalo==0) {
        return 0.0;
    }

    u = scale_plummer * (pi/scale_mwhalo); /** LMC scale radius times smallest wavenumber (k_min) **/
    c_log = u*u*( besselK0(u)*besselK2(u) - mw_pow(besselK1(u),2.0) )/2.0;

    return c_log;
}

/** Formula for Dynamical Friction from Merritt 2013 Eqn. 5.23 **/
mwvector dynamicalFriction_LMC(const Potential* pot, mwvector pos, mwvector vel, real mass_LMC, real scaleLength_LMC){
    mwvector result;        //Vector with acceleration due to DF

    real const G_CONST = 1; //(Time: Gyrs, Distance: kpc, Mass: SMU = 222288.47 solar masses)
    real thresh = mw_pow(2,-8);

    //object velocity where center of gravity is initially at rest
    real objectVel = mw_pow(mw_dotv(vel,vel), 0.5);
    objectVel = (objectVel >= thresh)*objectVel + (objectVel < thresh)*thresh; // To avoid divide-by-zero error.

    //Coloumb Logarithm
    real scaleLength_halo = (pot->halo)->scaleLength;
    real ln_lambda = CoulombLogPlummer(scaleLength_LMC, scaleLength_halo);

    //Calculate densities from each individual component
    real density = nbExtDensity(&pot, pos);

    //Get velocity dispersion of MW galaxy assuming isotropy
    real sigma2 = velDispersion(&pot, pos, 50.0*scaleLength_halo);
 
    //ratio of the velocity of the object to the modal velocity
    real X = objectVel / (mw_pow(2*sigma2, 0.5));

    //Acceleration from DF
    real acc = (-4 * pi * mw_pow(G_CONST, 2) * mass_LMC * ln_lambda * density / mw_pow(objectVel, 2)) * (erf(X) - 2*X/mw_pow(pi, 0.5)*exp(-1.0*mw_pow(X, 2)));
    
    result.x = (acc * vel.x / objectVel);
    result.y = (acc * vel.y / objectVel);
    result.z = (acc * vel.z / objectVel);
    return result;
}
