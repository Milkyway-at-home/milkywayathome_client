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

#include "nbody_friction.h"

/** Isotropic Velocity Dispersion Formulas **/

static inline real dispIntegrand(const Potential* pot, const mwvector pos, real_0 time){
    //mw_printf("Getting External Acceleration\n");
    mwvector acc = nbExtAcceleration(pot, pos, time);

    real mag = mw_sqrt(mw_dotv(pos,pos));
    mwvector unit_v = mw_mulvs(pos, inv(mag));

    real a_r = mw_abs(mw_dotv(acc,unit_v));

    //mw_printf("Getting External Density\n");
    real rho = nbExtDensity(pot, pos, time);
    return mw_mul(rho, a_r);
}

static inline real velDispersion(const Potential* pot, const mwvector pos, real_0 upperlimit, real_0 time){
    /** Use 5-point Gaussian Quadrature to calculate the RADIAL (1-dimensional) velocity dispersion integral **/
    int i,j;
    real a,b,integral,r;
    mwvector input_vec;
    real_0 weights[5];
    weights[0] = 0.2369268850561891;
    weights[1] = 0.4786286704993665;
    weights[2] = 0.5688888888888889;
    weights[3] = 0.4786286704993665;
    weights[4] = 0.2369268850561891;

    real_0 points[5];
    points[0] = -0.9061798459386640;
    points[1] = -0.5384693101056831;
    points[2] = 0.0;
    points[3] = 0.5384693101056831;
    points[4] = 0.9061798459386640;

    real dist = mw_sqrt(mw_dotv(pos,pos));
    int nDivs = 10;
    real width = mw_mul_s(mw_sub(mw_real_const(upperlimit), dist), inv_0(nDivs*1.0));
    real rho0 = nbExtDensity(pot, pos, time);
    if ((showRealValue(width) <= 0)||(showRealValue(rho0) == 0.0)) /** To avoid divide by zero later **/
    {
        return mw_real_const(-1);       /** Want to ignore any contribution more than 50 scale radii from galactic center. Chose this value because dispersion is never negative. **/
    }
    integral = ZERO_REAL;
    for (i = 0; i < nDivs; i++) {
        a = mw_add(dist, mw_mul_s(width,i));
        b = mw_add(dist, mw_mul_s(width,i+1));
        for (j = 0; j < 5; j++) {
            r = mw_add(mw_mul_s(width, points[j]/2.0), mw_mul_s(mw_add(a,b), 0.5));
            input_vec = mw_mulvs(pos, mw_div(r, dist));
            integral = mw_add(integral, mw_mul_s(mw_mul(dispIntegrand(pot, input_vec, time), width), 0.5 * weights[j]));
        }
    }

    return mw_div(integral, rho0); /** Reutrns a velocity squared **/
}

/** Coulomb Logarithm for plummer spheres derived from Esquivel and Fuchs (2018) **/
static inline real CoulombLogPlummer(real scale_plummer, real scale_mwhalo){
    real u, c_log;

    u = mw_mul(scale_plummer, mw_mul_s(inv(scale_mwhalo), M_PI)); /** LMC scale radius times smallest wavenumber (k_min) **/
    c_log = mw_mul_s(mw_mul(sqr(u), mw_sub(mw_mul(besselK0(u), besselK2(u)), sqr(besselK1(u)))), 0.5);

    return c_log; /** Normally, the Coulomb Log is between 3 and 30. But this value seems closer to 0.0003. **/
}

static inline real getHaloScaleLength(const Halo* halo){
    real scale = mw_real_var(halo->scaleLength, HALO_RADIUS_POS);
    if (showRealValue(scale) == 0.0)
    {
        scale = mw_real_var(1.0,HALO_RADIUS_POS);    /** This is the case when there is no halo **/
    }
    return scale;
}

/** Formula for Dynamical Friction using Chandrasekhar's formula and assuming an isotropic Maxwellian velocity distribution **/
mwvector dynamicalFriction_LMC(const Potential* pot, mwvector pos, mwvector vel, real mass_LMC, real scaleLength_LMC, mwbool dynaFric, real_0 time)
{
    mwvector result = mw_vec(ZERO_REAL, ZERO_REAL, ZERO_REAL);        //Vector with acceleration due to DF
    if (!dynaFric) {
        return result;
    }

    real X;
    Halo *mw_halo;

    const real_0 G_CONST = 1; //(Time: Gyrs, Distance: kpc, Mass: SMU = 222288.47 solar masses)
    const real_0 thresh_val = mw_pow_0(2,-8);

    //object velocity where center of gravity is initially at rest
    real objectVel = mw_absv(vel);
    objectVel = mw_add(objectVel, mw_real_const((showRealValue(objectVel) <= thresh_val)*thresh_val)); // To avoid divide-by-zero error.

    //Coloumb Logarithm
    mw_halo = &(pot->halo);
    real scaleLength_halo = getHaloScaleLength(mw_halo);
    //mw_printf("a = %.15f\n", scaleLength_halo);
    real ln_lambda = CoulombLogPlummer(scaleLength_LMC, scaleLength_halo);
    //mw_printf("ln(L) = %.15f\n", ln_lambda);

    //Calculate densities from each individual component
    real density = nbExtDensity(pot, pos, time);
    //mw_printf("rho = %.15f\n", density);

    //Get velocity dispersion of MW galaxy assuming isotropy
    real sigma2 = velDispersion(pot, pos, 50.0*showRealValue(scaleLength_halo), time);
    //mw_printf("sig2 = %.15f\n", sigma2);
 
    //ratio of the velocity of the object to the MW modal velocity
    if (showRealValue(sigma2) <= 0.0) {
        X = ZERO_REAL;
    }
    else {
        X = mw_div(objectVel, mw_sqrt(mw_mul_s(sigma2, 2)));
    }

    real factor = mw_mul_s(mw_mul(mass_LMC, mw_mul(ln_lambda, mw_div(density, sqr(objectVel)))), -4.0 * M_PI);
    real erfStuff = mw_sub(mw_erf(X), mw_mul_s(mw_mul(X, mw_exp(mw_neg(sqr(X)))), 2.0/mw_sqrt_0(M_PI)));

    //Acceleration from DF
    real acc = mw_mul(factor, erfStuff);
    
    result = mw_mulvs(vel, mw_div(acc, objectVel));
    return result;
}
