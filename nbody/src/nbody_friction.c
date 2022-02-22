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

static inline real dispIntegrand(const Potential* pot, const mwvector* pos, real_0 time){
    //mw_printf("Getting External Acceleration\n");
    mwvector unit_v;
    mwvector acc = nbExtAcceleration(pot, pos, time);

    real tmp = mw_dotv(pos,pos);
    real mag = mw_sqrt(&tmp);
    tmp = inv(&mag);
    unit_v.x = mw_mul(&X(pos),&tmp);
    unit_v.y = mw_mul(&Y(pos),&tmp);
    unit_v.z = mw_mul(&Z(pos),&tmp);

    tmp = mw_dotv(&acc, &unit_v);
    real a_r = mw_abs(&tmp);

    //mw_printf("Getting External Density\n");
    real rho = nbExtDensity(pot, pos, time);
    return mw_mul(&rho, &a_r);
}

static inline real velDispersion(const Potential* pot, const mwvector* pos, real_0 upperlimit, real_0 time){
    /** Use 5-point Gaussian Quadrature to calculate the RADIAL (1-dimensional) velocity dispersion integral **/
    int i,j;
    real a,b,integral,r;
    real tmp1, tmp2;
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

    real dist = mw_length(pos);
    //mw_printf("POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(pos)), showRealValue(&Y(pos)), showRealValue(&Z(pos)) );
    //mw_printf("dist = %.15f\n", showRealValue(&dist));
    int nDivs = 10;

    tmp1 = mw_add_s(&dist, -upperlimit);
    tmp1 = mw_neg(&tmp1);
    real width = mw_mul_s(&tmp1, inv_0(nDivs*1.0));
    real rho0 = nbExtDensity(pot, pos, time);
    //mw_printf("rho0 = %.15f\n", showRealValue(&rho0));
    if ((showRealValue(&width) <= 0)||(showRealValue(&rho0) == 0.0)) /** To avoid divide by zero later **/
    {
        return mw_real_const(-1);       /** Want to ignore any contribution more than 50 scale radii from galactic center. Chose this value because dispersion is never negative. **/
    }
    integral = ZERO_REAL;
    for (i = 0; i < nDivs; i++) {
        tmp1 = mw_mul_s(&width,(real_0) i);
        a = mw_add(&dist, &tmp1);

        tmp1 = mw_mul_s(&width,(real_0) (i+1));
        b = mw_add(&dist, &tmp1);

        for (j = 0; j < 5; j++) {
            tmp1 = mw_mul_s(&width, points[j]/2.0);
            tmp2 = mw_add(&a,&b);
            tmp2 = mw_mul_s(&tmp2, 0.5);
            r = mw_add(&tmp1, &tmp2);

            tmp1 = mw_div(&r, &dist);
            input_vec.x = mw_mul(&X(pos), &tmp1);
            input_vec.y = mw_mul(&Y(pos), &tmp1);
            input_vec.z = mw_mul(&Z(pos), &tmp1);

            //mw_printf("input_vec = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&input_vec)), showRealValue(&Y(&input_vec)), showRealValue(&Z(&input_vec)) );

            tmp1 = dispIntegrand(pot, &input_vec, time);
            tmp1 = mw_mul(&tmp1, &width);
            tmp1 = mw_mul_s(&tmp1, 0.5 * weights[j]);
            integral = mw_add(&integral, &tmp1);
        }
    }

    //mw_printf("integral = %.15f\n", showRealValue(&integral));
    real sig2 = mw_div(&integral, &rho0);
 
    return sig2; /** Returns a velocity squared **/
}

/** Coulomb Logarithm for plummer spheres derived from Esquivel and Fuchs (2018) **/
static inline real CoulombLogPlummer(real* scale_plummer, real* scale_mwhalo){
    real u, c_log;
    real tmp1, tmp2, tmp3;

    tmp1 = inv(scale_mwhalo);
    tmp1 = mw_mul_s(&tmp1, M_PI);
    u = mw_mul(scale_plummer, &tmp1); /** LMC scale radius times smallest wavenumber (k_min) **/

    tmp1 = sqr(&u);
    tmp2 = mw_besselK0(&u);
    tmp3 = mw_besselK2(&u);
    tmp2 = mw_mul(&tmp2, &tmp3);
    tmp3 = mw_besselK1(&u);
    tmp3 = sqr(&tmp3);
    tmp2 = mw_sub(&tmp2, &tmp3);
    tmp1 = mw_mul(&tmp1, &tmp2);
    c_log = mw_mul_s(&tmp1, 0.5);

    return c_log; /** Normally, the Coulomb Log is between 3 and 30. But this value seems closer to 0.0003. **/
}

static inline real getHaloScaleLength(const Halo* halo){
    real scale = mw_real_var(halo->scaleLength, HALO_RADIUS_POS);
    if (showRealValue(&scale) == 0.0)
    {
        scale = mw_real_var(1.0,HALO_RADIUS_POS);    /** This is the case when there is no halo **/
    }
    return scale;
}

/** Formula for Dynamical Friction using Chandrasekhar's formula and assuming an isotropic Maxwellian velocity distribution **/
mwvector dynamicalFriction_LMC(const Potential* pot, mwvector* pos, mwvector* vel, real* mass_LMC, real* scaleLength_LMC, mwbool dynaFric, real_0 time)
{
    mwvector result;        //Vector with acceleration due to DF
    if (!dynaFric) {
        result.x = ZERO_REAL;
        result.y = ZERO_REAL;
        result.z = ZERO_REAL;
        return result;
    }

    real X, tmp1, tmp2;

    const real_0 G_CONST = 1; //(Time: Gyrs, Distance: kpc, Mass: SMU = 222288.47 solar masses)
    const real_0 thresh_val = mw_pow_0(2,-8);

    //object velocity where center of gravity is initially at rest
    real objectVel = mw_absv(vel);
    real_0 small_amount = (showRealValue(&objectVel) <= thresh_val)*thresh_val;
    objectVel = mw_add_s(&objectVel, small_amount); // To avoid divide-by-zero error.
    //mw_printf("objectVel = %.15f\n", showRealValue(&objectVel));

    //Coloumb Logarithm
    real scaleLength_halo = getHaloScaleLength(&(pot->halo));
    //mw_printf("a = %.15f\n", showRealValue(&scaleLength_halo));
    real ln_lambda = CoulombLogPlummer(scaleLength_LMC, &scaleLength_halo);
    //mw_printf("ln_lambda = %.15f\n", showRealValue(&ln_lambda));

    //Calculate densities from each individual component
    real density = nbExtDensity(pot, pos, time);
    //mw_printf("density   = %.15f\n", showRealValue(&density));

    //Get velocity dispersion of MW galaxy assuming isotropy
    real sigma2 = velDispersion(pot, pos, 50.0*showRealValue(&scaleLength_halo), time);
    //mw_printf("sig2 = %.15f\n", showRealValue(&sigma2));
 
    //ratio of the velocity of the object to the MW modal velocity
    if (showRealValue(&sigma2) <= 0.0) {
        X = ZERO_REAL;
    }
    else {
        tmp1 = mw_mul_s(&sigma2, 2);
        tmp1 = mw_sqrt(&tmp1);
        X = mw_div(&objectVel, &tmp1);
    }

    tmp1 = sqr(&objectVel);
    tmp1 = mw_div(&density, &tmp1);
    tmp1 = mw_mul(&ln_lambda, &tmp1);
    tmp1 = mw_mul(mass_LMC, &tmp1);
    real factor = mw_mul_s(&tmp1, -4.0 * M_PI);
    //mw_printf("mass_LMC  = %.15f\n", showRealValue(mass_LMC));

    tmp1 = mw_erf(&X);
    tmp2 = sqr(&X);
    tmp2 = mw_neg(&tmp2);
    tmp2 = mw_exp(&tmp2);
    tmp2 = mw_mul(&X, &tmp2);
    tmp2 = mw_mul_s(&tmp2, 2.0/mw_sqrt_0(M_PI));
    real erfStuff = mw_sub(&tmp1, &tmp2);

    //mw_printf("erfStuff = %.15f\n", showRealValue(&erfStuff));
    //mw_printf("factor   = %.15f\n", showRealValue(&factor));

    //Acceleration from DF
    real acc = mw_mul(&factor, &erfStuff);

    tmp1 = mw_div(&acc, &objectVel);
    result.x = mw_mul(&vel->x, &tmp1);
    result.y = mw_mul(&vel->y, &tmp1);
    result.z = mw_mul(&vel->z, &tmp1);

    //mw_printf("ACC_DF = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&result)), showRealValue(&Y(&result)), showRealValue(&Z(&result)) );

    return result;
}
