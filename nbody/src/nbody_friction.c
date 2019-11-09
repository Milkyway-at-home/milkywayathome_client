//File for Dynamical Friction and potentials
#include "nbody_priv.h"
#include "nbody_potential.h"
#include "nbody_density.h"
#include "nbody_potential_types.h"
#include "milkyway_util.h"
#include "nbody_caustic.h"
#include "nbody_bessel.h"

#include "milkway_math.h"
#include "friction.h"

/** Helper function
	Calculates the integral of Maxwell's distrubution DF formula**/
real isotropic_v_dispersion(real a, real v0, real r){
    return (v0*v0) / 2 * (4 * r * (a+r) * (a+r) * mw_log(1 + r/a) - a * r * (5*a + 4*r)) / (a*a * (2*a + r));
}

/** Adds the change in velocity that is the result of DF from a massive body. 
	Uses the Maxwell's distribution dynamical friction formula.**/
mwvector dynamicalFriction(mwvector pos, mwvector vel, real mass, const Potential* pot){
    mwvector result; //Vector with change in velocity due to DF

    real const G_CONST = 6.67*mw_pow(10, -11); //figure out better place to get this from
    real k = 1.428; 
    real r = mw_pow(mw_dotv(pos,pos), 0.5); 
    int a = 12; 
    int v0 = 73; 

    //object velocity where center of gravity is initially at rest
    real objectVel = mw_pow(mw_dotv(vel,vel), 0.5); 
    //for Coulomb logarithm
	real lambda = r / 1.6 / k; 
    //dispersion
    real sigma = isotropic_v_dispersion(a, v0, r); 
    //ratio of the velocity of the object to the modal velocity
    real X = objectVel / (mw_pow(2, 0.5) * sigma);
    //density of the field
    real density = NFWMHaloDensity(&(pot->halo), r);//milky halo density LMC_mass * temp_N;
    //force from DF
    real F = (-4*pi*mw_pow(G_CONST, 2)*mw_pow(mass, 2)* mw_log(lambda)*density 
    	/ mw_pow(objectVel, 3)) * (erf(X) - 2*X/mw_pow(pi, 0.5)*exp(mw_pow(-1*X, 2)));
    
    result.x = (F * vel.x / pow(mw_dotv(vel,vel), 0.5));
    result.y = (F * vel.y / pow(mw_dotv(vel,vel), 0.5));
    result.z = (F * vel.z / pow(mw_dotv(vel,vel), 0.5));
    return result;
}
