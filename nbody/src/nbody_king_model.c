#include "milkyway_util.h"
#include "milkyway_math.h"
#include "milkyway_lua.h"
#include "nbody_math_funcs.h"
#include "nbody_potential_types.h"
#include "nbody_king_model.h"
/*NOTE: These functions use mw_erf which needs further testing to determine whether there are differences between different OS.*/
/*UPDATE: mw_erf uses the CORE-MATH function for consistent rounding from version 1.97 on. This should solve any issues*/


real kingDimlessRho(real W, real W0) {
    real rho = mw_exp(W)*mw_erf(mw_sqrt(W)) - mw_sqrt(4.0*W/M_PI)*(1.0 + (2.0/3.0)*W);
    real rho0 = mw_exp(W0)*mw_erf(mw_sqrt(W0)) - mw_sqrt(4.0*W0/M_PI)*(1.0 + (2.0/3.0)*W0);
    return rho/rho0;
}

// King model dimensionless Poisson equation, evaluated for 2nd derivative term
real kingDimless2ndDeriv(real R, real W, real dWdR, Dwarf *model) {
    real W0 = model->W0; 
    real dimlessRho = kingDimlessRho(W, W0);

    return (-9.0) * dimlessRho + (-2.0/R) * (dWdR);  
}

// This function is formatted in such a way that gauss_quad() will accept it, the integrand for mu parameter
// parameters are radius, Dwarf (used), Dwarf (unused), energy (unused), isDark (unused)
real kingDimlessMass(real R, Dwarf* model, Dwarf* unusedModel, real unusedEnergy, mwbool unusedIsDark) {
    real W_R = ODE2ndOrderSolver(R, 1000, model->W0, 0.0, kingDimless2ndDeriv, model, 0);
    return kingDimlessRho(W_R, model->W0) * 4.0 * M_PI * R * R;
}

// Equation 4.111 from Binney & Tremaine 2nd ed.
real kingDensityFromPsi(real psi, real sig, real rho1) {
    real erfTerm = mw_exp(psi/(sig*sig)) * mw_erf(mw_sqrt(psi/(sig*sig)));
    return rho1 * (erfTerm - mw_sqrt(4.0*psi/(M_PI*sig*sig))*(1.0 + (2.0*psi)/(3.0*sig*sig)));
}

// Equation 4.112 from Binney & Tremaine 2nd ed.
real kingRelPot2ndDeriv(real r, real psi, real dPsidr, Dwarf *model) {
    real rho1 = model->rho1;
    real sigma = model->sigma;
    real rhs = -4.0*M_PI*kingDensityFromPsi(psi, sigma, rho1);
    return rhs + (-2.0/r)*dPsidr;
}