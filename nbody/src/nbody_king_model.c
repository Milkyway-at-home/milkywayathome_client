#include "milkyway_util.h"
#include "milkyway_math.h"
#include "milkyway_lua.h"
#include "nbody_potential_types.h"
#include "nbody_king_model.h"
/*NOTE: These functions use mw_erf which needs further testing to determine whether there are differences between different OS.*/

// Generic function that numerically solves 2nd order ODEs of the form y'' = f(x, y(x), y'(x)), where y' is dy/dx
// Uses 4th order Runge-Kutta numerical method, function input is of the form of ODE2ndDeriv
// The last parameter returnXWhen0 is a special case boolean where the function will instead return the x value where the otherwise positive y(x) function hits zero/negative
real ODE2ndOrderSolver(real xEval, int stepsPerx, real yInit, real yPrimeInit, ODE2ndDeriv f, Dwarf* params, int returnXWhen0) {
    real nSteps = floor(stepsPerx * xEval);
    real stepRes = stepsPerx*xEval - nSteps;
    real deltax = xEval/nSteps;
    
    real y = yInit; // current value for y(x)
    real z = yPrimeInit; // let z be y'(x)
    real* yCurr = &y;
    real* zCurr = &z;
    real k1, k2, k3, k4, l1, l2, l3, l4;
    real x = deltax;
    // the the last loop needs to be one before the desired number of steps, since the end of a loop gives the values for the next step
    for (int n = 1; n < nSteps + 1; n++) {
        // l factors adjust z, k factors adjust y
        x = n * deltax;

        l1 = deltax * f(x, y, z, params);
        k1 = deltax * z;

        l2 = deltax * f(x + 0.5*deltax, y + 0.5*k1, z + 0.5*l1, params);
        k2 = deltax * (z + 0.5*l1);

        l3 = deltax * f(x + 0.5*deltax, y + 0.5*k2, z + 0.5*l2, params);
        k3 = deltax * (z + 0.5*l2);

        l4 = deltax * f(x + deltax, y + k3, z + l3, params);
        k4 = deltax * (z + l3);

        *yCurr = *yCurr + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        *zCurr = *zCurr + (1.0/6.0)*(l1 + 2.0*l2 + 2.0*l3 + l4);

        if ((*yCurr <= 0.0 || isnan(*yCurr) || !isfinite(*yCurr)) && returnXWhen0 == 1) {
            stepRes = 0.0;
            *yCurr = x;
            break;
        }

    }

    // run the iteration one more time if there is a fractional step at the end
    if (stepRes > 0.0) {
        deltax = stepRes/stepsPerx; // smaller than regular deltax
        l1 = deltax * f(x, y, z, params);
        k1 = deltax * z;
        l2 = deltax * f(x + 0.5*deltax, y + 0.5*k1, z + 0.5*l1, params);
        k2 = deltax * (z + 0.5*l1);
        l3 = deltax * f(x + 0.5*deltax, y + 0.5*k2, z + 0.5*l2, params);
        k3 = deltax * (z + 0.5*l2);
        l4 = deltax * f(x + deltax, y + k3, z + l3, params);
        k4 = deltax * (z + l3);
        *yCurr = *yCurr + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        *zCurr = *zCurr + (1.0/6.0)*(l1 + 2.0*l2 + 2.0*l3 + l4);
    }
    return *yCurr;
}


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