

#ifndef _NBODY_KING_MODEL_H_
#define _NBODY_KING_MODEL_H_

#include "nbody_types.h"
#include "milkyway_math.h"
#include "nbody_potential_types.h"

typedef real (*ODE2ndDeriv)(real, real, real, Dwarf *params);

real ODE2ndOrderSolver(real xEval, int stepsPerx, real yInit, real yPrimeInit, ODE2ndDeriv f, Dwarf* params, int returnXWhen0);
real kingDimlessRho(real W, real W0);
real kingDimless2ndDeriv(real R, real W, real dWdR, Dwarf *model);
real kingDimlessMass(real R, Dwarf* model, Dwarf* unusedModel, real unusedEnergy, mwbool unusedIsDark);
real kingDensityFromPsi(real psi, real sig, real rho1);
real kingRelPot2ndDeriv(real r, real psi, real dPsidr, Dwarf *model);


#endif