

#ifndef _NBODY_KING_MODEL_H_
#define _NBODY_KING_MODEL_H_

#include "nbody_types.h"
#include "milkyway_math.h"
#include "nbody_potential_types.h"

typedef real (*ODE2ndDeriv)(real, real, real, Dwarf *params);

typedef struct {
    real x;
    real y;
    real dydx;
} ODE2ndOrderVals;

ODE2ndOrderVals ODE2ndOrderSolver(real xEval, int stepsPerx, real xInit, real yInit, real yPrimeInit, ODE2ndDeriv f, Dwarf* params, int stopWhenZero);
real interpolateLinear(real x, real x0, real x1, real f0, real f1);
real kingDimlessRho(real W, real W0);
real kingDimless2ndDeriv(real R, real W, real dWdR, Dwarf *model);
real kingDimlessMass(real R, Dwarf* model, Dwarf* unusedModel, real unusedEnergy, mwbool unusedIsDark);



#endif