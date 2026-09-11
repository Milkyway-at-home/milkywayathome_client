#ifndef _NBODY_MATH_FUNCS_H_
#define _NBODY_MATH_FUNCS_H_

#include "milkyway_math.h"
#include "nbody_types.h"
#include "nbody_potential_types.h"

typedef real (*ODE2ndDeriv)(real, real, real, const Dwarf *params);
real ODE2ndOrderSolver(real xEval, int stepsPerx, real yInit, real yPrimeInit, ODE2ndDeriv f, const Dwarf* params, int returnXWhen0);

real GammaFunc(const real z);

real UpperIncompleteGammaFunc(real a, real x);
real LowerIncompleteGammaFunc(real a, real x);
real ErrorFunc(real x);
real ComplementaryErrorFunc(real x);
real ComplementaryErrorFuncApprox(real x);

real first_derivative(real (*func)(const Dwarf*, real), real x, const Dwarf* comp1);
real second_derivative(real (*func)(const Dwarf*, real), real x, const Dwarf* comp1);
real gauss_quad(real (*func)(real, const Dwarf*, const Dwarf*, real, mwbool), real lower, real upper, const Dwarf* comp1, const Dwarf* comp2, real energy, mwbool isDark);

#endif /* _NBODY_MATH_FUNCS_H_ */