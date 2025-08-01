
#ifndef _NBODY_PLUMMER_MOMENTUM_H_
#define _NBODY_PLUMMER_MOMENTUM_H_

#include <lua.h>

//real enclosedMass (real M0, real k);

real velocityAdj (real a, 
                  mwvector pos,
                  mwvector vel,
                  real M0,
                  real k);

#endif /* _NBODY_PLUMMER_MOMENTUM_H_ */