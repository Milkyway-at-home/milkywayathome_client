#ifndef _NBODY_BAR_TIME_H_
#define _NBODY_BAR_TIME_H_

#include "nbody_types.h"

//using XYZ coordinates and rotation angles, fills st->backwardOrbitAngles
void fillBackwardOrbitAngles(NBodyState* st);

//calculates the orbit rotation needed to align our coordinates
//with the plane of the star stream
void setBackwardOrbitRotation(NBodyState* st);

real meanBodyAngle(Body *bodies, int size, NBodyState* st);

//also updates the state with the new bar timestep
int getBarTime(Body* bodies, int nbody, NBodyState* st, NBodyCtx* ctx);

real getAngleDiff(real a1, real a2);

mwbool angleIsBetween(real start, real end, real mid);


mwvector getStreamCenter(NBodyState* st, mwvector* meanBinCenter, 
mwvector* histCenterVelocity, mwvector* meanBinVelocity);
#endif