#ifndef _BLENDER_VISUALIZER_H_
#define _BLENDER_VISUALIZER_H_
#include <math.h>
#include "nbody_util.h"
#include "nbody_types.h"


int nbFindCenterOfMass(mwvector* cmPos, const NBodyState* st);
NBodyStatus deleteOldFiles(const NBodyState* st);
NBodyStatus blenderPrintBodies(const NBodyState* st, const NBodyCtx* ctx);
NBodyStatus blenderPrintCOM(const NBodyState* st);
NBodyStatus blenderPrintMisc(const NBodyState* st,const NBodyCtx* ctx, mwvector cmPos1, mwvector cmPos2);
NBodyStatus blenderPossiblyChangePerpendicularCmPos(mwvector* next, mwvector* perp, mwvector* start);
NBodyStatus blenderWriteNewEvents(const NBodyState* st);

#endif