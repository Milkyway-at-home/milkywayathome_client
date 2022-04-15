/*
 * Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
 * Copyright (c) 2010, 2011 Matthew Arsenault
 * Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute.
 * Copyright (c) 2002-2006 John M. Fregeau, Richard Campbell, Jeff Molofee
 * Copyright (c) 2016-2018 Siddhartha Shelton
 * 
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NBODY_H_
#define _NBODY_H_

#include "nbody_config.h"
#include "nbody_types.h"
#include "milkyway_util.h"

#include <time.h>


#ifdef _cplusplus
extern "C" {
#endif

/* Command line arguments */
typedef struct
{
    char* inputFile;
    char* outFileName;
    char* checkpointFileName;
    char* histogramFileName;
    char* histoutFileName;
    char* matchHistogram;   /* Just match this histogram to other histogram, no simulation */
    char* matchHistVelDisp;   /* Just match this histogram to other histogram, no simulation -- with vel dispersion calc*/
    char* matchHistBetaDisp;  /* Just match this histogram to other histogram, no simulation -- with beta dispersion calc*/
    char* matchHistBetaVelDisp; /* Just match this histogram to other histogram, no simulation -- with beta and vel dispersion calc*/
    char* matchHistBetaVlosDist; /* Just match this histogram to other histogram, no simulation -- with avg beta, avg vlos, dist calc*/
    char* matchHistAll;          /* Match this histogram to other histogram, no simulation -- with veta and vel disp, avg beta/vlos/dist */
    char* graphicsBin;
    char* visArgs;

    const char** forwardedArgs;
    unsigned int numForwardedArgs;

    int setSeed;  /* If the seed was specified or not */
    uint32_t seed;   /* Seed value */

    int numThreads;

    time_t checkpointPeriod;
    unsigned int platform;
    unsigned int devNum;

    /* These all must be int since that's the type popt expects them to be */
    int visualizer;
    int debugBOINC;
    int outputCartesian;
    int outputlbrCartesian;
    int printTiming;
    int verifyOnly;
    int printHistogram;  /* Print histogram at end */
    int outputBinary;
    int ignoreCheckpoint;

    int debugLuaLibs;   /* Open IO libraries etc. */
    int noCL;
    int reportProgress;
    int ignoreResponsive;
    int noCleanCheckpoint;
    int disableGPUCheckpointing;
    int verbose;
} NBodyFlags;

#define EMPTY_NBODY_FLAGS { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }

NBodyStatus nbStepSystem(const NBodyCtx* ctx, NBodyState* st);
NBodyStatus nbRunSystem(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf);
int nbVerifyFile(const NBodyFlags* nbf);
int nbMain(const NBodyFlags* nbf);
static NBodyCtx _ctx = EMPTY_NBODYCTX;
static NBodyState _st = EMPTY_NBODYSTATE;

#ifdef _cplusplus
}
#endif

#endif /* _NBODY_H_ */

