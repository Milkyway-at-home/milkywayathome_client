/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
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
    char* visArgs;

    const char** forwardedArgs;
    unsigned int numForwardedArgs;

    mwbool setSeed;  /* If the seed was specified or not */
    uint32_t seed;   /* Seed value */

    int numThreads;

    time_t checkpointPeriod;
    unsigned int platform;
    unsigned int devNum;

    mwbool visualizer;
    mwbool debugBOINC;
    mwbool outputCartesian;
    mwbool printTiming;
    mwbool verifyOnly;
    mwbool printBodies;
    mwbool printHistogram;
    mwbool cleanCheckpoint;
    mwbool ignoreCheckpoint;

    mwbool debugLuaLibs;   /* Open IO libraries etc. */
    mwbool noCL;
    mwbool reportProgress;
    mwbool ignoreResponsive;
} NBodyFlags;

#define EMPTY_NBODY_FLAGS { NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }

int verifyFile(const NBodyFlags* nbf);
int runNBodySimulation(const NBodyFlags* nbf);

#ifdef _cplusplus
}
#endif

#endif /* _NBODY_H_ */

