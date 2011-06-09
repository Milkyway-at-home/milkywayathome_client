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

#ifndef _SEPARATION_H_
#define _SEPARATION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "separation_config.h"
#include "separation_types.h"
#include "milkyway_util.h"
#include "milkyway_math.h"
#include "evaluation.h"
#include "evaluation_state.h"
#include "parameters.h"
#include "star_points.h"
#include "calculated_constants.h"
#include "separation_constants.h"
#include "r_points.h"
#include "gauss_legendre.h"
#include "io_util.h"
#include "coordinates.h"
#include "integrals.h"
#include "likelihood.h"
#include "separation_utils.h"

#if SEPARATION_OPENCL
  #include "run_cl.h"
#endif /* SEPARATION_OPENCL */

#include <limits.h>

typedef struct
{
    char* star_points_file;
    char* ap_file;  /* astronomy parameters */
    char* separation_outfile;
    char* preferredPlatformVendor;
    const char** forwardedArgs;
    real* numArgs;   /* Temporary */
    unsigned int nForwardedArgs;
    int debugBOINC;
    int do_separation;
    int setSeed;
    int separationSeed;
    int cleanupCheckpoint;
    int ignoreCheckpoint;  /* Ignoring checkpoint is not the same as disabling GPU checkpoints */
    unsigned int usePlatform;
    unsigned int useDevNumber;  /* Choose CL platform and device */
    int nonResponsive;
    double targetFrequency;
    int pollingMode;
    int disableGPUCheckpointing;

    MWPriority processPriority;
    int setPriority;

    /* Force between normal, SSE2, SSE3 paths */
    int forceNoIntrinsics;
    int forceX87;
    int forceSSE2;
    int forceSSE3;

    int verbose;
    int printVersion;
} SeparationFlags;

/* Process priority to use for GPU version */
#ifndef _WIN32
  #define DEFAULT_GPU_PRIORITY 0
#else
  #define DEFAULT_GPU_PRIORITY MW_PRIORITY_NORMAL
#endif /* _WIN32 */

#define DEFAULT_POLLING_MODE 1
#define DEFAULT_TARGET_FREQUENCY 30.0
#define DEFAULT_DISABLE_GPU_CHECKPOINTING 0


#define EMPTY_SEPARATION_FLAGS { NULL, NULL, NULL, NULL, NULL, NULL,           \
                                 0, FALSE, FALSE, FALSE, 0, FALSE, FALSE,      \
                                 UINT_MAX, 0, FALSE,                           \
                                 DEFAULT_TARGET_FREQUENCY,                     \
                                 DEFAULT_POLLING_MODE,                         \
                                 DEFAULT_DISABLE_GPU_CHECKPOINTING,            \
                                 0, FALSE,                                     \
                                 FALSE, FALSE, FALSE, FALSE, FALSE, FALSE      \
                               }

#ifdef __cplusplus
}
#endif


#endif /* _SEPARATION_H_ */

