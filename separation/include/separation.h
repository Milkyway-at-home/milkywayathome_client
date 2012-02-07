/*
 *  Copyright (c) 2008-2010 Travis Desell, Nathan Cole
 *  Copyright (c) 2008-2010 Boleslaw Szymanski, Heidi Newberg
 *  Copyright (c) 2008-2010 Carlos Varela, Malik Magdon-Ismail
 *  Copyright (c) 2008-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
  #include "milkyway_cl_util.h"
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
    int usePlatform;
    int useDevNumber;  /* Choose CL platform and device */
    int magicFactor;
    int nonResponsive;
    double targetFrequency;
    double waitFactor;  /* When using high CPU CL workarounds, factor for initial wait */
    int pollingMode;
    int disableGPUCheckpointing;

    MWPriority processPriority;

    int forceNoOpenCL;
    int forceNoILKernel;

    /* Force between normal, SSE2, SSE3 paths */
    int forceNoIntrinsics;
    int forceX87;
    int forceSSE2;
    int forceSSE3;
    int forceSSE41;
    int forceAVX;

    int verbose;
} SeparationFlags;

#define DEFAULT_GPU_PRIORITY MW_PRIORITY_NORMAL
#define DEFAULT_NON_RESPONSIVE FALSE
#define DEFAULT_TARGET_FREQUENCY 60.0
#define DEFAULT_WAIT_FACTOR 0.75
#define DEFAULT_DISABLE_GPU_CHECKPOINTING FALSE
#define DEFAULT_DISABLE_OPENCL FALSE
#define DEFAULT_DISABLE_IL_KERNEL FALSE

#if SEPARATION_OPENCL
  #define DEFAULT_POLLING_MODE MW_POLL_WORKAROUND_CL_WAIT_FOR_EVENTS
#else
  #define DEFAULT_POLLING_MODE 0
#endif


#ifdef __cplusplus
}
#endif


#endif /* _SEPARATION_H_ */

