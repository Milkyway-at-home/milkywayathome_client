/*
 * Copyright (c) 2011 Matthew Arsenault
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

#ifndef _NBODY_DEFAULTS_H_
#define _NBODY_DEFAULTS_H_

#include "nbody_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DEFAULT_CHECKPOINT_FILE "nbody_checkpoint"
#define DEFAULT_HISTOGRAM_FILE  "histogram"

#define DEFAULT_LIKELIHOOD_METHOD NBODY_EMD

/* 15 minutes */
#define NOBOINC_DEFAULT_CHECKPOINT_PERIOD 900


#define DEFAULT_SUN_GC_DISTANCE ((real) 8.0)
#define DEFAULT_CRITERION NewCriterion
#define DEFAULT_TREE_ROOT_SIZE ((real) 4.0)

#define DEFAULT_USE_QUADRUPOLE_MOMENTS TRUE
#define DEFAULT_ALLOW_INCEST FALSE
#define DEFAULT_QUIET_ERRORS FALSE

  /*
    Return this when a big likelihood is needed.
  */
#define DEFAULT_WORST_CASE ((real) 9999999.9)
 /*
    Return this when a small likelihood is needed.
  */
#define DEFAULT_BEST_CASE ((real) 1e-9)

#define histogramPhi 128.79
#define histogramTheta 54.39
#define histogramPsi 90.70
#define histogramlambdaStart ((real) -50.0)
#define histogramlambdaEnd ((real) 50.0)
#define histogramlambdaBins ((unsigned int) 34)
#define histogrambetaStart ((real) -25.0)
#define histogrambetaEnd ((real) 25.0)
#define histogrambetaBins ((unsigned int) 10)


extern const NBodyCtx defaultNBodyCtx;
extern const HistogramParams defaultHistogramParams;

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_DEFAULTS_H_ */

