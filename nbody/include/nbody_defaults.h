/*
Copyright (C) 2011  Matthew Arsenault

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

#ifndef _NBODY_DEFAULTS_H_
#define _NBODY_DEFAULTS_H_

#include "nbody_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DEFAULT_CHECKPOINT_FILE "nbody_checkpoint"
#define DEFAULT_HISTOGRAM_FILE  "histogram"

#define DEFAULT_LIKELIHOOD_METHOD NBODY_ORIG_CHISQ

/* 15 minutes */
#define NOBOINC_DEFAULT_CHECKPOINT_PERIOD 900


#define DEFAULT_SUN_GC_DISTANCE ((real) 8.0)
#define DEFAULT_CRITERION NewCriterion
#define DEFAULT_TREE_ROOT_SIZE ((real) 4.0)

#define DEFAULT_USE_QUADRUPOLE_MOMENTS TRUE
#define DEFAULT_ALLOW_INCEST FALSE
#define DEFAULT_QUIET_ERRORS FALSE

#define histogramPhi 128.79
#define histogramTheta 54.39
#define histogramPsi 90.70
#define histogramStartRaw ((real) -50.0)
#define histogramEndRaw ((real) 50.0)
#define histogramBinSize ((real) 2.9411764705882355)
#define histogramCenter ((real) 0.0)

#define N_ORBIT_TRACE_POINTS 50


extern const NBodyCtx defaultNBodyCtx;
extern const HistogramParams defaultHistogramParams;

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_DEFAULTS_H_ */

