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

#ifdef _cplusplus
extern "C" {
#endif

#include <json/json.h>
#include "nbody_config.h"
#include "nbody_util.h"


#if !DYNAMIC_PRECISION
  __attribute__ ((visibility("default"))) void runNBodySimulation(json_object* obj,
                                                                  const char* outFileName,
                                                                  const char* checkpointFileName,
                                                                  const char* histogramFileName,
                                                                  const char* histoutFileName,
                                                                  const int outputCartesian,
                                                                  const int printTiming,
                                                                  const int verifyOnly);
#else
  __attribute__ ((visibility("default"))) void runNBodySimulation_float(json_object* obj,
                                                                        const char* outFileName,
                                                                        const char* checkpointFileName,
                                                                        const char* histogramFileName,
                                                                        const char* histoutFileName,
                                                                        const int outputCartesian,
                                                                        const int printTiming,
                                                                        const int verifyOnly);

  __attribute__ ((visibility("default"))) void runNBodySimulation_double(json_object* obj,
                                                                         const char* outFileName,
                                                                         const char* checkpointFileName,
                                                                         const char* histogramFileName,
                                                                         const char* histoutFileName,
                                                                         const int outputCartesian,
                                                                         const int printTiming,
                                                                         const int verifyOnly);
#endif /* DYNAMIC_PRECISION */

#ifdef _cplusplus
}
#endif

#endif /* _NBODY_H_ */

