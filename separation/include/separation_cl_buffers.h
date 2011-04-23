/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#ifndef _SEPARATION_CL_BUFFERS_H_
#define _SEPARATION_CL_BUFFERS_H_

#include "milkyway_cl.h"
#include "mw_cl.h"
#include "separation_types.h"
#include "setup_cl.h"


#ifdef __cplusplus
extern "C" {
#endif

cl_int createSeparationBuffers(CLInfo* ci,
                               SeparationCLMem* cm,
                               const AstronomyParameters* ap,
                               const IntegralArea* ia,
                               const StreamConstants* sc,
                               const StreamGauss sg,
                               const SeparationSizes* sizes,
                               cl_bool useImages);

void releaseSeparationBuffers(SeparationCLMem* cm);

void calculateSizes(SeparationSizes* sizes, const AstronomyParameters* ap, const IntegralArea* ia);

real* mapIntegralResults(CLInfo* ci, SeparationCLMem* cm, size_t resultsSize);

real* mapStreamsResults(CLInfo* ci, SeparationCLMem* cm, size_t streamsResultsSize);


#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_CL_BUFFERS_H_ */

