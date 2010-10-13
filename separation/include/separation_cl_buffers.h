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

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_cl.h"
#include "separation_types.h"
#include "setup_cl.h"
#include "build_cl.h"

typedef struct
{
    size_t outMu;
    size_t outProbs;

    size_t ap;        /* Constants */
    size_t sc;
    size_t ia;
    size_t ic;
    size_t rc;
    size_t rPts;
} SeparationSizes;


cl_int createSeparationBuffers(CLInfo* ci,
                               SeparationCLMem* cm,
                               const ASTRONOMY_PARAMETERS* ap,
                               const INTEGRAL_AREA* ia,
                               const STREAM_CONSTANTS* sc,
                               const STREAM_GAUSS* sg,
                               const SeparationSizes* sizes);

void releaseSeparationBuffers(SeparationCLMem* cm);

void swapOutputBuffers(SeparationCLMem* cm);
cl_int separationSetOutputBuffers(CLInfo* ci, SeparationCLMem* cm);
void calculateSizes(SeparationSizes* sizes, const ASTRONOMY_PARAMETERS* ap, const INTEGRAL_AREA* ia);


#ifdef __cplusplus
}
#endif

#endif /* _SEPARATION_CL_BUFFERS_H_ */

