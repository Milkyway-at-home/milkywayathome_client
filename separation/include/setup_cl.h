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

#ifndef _SETUP_CL_H_
#define _SETUP_CL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_cl.h"
#include "separation_types.h"
#include "build_cl.h"

/* The various buffers needed by the integrate function. */
typedef struct
{
    /* Write only buffers */
    cl_mem outMu;     /* Output from each mu_sum done in parallel */
    cl_mem outProbs;

    /* constant, read only buffers */
    cl_mem ap;
    cl_mem ia;
    cl_mem sc;        /* Stream Constants */
    cl_mem rPts;
} SeparationCLMem;

#define EMPTY_SEPARATION_CL_MEM { NULL, NULL, NULL, NULL, NULL, NULL }


cl_int setupSeparationCL(const ASTRONOMY_PARAMETERS* ap,
                         const INTEGRAL_AREA* ia,
                         const STREAM_CONSTANTS* sc,
                         const R_POINTS* r_pts_all,
                         CLInfo* ci,
                         SeparationCLMem* cm);

#ifdef __cplusplus
}
#endif

#endif /* _SETUP_CL_H_ */

