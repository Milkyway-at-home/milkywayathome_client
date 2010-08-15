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

#include "separation.h"
#include <OpenCL/cl.h>
#include <OpenCL/cl_ext.h>


    /* Host OpenCL stuff */
typedef struct
{
    cl_device_id dev;
    cl_device_type devType;
    unsigned int devCount;
    cl_context clctx;
    cl_command_queue queue;
    cl_program prog;
    cl_kernel kern;
} SeparationCLInfo;

#define EMPTY_SEPARATION_CL_INFO { -1, -1, 0, NULL, NULL, NULL, NULL }

/* The various buffers needed by the integrate function. */
typedef struct
{
    /* Write only buffers */
    cl_mem outNu;     /* Output from each nu_sum done in parallel */

    /* constant, read only buffers */
    cl_mem ap;        /* Astronomy parameters */
    cl_mem sc;        /* Stream Constants */
    cl_mem rConsts;   /* r step constants */
    cl_mem rPoints;   /* r points */
    cl_mem nuConsts;  /* nu step constants */
    cl_mem ia;        /* Integral areas */
} SeparationCLMem;

#define EMPTY_SEPARATION_CL_MEM { NULL, NULL, NULL, NULL, NULL, NULL, NULL }

int setupSeparationCL(const ASTRONOMY_PARAMETERS* ap,
                      const STREAM_CONSTANTS* sc,
                      const R_CONSTANTS* r_consts,
                      const R_POINTS* r_points,
                      const NU_CONSTANTS* nu_st,
                      const INTEGRAL_AREA* ia);

#ifdef __cplusplus
}
#endif

#endif /* _SETUP_CL_H_ */

