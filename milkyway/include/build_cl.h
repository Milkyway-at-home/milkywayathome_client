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

#ifndef _BUILD_CL_H_
#define _BUILD_CL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <OpenCL/cl.h>
#include <OpenCL/cl_ext.h>

typedef struct
{
    cl_device_id dev;
    cl_device_type devType;
    unsigned int devCount;
    cl_context clctx;
    cl_command_queue queue;
    cl_program prog;
    cl_kernel kern;
} CLInfo;

#define EMPTY_CL_INFO { -1, -1, 0, NULL, NULL, NULL, NULL }


cl_int getCLInfo(CLInfo* ci,
                 cl_device_type type,
                 const char* kernName,
                 const char** src,
                 const char* compileDefs);

void destroyCLInfo(CLInfo* ci);
cl_int printCLExtensions(cl_device_id dev);

#ifdef __cplusplus
}
#endif

#endif /* _BUILD_CL_H_ */

