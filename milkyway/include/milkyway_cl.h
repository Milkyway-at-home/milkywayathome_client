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

#ifndef _MILKYWAY_CL_H_
#define _MILKYWAY_CL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_config.h"

#if !defined(__OPENCL_VERSION__) && MILKYWAY_OPENCL
  #ifdef __APPLE__
    #include <OpenCL/cl.h>
    #include <OpenCL/cl_platform.h>
    #include <OpenCL/cl_ext.h>
  #else
    #include <CL/cl.h>
    #include <CL/cl_platform.h>
    #include <CL/cl_ext.h>
  #endif /* __APPLE__ */
#endif  /* !defined(__OPENCL_VERSION__) && MILKYWAY_OPENCL */

#ifdef __OPENCL_VERSION__
  #define __MW_LOCAL __local
  #define __MW_PRIVATE __private
  #define __MW_GLOBAL __global
  #define __MW_CONSTANT __constant
#else
  #define __MW_LOCAL
  #define __MW_PRIVATE
  #define __MW_GLOBAL
  #define __MW_CONSTANT const
#endif /* __OPENCL_VERSION__ */

/* The ATI CL compiler is horribly broken, especially when using
 * doubles on the GPU. */
#if defined(__ATI_CL__) && !defined(__CPU__)
  #define BROKEN_CL_MATH 1
#else
  #define BROKEN_CL_MATH 0
#endif

#if defined(__ATI_RV770__) || defined(__ATI_RV730__) || defined(__ATI_RV710__)
  #define MW_RADEON_4XXX 1
#endif

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_CL_H_ */

