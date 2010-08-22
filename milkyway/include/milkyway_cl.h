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

#ifndef __OPENCL_VERSION__
  #ifdef __APPLE__
    #include <OpenCL/cl.h>
    #include <OpenCL/cl_ext.h>
  #else
    #include <CL/cl.h>
    #include <CL/cl_ext.h>
  #endif /* __APPLE__ */
#endif  /* __OPENCL_VERSION__ */

#ifdef __OPENCL_VERSION__
  #define __MW_LOCAL __local
  #define __MW_PRIVATE __private
  #define __MW_GLOBAL __global
  #define __MW_CONSTANT __constant
#else
  #define __MW_LOCAL
  #define __MW_PRIVATE
  #define __MW_GLOBAL
  #define __MW_CONSTANT
#endif /* __OPENCL_VERSION__ */

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_CL_H_ */

