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

#ifndef _SHOW_CL_TYPES_H_
#define _SHOW_CL_TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_cl.h"

/* TODO: clint and memflags are usually or'd, so most of the time these won't work right */
const char* showCLDeviceType(const cl_device_type x) __attribute__ ((const));
const char* showCLBuildStatus(const cl_build_status x) __attribute__ ((const));
const char* showCLInt(const cl_int x) __attribute__ ((const));
const char* showCLMemFlags(const cl_mem_flags x) __attribute__ ((const));

#ifdef __cplusplus
}
#endif

#endif /* _SHOW_CL_TYPES_H_ */

