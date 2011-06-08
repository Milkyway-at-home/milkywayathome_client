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

#if !defined(_MILKYWAY_CL_H_INSIDE_) && !defined(MILKYWAY_CL_COMPILATION)
  #error "Only milkyway_cl.h can be included directly."
#endif


#ifndef _MW_CL_SHOW_TYPES_H_
#define _MW_CL_SHOW_TYPES_H_

#include "mw_cl_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: clint and memflags are usually or'd, so most of the time these won't work right */
const char* showCLDeviceType(const cl_device_type x) CONST_F;
const char* showCLBuildStatus(const cl_build_status x) CONST_F;
const char* showCLInt(const cl_int x) CONST_F;
const char* showCLMemFlags(const cl_mem_flags x) CONST_F;
const char* showCLDeviceFPConfig(const cl_device_fp_config x) CONST_F;
const char* showCLDeviceLocalMemType(const cl_device_local_mem_type x) CONST_F;
const char* showCLDeviceExecCapabilities(const cl_device_exec_capabilities x) CONST_F;
const char* showCLCommandQueueProperties(const cl_command_queue_properties x) CONST_F;
const char* showCLBool(const cl_bool x) CONST_F;
const char* showCLDeviceMemCacheType(const cl_device_mem_cache_type x) CONST_F;
const char* showCLKernelInfo(const cl_kernel_info x) CONST_F;
const char* showMWDoubleExts(const MWDoubleExts x)  CONST_F;

#ifdef __cplusplus
}
#endif

#endif /* _MW_CL_SHOW_TYPES_H_ */

