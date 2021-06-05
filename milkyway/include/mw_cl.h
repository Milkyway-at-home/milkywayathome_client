/*
 *  Copyright (c) 2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MW_CL_H_
#define _MW_CL_H_

#include "milkyway_config.h"

#ifdef __APPLE__
  #include <OpenCL/cl.h>
  #include <OpenCL/cl_platform.h>
  #include <OpenCL/cl_ext.h>
#else
/* FIXME: Sometimes when building the Windows executable, we need to use the full pathway.
   Would be nice if this were automated or if the setup was explained better in README.    */
  #include "CL/cl.h"
  #include "CL/cl_platform.h"
  #include "CL/cl_ext.h"

//  #include "/usr/include/CL/cl.h"
//  #include "/usr/include/CL/cl_platform.h"
//  #include "/usr/include/CL/cl_ext.h"
#endif /* __APPLE__ */

/* This doesn't seem to exist on OS X, but the callback on ATI on

/* This doesn't seem to exist on OS X, but the callback on ATI on
 * Linux/Windows dies without it */
#ifndef CL_CALLBACK
  #define CL_CALLBACK
#endif


#ifndef CL_PLATFORM_NOT_FOUND_KHR
  #define CL_PLATFORM_NOT_FOUND_KHR  -1001
#endif

#ifndef CL_MEM_USE_PERSISTENT_MEM_AMD
  #define CL_MEM_USE_PERSISTENT_MEM_AMD (1 << 6)
#endif

#ifndef CL_DEVICE_BOARD_NAME_AMD
  #define CL_DEVICE_BOARD_NAME_AMD 0x4038
#endif


#endif /* _MW_CL_H_ */

