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

#ifndef _MW_CL_TYPES_H_
#define _MW_CL_TYPES_H_

#include "mw_cl.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MW_CL_ERROR ((cl_int) 1)


typedef enum
{
    MW_NONE_DOUBLE = 0,
    MW_CL_AMD_FP64 = 1 << 1,
    MW_CL_KHR_FP64 = 1 << 2
} MWDoubleExts;

typedef enum
{
    MW_AMD_ATI = 0x1002,
    MW_NVIDIA = 0x10de
} MW_VENDOR_ID;

typedef struct
{
    cl_device_id devID;
    cl_device_type devType;
    cl_uint maxCompUnits, clockFreq;
    cl_ulong memSize;
    cl_ulong gMemCache;
    cl_ulong localMemSize;
    cl_device_local_mem_type localMemType;
    //cl_bool unifiedMem;
    cl_uint maxConstArgs;
    cl_ulong maxConstBufSize;
    cl_ulong maxMemAlloc;
    size_t maxWorkGroupSize;
    size_t maxParamSize;
    size_t timerRes;
    cl_uint maxWorkItemDim;
    cl_uint memBaseAddrAlign;
    cl_uint minAlignSize;
    cl_uint vendorID;
    cl_uint addrBits;
    cl_uint cachelineSize;
    cl_bool littleEndian;
    cl_bool errCorrect;
    cl_bool imgSupport;
    cl_bool nonOutput;  /* A check for this should maybe be in OpenCL. Something like a Tesla GPU without video output */
    char devName[128];
    char vendor[128];
    char version[128];
    char driver[128];
    cl_uint computeCapabilityMajor; /* Nvidia only */
    cl_uint computeCapabilityMinor;

    //char clCVer[128];

    size_t maxWorkItemSizes[3];
    MWDoubleExts doubleExts;
    char exts[1024];
} DevInfo;


typedef struct
{
    cl_device_id dev;
    cl_device_type devType;
    cl_uint devCount;
    cl_context clctx;
    cl_command_queue queue;
    cl_command_queue bufQueue; /* Queue for buffer ops when double buffering */
    cl_program prog;
    cl_kernel kern;

    DevInfo di;
} CLInfo;


typedef struct
{
    char name[128];
    char vendor[128];
    char version[128];
    char extensions[512];
    char profile[128];
} PlatformInfo;

#define EMPTY_PLATFORM_INFO { "", "", "", "", "" }


typedef struct
{
    size_t wgs;      /* CL_KERNEL_WORK_GROUP_SIZE */
    size_t cwgs[3];  /* CL_KERNEL_COMPILE_WORK_GROUP_SIZE */
    cl_ulong lms;    /* CL_KERNEL_LOCAL_MEM_SIZE */
} WGInfo;

#define EMPTY_WG_INFO { 0, { 0, 0, 0 }, 0 }



#ifdef __cplusplus
}
#endif

#endif /* _MW_CL_TYPES_H_ */

