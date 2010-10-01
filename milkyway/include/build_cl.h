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

#include "milkyway_cl.h"

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

#define EMPTY_CL_INFO { NULL, -1, 0, NULL, NULL, NULL, NULL }

typedef struct
{
    size_t wgs;      /* CL_KERNEL_WORK_GROUP_SIZE */
    size_t cwgs[3];  /* CL_KERNEL_COMPILE_WORK_GROUP_SIZE */
    cl_ulong lms;    /* CL_KERNEL_LOCAL_MEM_SIZE */
} WGInfo;

#define EMPTY_WG_INFO { 0, { 0, 0, 0 }, 0 }

typedef struct
{
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
    char devName[128];
    char vendor[128];
    char version[128];
    char driver[128];
    //char clCVer[128];

    size_t maxWorkItemSizes[3];
} DevInfo;

cl_int mwSetupCL(CLInfo* ci,
                 cl_device_type type,
                 const char* kernName,
                 const char** src,
                 const cl_uint srcCount,
                 const char* compileDefs);

cl_int destroyCLInfo(CLInfo* ci);
cl_int printCLExtensions(cl_device_id dev);
cl_int getWorkGroupInfo(CLInfo* ci, WGInfo* wgi);
void printWorkGroupInfo(const WGInfo* wgi);

cl_int getDevInfo(DevInfo* di, cl_device_id dev);
void printDevInfo(const DevInfo* di);

#ifdef __cplusplus
}
#endif

#endif /* _BUILD_CL_H_ */

