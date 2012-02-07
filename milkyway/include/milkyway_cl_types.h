/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
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

typedef enum MWCALtargetEnum {
    MW_CAL_TARGET_UNKNOWN   = -2,
    MW_CAL_TARGET_INVALID   = -1,
    MW_CAL_TARGET_600       = 0,      /**< R600 GPU ISA */
    MW_CAL_TARGET_610       = 1,      /**< RV610 GPU ISA */
    MW_CAL_TARGET_630       = 2,      /**< RV630 GPU ISA */
    MW_CAL_TARGET_670       = 3,      /**< RV670 GPU ISA */
    MW_CAL_TARGET_7XX       = 4,      /**< R700 class GPU ISA */
    MW_CAL_TARGET_770       = 5,      /**< RV770 GPU ISA */
    MW_CAL_TARGET_710       = 6,      /**< RV710 GPU ISA */
    MW_CAL_TARGET_730       = 7,      /**< RV730 GPU ISA */
    MW_CAL_TARGET_CYPRESS   = 8,      /**< CYPRESS GPU ISA */
    MW_CAL_TARGET_JUNIPER   = 9,      /**< JUNIPER GPU ISA */
    MW_CAL_TARGET_REDWOOD   = 10,     /**< REDWOOD GPU ISA */
    MW_CAL_TARGET_CEDAR     = 11,     /**< CEDAR GPU ISA */

    MW_CAL_TARGET_SUMO      = 12,
    MW_CAL_TARGET_SUPERSUMO = 13,

    MW_CAL_TARGET_WRESTLER  = 14,     /**< WRESTLER GPU ISA */
    MW_CAL_TARGET_CAYMAN    = 15,     /**< CAYMAN GPU ISA */
    MW_CAL_TARGET_RESERVED2 = 16,
    MW_CAL_TARGET_BARTS     = 17,     /**< BARTS GPU ISA */

    MW_CAL_TARGET_TURKS     = 18,
    MW_CAL_TARGET_CAICOS    = 19,

    MW_CAL_TARGET_TAHITI    = 20,
    MW_CAL_TARGET_THAMES    = 21,
    MW_CAL_TARGET_LOMBOK    = 22
} MWCALtargetEnum;

typedef struct
{
    cl_device_id devID;
    cl_device_type devType;
    cl_uint maxCompUnits, clockFreq;
    cl_ulong memSize;
    cl_ulong gMemCache;
    cl_ulong localMemSize;
    cl_device_local_mem_type localMemType;
    cl_device_fp_config doubleFPConfig;
    cl_device_fp_config floatFPConfig;
    //cl_bool unifiedMem;
    cl_uint maxConstArgs;
    cl_ulong maxConstBufSize;
    cl_ulong maxMemAlloc;
    size_t maxWorkGroupSize;
    size_t maxParamSize;
    size_t timerRes;
    cl_uint warpSize;
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
    cl_bool hasGraphicsQOS;
    char devName[128];
    char vendor[128];
    char version[128];
    char driver[128];
    cl_uint computeCapabilityMajor; /* Nvidia only */
    cl_uint computeCapabilityMinor;
    MWCALtargetEnum calTarget;       /* AMD Only */

    cl_uint doubleFrac; /* Estimated speed of doubles relative to float */
    cl_uint aluPerCU;

    size_t maxWorkItemSizes[3];
    MWDoubleExts doubleExts;
    char exts[1024];
} DevInfo;


typedef struct
{
    cl_platform_id plat;
    cl_device_id dev;
    cl_device_type devType;
    cl_uint devCount;
    cl_context clctx;
    cl_command_queue queue;
    cl_command_queue bufQueue; /* Queue for buffer ops when double buffering */

    cl_int pollingMode;  /* Hint for how to poll for completion of event */

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

