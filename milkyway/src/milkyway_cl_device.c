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

#include "milkyway_cl_device.h"
#include "milkyway_cl_util.h"
#include "milkyway_util.h"
#include "milkyway_cl_show_types.h"

/* These are missing from the current OS X headers */
#ifndef CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
  #define CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV 0x4000
#endif
#ifndef CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV
  #define CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV 0x4001
#endif
#ifndef CL_DEVICE_WARP_SIZE_NV
  #define CL_DEVICE_WARP_SIZE_NV 0x4003
#endif

typedef struct
{
    const char* name;
    MWCALtargetEnum target;
    cl_uint aluPerCU; /* Needed for calculating estimated peak flops */
    cl_uint doubleFrac;
    cl_uint wavefrontSize;
} AMDGPUData;

/* A reasonable guess */
static const AMDGPUData invalidAMDGPUData = { "Invalid", MW_CAL_TARGET_INVALID, 5, 5, 64 };

/* Table of more detailed info we want that OpenCL won't give us */
static const AMDGPUData amdGPUData[] =
{
     /* R600 */
    { "ATI RV670",  MW_CAL_TARGET_670,      16 * 5, 0, 64 },
    { "ATI RV630",  MW_CAL_TARGET_630,      16 * 5, 0, 32 },
    { "ATI RV610",  MW_CAL_TARGET_610,      16 * 5, 0, 32 },
    { "ATI RV600",  MW_CAL_TARGET_600,      16 * 5, 0, 16 },

    /* R700 */
    { "ATI RV770",  MW_CAL_TARGET_770,      16 * 5, 5, 64 },
    { "ATI RV730",  MW_CAL_TARGET_730,      16 * 5, 0, 32 },
    { "ATI RV710",  MW_CAL_TARGET_710,      16 * 5, 0, 16 },

    /* Evergreen */
    { "Cypress",    MW_CAL_TARGET_CYPRESS,  16 * 5, 5, 64 },
    { "Juniper",    MW_CAL_TARGET_JUNIPER,  16 * 5, 0, 64 },
    { "Redwood",    MW_CAL_TARGET_REDWOOD,  16 * 5, 0, 64 },
    { "Cedar",      MW_CAL_TARGET_CEDAR,    16 * 5, 0, 32 },

    /* Northern Islands */
    { "Cayman",     MW_CAL_TARGET_CAYMAN,   16 * 4, 4, 64 },
    { "Barts",      MW_CAL_TARGET_BARTS,    16 * 5, 0, 64 },
    { "Turks",      MW_CAL_TARGET_TURKS,    16 * 5, 0, 64 },
    { "Caicos",     MW_CAL_TARGET_CAICOS,   16 * 5, 0, 64 },

    /* Southern Islands */
    { "Tahiti",     MW_CAL_TARGET_TAHITI,   4 * 16, 4, 64 },
    { "Thames",     MW_CAL_TARGET_THAMES,   16 * 4, 4, 64 },
    { "Lombok",     MW_CAL_TARGET_LOMBOK,   16 * 4, 4, 64 },

#if 0
    /* These are there, but I don't know about them */
    { "WinterPark",                        16 * 5, 0, 64 },
    { "BeaverCreek",                       16 * 5, 0, 64 },
    { "Loveland",                          16 * 5, 0, 64 },

    { "Lions",                             16 * 5, 0, 64 },
    { "Tigers",                            16 * 5, 0, 64 },
    { "Bears",                             16 * 5, 0, 64 },
#endif

    { NULL,          MW_CAL_TARGET_INVALID, 4 * 16, 0, 64 }
};

cl_int mwUAVIdFromMWCALtargetEnum(MWCALtargetEnum x)
{
    return (x >= MW_CAL_TARGET_CYPRESS) ? 11 : 1;
}

static const AMDGPUData* mwLookupAMDGPUInfo(const DevInfo* di)
{
    const AMDGPUData* p = amdGPUData;

    while (p->name)
    {
        if (!strncasecmp(di->devName, p->name, sizeof(di->devName)))
        {
            return p;
        }

        ++p;
    }

    return &invalidAMDGPUData;
}

cl_double mwAMDEstimateGFLOPs(const DevInfo* di, cl_bool useDouble)
{
    cl_ulong flops, flopsFloat, flopsDouble;
    cl_double gflops;

    flopsFloat = 2 * (di->maxCompUnits * di->aluPerCU) * (cl_ulong) di->clockFreq * 1000000;
    flopsDouble = flopsFloat / di->doubleFrac;

    mw_printf("Estimated AMD GPU GFLOP/s: %.0f SP GFLOP/s, %.0f DP FLOP/s\n",
              1.0e-9 * (cl_double) flopsFloat,
              1.0e-9 * (cl_double) flopsDouble);

    flops = useDouble ? flopsDouble : flopsFloat;

    gflops = floor(1.0e-9 * (cl_double) flops);

    /* At different times the AMD drivers have reported 0 as the clock
     * speed, so try to catch that. We could test the GPU and figure
     * out what the FLOPs should be to get a better estimate.
     */
    if (gflops <= 100.0)
    {
        mw_printf("Warning: Bizarrely low flops (%.0f). Defaulting to %.0f\n", gflops, 100.0);
        gflops = 100.0;
    }

    return gflops;
}

cl_bool mwHasNvidiaCompilerFlags(const DevInfo* di)
{
    return strstr(di->exts, "cl_nv_compiler_options") != NULL;
}

cl_bool mwDeviceVendorIsAMD(const DevInfo* di)
{
    return (strncmp(di->vendor, "Advanced Micro Devices, Inc.", sizeof(di->vendor)) == 0);
}

cl_bool mwDeviceVendorIsNvidia(const DevInfo* di)
{
    return (strncmp(di->vendor, "Nvidia Corporation", sizeof(di->vendor)) == 0);
}

/* True if devices compute capability >= requested version */
cl_bool mwMinComputeCapabilityCheck(const DevInfo* di, cl_uint major, cl_uint minor)
{
    return    di->computeCapabilityMajor > major
           || (   di->computeCapabilityMajor == major
               && di->computeCapabilityMinor >= minor);
}

/* Exact check on compute capability version */
cl_bool mwComputeCapabilityIs(const DevInfo* di, cl_uint major, cl_uint minor)
{
    return di->computeCapabilityMajor == major && di->computeCapabilityMinor == minor;
}

/* approximate ratio of float : double flops */
static cl_uint mwCUDAEstimateDoubleFrac(const DevInfo* di)
{
    /* FIXME: This also differs with generation.
       Is there a better way to find out the generation and if
     */

    if (mwMinComputeCapabilityCheck(di, 2, 0))
    {
        if (strstr(di->devName, "Tesla") != NULL)
        {
            return 2;
        }
        else
        {
            return 8;
        }
    }
    else
    {
        if (strstr(di->devName, "Tesla") != NULL)
        {
            return 2;
        }
        else
        {
            return 8;
        }
    }
}

/* Different on different Nvidia architectures.
   Uses numbers from appendix of Nvidia OpenCL programming guide. */
static cl_uint mwCUDACoresPerComputeUnit(const DevInfo* di)
{
    if (mwMinComputeCapabilityCheck(di, 2, 0))
        return 32;

    return 8;     /* 1.x is 8 */
}

cl_double mwCUDAEstimateGFLOPs(const DevInfo* di, cl_bool useDouble)
{
    cl_ulong flopsFloat, flopsDouble, flops;
    cl_double gflops;

    flopsFloat = 2 * di->maxCompUnits * di->aluPerCU * (cl_ulong) di->clockFreq * 1000000;
    flopsDouble = flopsFloat / di->doubleFrac;

    mw_printf("Estimated Nvidia GPU GFLOP/s: %.0f SP GFLOP/s, %.0f DP FLOP/s\n",
              1.0e-9 * (cl_double) flopsFloat, 1.0e-9 * (cl_double) flopsDouble);

    flops = useDouble ? flopsDouble : flopsFloat;
    gflops = 1.0e-9 * (cl_double) flops;


    if (gflops <= 50.0)
    {
        mw_printf("Warning: Bizarrely low flops (%.0f). Defaulting to %.0f\n", gflops, 50.0);
        gflops = 50.0;
    }

    return gflops;
}

cl_bool mwIsNvidiaGPUDevice(const DevInfo* di)
{
    return (di->vendorID == MW_NVIDIA) && (di->devType == CL_DEVICE_TYPE_GPU);
}

cl_bool mwIsAMDGPUDevice(const DevInfo* di)
{
    /* Not sure if the vendor ID for AMD is the same with their
       CPUs.  Also something else weird was going on with the
       vendor ID, so check the name just in case.
    */

    return (di->vendorID == MW_AMD_ATI || mwDeviceVendorIsAMD(di));
}

cl_double mwDeviceEstimateGFLOPs(const DevInfo* di, cl_bool useDouble)
{
    cl_double gflops = 0.0;

    if (di->devType == CL_DEVICE_TYPE_GPU)
    {
        if (mwIsNvidiaGPUDevice(di))
        {
            gflops = mwCUDAEstimateGFLOPs(di, useDouble);
        }
        else if (mwIsAMDGPUDevice(di))
        {
            gflops = mwAMDEstimateGFLOPs(di, useDouble);
        }
        else
        {
            mw_printf("Unhandled GPU vendor '%s' (0x%x)\n", di->vendor, di->vendorID);
            gflops = 100.0;
        }
    }
    else
    {
        mw_printf("Missing flops estimate for device type %s\n", showCLDeviceType(di->devType));
        return 1.0;
    }

    return gflops;
}

/* Read the double supported extensions; i.e. AMD's subset or the actual Khronos one. */
MWDoubleExts mwGetDoubleExts(const char* exts)
{
    MWDoubleExts found = MW_NONE_DOUBLE;

    if (strstr(exts, "cl_amd_fp64"))
        found |= MW_CL_AMD_FP64;
    if (strstr(exts, "cl_khr_fp64"))
        found |= MW_CL_KHR_FP64;

    return found;
}

cl_bool mwSupportsDoubles(const DevInfo* di)
{
    return di->doubleExts != MW_NONE_DOUBLE;
}

cl_bool mwDeviceHasDenormals(const DevInfo* di, cl_bool doublePrec)
{
    cl_device_fp_config config = doublePrec ? di->doubleFPConfig : di->floatFPConfig;
    return !!(config & CL_FP_DENORM);
}

cl_bool mwDeviceHasFMA(const DevInfo* di, cl_bool doublePrec)
{
    cl_device_fp_config config = doublePrec ? di->doubleFPConfig : di->floatFPConfig;
    return !!(config & CL_FP_FMA);
}

cl_bool mwDeviceHasInfNan(const DevInfo* di, cl_bool doublePrec)
{
    cl_device_fp_config config = doublePrec ? di->doubleFPConfig : di->floatFPConfig;
    return !!(config & CL_FP_INF_NAN);
}

cl_bool mwDeviceHasRTN(const DevInfo* di, cl_bool doublePrec)
{
    cl_device_fp_config config = doublePrec ? di->doubleFPConfig : di->floatFPConfig;
    return !!(config & CL_FP_ROUND_TO_NEAREST);
}

static cl_bool mwDeviceIsNonOutput(const DevInfo* di)
{
    /* TODO: Correct way to test for Tesla type things?
       Is there a way we can find if a GPU is connected to a display in general?
     */
    return ((di->devType != CL_DEVICE_TYPE_GPU) || (strstr(di->devName, "Tesla") != NULL));
}

static cl_bool mwDeviceHasGraphicsQOS(const DevInfo* di)
{
    /* Tahiti has the capability but it hasn't been enabled in current drivers yet */
    return CL_FALSE;
}

cl_int mwGetDevInfo(DevInfo* di, cl_device_id dev)
{
    const AMDGPUData* amdData;
    cl_int err = CL_SUCCESS;

    di->devID = dev;

    err |= clGetDeviceInfo(dev, CL_DEVICE_TYPE,                     sizeof(di->devType),  &di->devType, NULL);

    err |= clGetDeviceInfo(dev, CL_DEVICE_NAME,                     sizeof(di->devName),  di->devName, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VENDOR,                   sizeof(di->vendor),   di->vendor, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VENDOR_ID,                sizeof(cl_uint),  &di->vendorID, NULL);
    err |= clGetDeviceInfo(dev, CL_DRIVER_VERSION,                  sizeof(di->driver),   di->driver, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VERSION,                  sizeof(di->version),  di->version, NULL);
  //err |= clGetDeviceInfo(dev, CL_DEVICE_OPENCL_C_VERSION,         sizeof(di->clCVer),   di->clCVer, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_ENDIAN_LITTLE,            sizeof(cl_bool),  &di->littleEndian, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(cl_bool),  &di->errCorrect, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_IMAGE_SUPPORT, sizeof(cl_bool),  &di->imgSupport, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_ADDRESS_BITS,             sizeof(cl_uint),  &di->addrBits, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS,        sizeof(cl_uint),  &di->maxCompUnits, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_CLOCK_FREQUENCY,      sizeof(cl_uint),  &di->clockFreq, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_SIZE,          sizeof(cl_ulong), &di->memSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_MEM_ALLOC_SIZE,       sizeof(cl_ulong), &di->maxMemAlloc, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,    sizeof(cl_ulong), &di->gMemCache, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, sizeof(cl_uint), &di->cachelineSize, NULL);

  //err |= clGetDeviceInfo(dev, CL_DEVICE_HOST_UNIFIED_MEMORY,      sizeof(cl_ulong), &unifiedMem, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(cl_device_local_mem_type), &di->localMemType, NULL);

    err |= clGetDeviceInfo(dev, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(cl_device_fp_config), &di->doubleFPConfig, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_SINGLE_FP_CONFIG, sizeof(cl_device_fp_config), &di->floatFPConfig, NULL);

    err |= clGetDeviceInfo(dev, CL_DEVICE_LOCAL_MEM_SIZE,           sizeof(cl_ulong), &di->localMemSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_CONSTANT_ARGS,        sizeof(cl_uint),  &di->maxConstArgs, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong), &di->maxConstBufSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_PARAMETER_SIZE, sizeof(size_t), &di->maxParamSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &di->maxWorkGroupSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &di->maxWorkItemDim, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(di->maxWorkItemSizes), di->maxWorkItemSizes, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MEM_BASE_ADDR_ALIGN, sizeof(cl_uint), &di->memBaseAddrAlign, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, sizeof(cl_uint), &di->minAlignSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_PROFILING_TIMER_RESOLUTION, sizeof(size_t), &di->timerRes, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_EXTENSIONS, sizeof(di->exts), &di->exts, NULL);

    di->computeCapabilityMajor = di->computeCapabilityMinor = 0;
    di->warpSize = 0;
    if (err == CL_SUCCESS)
    {
        if (strstr(di->exts, "cl_nv_device_attribute_query") != NULL)
        {
            err |= clGetDeviceInfo(dev, CL_DEVICE_WARP_SIZE_NV,
                                   sizeof(di->warpSize), &di->warpSize, NULL);
            err |= clGetDeviceInfo(di->devID, CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV,
                                   sizeof(cl_uint), &di->computeCapabilityMajor, NULL);
            err |= clGetDeviceInfo(di->devID, CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV,
                                   sizeof(cl_uint), &di->computeCapabilityMinor, NULL);
        }
        else
        {
            if (di->devType == CL_DEVICE_TYPE_CPU)
            {
                di->warpSize = 1;
            }
            else if (di->devType == CL_DEVICE_TYPE_GPU)
            {
                /* FIXME: How do I get this on AMD? It's 64 for all of
                 * the high end stuff, but 32 for lower. I think it's
                 * 64 for all the GPUs that do have doubles */
                di->warpSize = 64;
            }
            else
            {
                mw_printf("Unknown device type, using warp size = 1\n");
                di->warpSize = 1;
            }
        }
    }

    di->nonOutput = mwDeviceIsNonOutput(di);
    di->hasGraphicsQOS = mwDeviceHasGraphicsQOS(di);


    if (mwIsNvidiaGPUDevice(di))
    {
        di->aluPerCU = mwCUDACoresPerComputeUnit(di);
        di->doubleFrac = mwCUDAEstimateDoubleFrac(di);
        di->calTarget = MW_CAL_TARGET_INVALID;

        if (strstr(di->exts, "cl_nv_device_attribute_query") != NULL)
        {
            err |= clGetDeviceInfo(dev, CL_DEVICE_WARP_SIZE_NV,
                                   sizeof(di->warpSize), &di->warpSize, NULL);
            err |= clGetDeviceInfo(di->devID, CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV,
                                   sizeof(cl_uint), &di->computeCapabilityMajor, NULL);
            err |= clGetDeviceInfo(di->devID, CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV,
                                   sizeof(cl_uint), &di->computeCapabilityMinor, NULL);
        }
    }
    else if (mwIsAMDGPUDevice(di))
    {
        amdData = mwLookupAMDGPUInfo(di);

        di->aluPerCU   = amdData->aluPerCU;
        di->doubleFrac = amdData->doubleFrac;
        di->calTarget  = amdData->target;
        di->warpSize   = amdData->wavefrontSize;
    }

    if (di->warpSize == 0)
    {
        mw_printf("Unknown device type, using warp size = 1\n");
        di->warpSize = 1;
    }

    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error getting device information");
    }
    else
    {
        di->doubleExts = mwGetDoubleExts(di->exts);
    }

    return err;
}

void mwPrintDevInfo(const DevInfo* di)
{
    mw_printf("Device '%s' (%s:0x%x) (%s)\n"
              "Driver version:      %s\n"
              "Version:             %s\n"
              "Compute capability:  %u.%u\n"
              "Little endian:       %s\n"
              "Error correction:    %s\n"
              "Image support:       %s\n"
              "Address bits:        %u\n"
              "Max compute units:   %u\n"
              "Clock frequency:     %u Mhz\n"
              "Global mem size:     "LLU"\n"
              "Max mem alloc:       "LLU"\n"
              "Global mem cache:    "LLU"\n"
              "Cacheline size:      %u\n"
              "Local mem type:      %s\n"
              "Local mem size:      "LLU"\n"
              "Max const args:      %u\n"
              "Max const buf size:  "LLU"\n"
              "Max parameter size:  "ZU"\n"
              "Max work group size: "ZU"\n"
              "Max work item dim:   %u\n"
              "Max work item sizes: { "ZU", "ZU", "ZU" }\n"
              "Mem base addr align: %u\n"
              "Min type align size: %u\n"
              "Timer resolution:    "ZU" ns\n"
              "Warp size:           %u\n"
              "ALU per CU:          %u\n"
              "Double extension:    %s\n"
              "Double fraction:     1/%u\n"
              "FP config:     float  | double\n"
              "  FMA:       %8s  %8s\n"
              "  Denormals: %8s  %8s\n"
              "  RTN:       %8s  %8s\n"
              "  Inf/Nan:   %8s  %8s\n"
              "Extensions:\n"
              "  %s\n",
              di->devName,
              di->vendor,
              di->vendorID,
              showCLDeviceType(di->devType),
              di->driver,
              di->version,
              di->computeCapabilityMajor, di->computeCapabilityMinor,
              showCLBool(di->littleEndian),
              showCLBool(di->errCorrect),
              showCLBool(di->imgSupport),
              di->addrBits,
              di->maxCompUnits,
              di->clockFreq,
              di->memSize,
              di->maxMemAlloc,
              di->gMemCache,
              di->cachelineSize,
              //di->showBool(unifiedMem),
              showCLDeviceLocalMemType(di->localMemType),
              di->localMemSize,
              di->maxConstArgs,
              di->maxConstBufSize,
              di->maxParamSize,
              di->maxWorkGroupSize,
              di->maxWorkItemDim,
              di->maxWorkItemSizes[0], di->maxWorkItemSizes[1], di->maxWorkItemSizes[2],
              di->memBaseAddrAlign,
              di->minAlignSize,
              di->timerRes,
              di->warpSize,
              di->aluPerCU,
              showMWDoubleExts(di->doubleExts),
              di->doubleFrac,
              showCLBool(mwDeviceHasDenormals(di, CL_FALSE)), showCLBool(mwDeviceHasDenormals(di, CL_TRUE)),
              showCLBool(mwDeviceHasFMA(di, CL_FALSE)), showCLBool(mwDeviceHasFMA(di, CL_TRUE)),
              showCLBool(mwDeviceHasRTN(di, CL_FALSE)), showCLBool(mwDeviceHasRTN(di, CL_TRUE)),
              showCLBool(mwDeviceHasInfNan(di, CL_FALSE)), showCLBool(mwDeviceHasInfNan(di, CL_TRUE)),
              di->exts
        );
}

void mwPrintDevInfoShort(const DevInfo* di)
{
    mw_printf("Device '%s' (%s:0x%x) (%s)\n"
              "Driver version:      %s\n"
              "Version:             %s\n"
              "Compute capability:  %u.%u\n"
              "Max compute units:   %u\n"
              "Clock frequency:     %u Mhz\n"
              "Global mem size:     "LLU"\n"
              "Local mem size:      "LLU"\n"
              "Max const buf size:  "LLU"\n"
              "Double extension:    %s\n",
              di->devName,
              di->vendor, di->vendorID,
              showCLDeviceType(di->devType),
              di->driver,
              di->version,
              di->computeCapabilityMajor, di->computeCapabilityMinor,
              di->maxCompUnits,
              di->clockFreq,
              di->memSize,
              di->localMemSize,
              di->maxConstBufSize,
              showMWDoubleExts(di->doubleExts)
        );
}

static void mwGetPlatformInfo(PlatformInfo* pInfo, cl_platform_id platform)
{
    cl_int err;
    size_t readSize;

    err = clGetPlatformInfo(platform, CL_PLATFORM_NAME,
                            sizeof(pInfo->name), pInfo->name, &readSize);
    if (readSize > sizeof(pInfo->name))
        mwPerrorCL(err, "Failed to read platform name");

    err = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR,
                            sizeof(pInfo->vendor), pInfo->vendor, &readSize);
    if (readSize > sizeof(pInfo->vendor))
        mwPerrorCL(err, "Failed to read platform vendor");

    err = clGetPlatformInfo(platform, CL_PLATFORM_VERSION,
                            sizeof(pInfo->version), pInfo->version, &readSize);
    if (readSize > sizeof(pInfo->version))
        mwPerrorCL(err, "Failed to read platform version");

    err = clGetPlatformInfo(platform, CL_PLATFORM_EXTENSIONS,
                            sizeof(pInfo->extensions), pInfo->extensions, &readSize);
    if (readSize > sizeof(pInfo->extensions))
        mwPerrorCL(err, "Failed to read platform extensions");

    err = clGetPlatformInfo(platform, CL_PLATFORM_PROFILE,
                            sizeof(pInfo->profile), pInfo->profile, &readSize);
    if (readSize > sizeof(pInfo->profile))
        mwPerrorCL(err, "Failed to read platform profile");
}

static void mwPrintPlatformInfo(PlatformInfo* pInfo, cl_uint n)
{
    mw_printf("Platform %u information:\n"
              "  Name:       %s\n"
              "  Version:    %s\n"
              "  Vendor:     %s\n"
              "  Extensions: %s\n"
              "  Profile:    %s\n",
              n,
              pInfo->name,
              pInfo->version,
              pInfo->vendor,
              pInfo->extensions,
              pInfo->profile
        );
}

void mwPrintPlatforms(cl_platform_id* platforms, cl_uint n_platforms)
{
    cl_uint i;
    PlatformInfo pInfo = EMPTY_PLATFORM_INFO;

    for (i = 0; i < n_platforms; ++i)
    {
        mwGetPlatformInfo(&pInfo, platforms[i]);
        mwPrintPlatformInfo(&pInfo, i);
    }
}

cl_platform_id* mwGetAllPlatformIDs(cl_uint* nPlatformsOut)
{
    cl_uint nPlatform = 0;
    cl_uint nPlatformActual = 0;
    cl_platform_id* ids = NULL;
    cl_int err;

    err = clGetPlatformIDs(0, NULL, &nPlatform);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error getting number of platform");
        return NULL;
    }

    if (nPlatform == 0)
    {
        mw_printf("No CL platforms found\n");
        return NULL;
    }

    ids = mwMalloc(nPlatform * sizeof(cl_platform_id));
    err = clGetPlatformIDs(nPlatform, ids, &nPlatformActual);
    if ((err != CL_SUCCESS) || (nPlatformActual != nPlatform))
    {
        mwPerrorCL(err,
                   "Error getting platform IDs or inconsistent platform count (expected %u, actual %u)\n",
                   nPlatform,
                   nPlatformActual
            );
        free(ids);
        return NULL;
    }

    mw_printf("Found %u platform%s\n", nPlatform, nPlatform > 1 ? "s" : "");

    *nPlatformsOut = nPlatform;
    return ids;
}

cl_device_id* mwGetAllDevices(cl_platform_id platform, cl_uint* numDevOut)
{
    cl_int err;
    cl_device_id* devs;
    cl_uint numDev;

    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &numDev);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to find number of devices");
        return NULL;
    }

    if (numDev == 0)
    {
        mw_printf("Didn't find any CL devices\n");
        return NULL;
    }

    mw_printf("Found %u CL device%s\n", numDev, numDev > 1 ? "s" : "");

    devs = (cl_device_id*) mwMalloc(sizeof(cl_device_id) * numDev);
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, numDev, devs, &numDev);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to get device IDs");
        return NULL;
    }

    *numDevOut = numDev;
    return devs;
}

static cl_int mwGetDeviceType(cl_device_id dev, cl_device_type* devType)
{
    cl_int err = CL_SUCCESS;

    err = clGetDeviceInfo(dev, CL_DEVICE_TYPE, sizeof(cl_device_type), devType, NULL);
    if (err != CL_SUCCESS)
        mwPerrorCL(err, "Failed to get device type");

    return err;
}

cl_int mwSelectDevice(CLInfo* ci, const cl_device_id* devs, const CLRequest* clr, const cl_uint nDev)
{
    cl_int err = CL_SUCCESS;

    if (clr->devNum >= nDev)
    {
        mw_printf("Requested device is out of range of number found devices\n");
        return MW_CL_ERROR;
    }

    ci->dev = devs[clr->devNum];
    err = mwGetDeviceType(ci->dev, &ci->devType);
    if (err != CL_SUCCESS)
        mw_printf("Failed to find type of device %u\n", clr->devNum);

    return err;
}

cl_bool mwPlatformSupportsAMDOfflineDevices(const CLInfo* ci)
{
    cl_int err;
    char exts[4096];
    size_t readSize = 0;

    err = clGetPlatformInfo(ci->plat, CL_PLATFORM_EXTENSIONS, sizeof(exts), exts, &readSize);
    if ((err != CL_SUCCESS) || (readSize >= sizeof(exts)))
    {
        mwPerrorCL(err, "Error reading platform extensions (readSize = "ZU")\n", readSize);
        return CL_FALSE;
    }

    return (strstr(exts, "cl_amd_offline_devices") != NULL);
}

cl_bool mwNvidiaDriverVersionGreaterEqual(const DevInfo* di, cl_uint minMajor, cl_uint minMinor)
{
    cl_uint minor = 0;
    cl_uint major = 0;

    if (!mwIsNvidiaGPUDevice(di) || (sscanf(di->driver, "%u.%u", &major, &minor) != 2))
    {
        return CL_FALSE;
    }

    return (major > minMajor) || (major == minMajor && minor >= minMinor);
}

cl_bool mwNvidiaInlinePTXAvailable(cl_platform_id platform)
{
    cl_int err;
    size_t readSize = 0;
    char version[128];
    char name[128];
    cl_uint clMajor = 0, clMinor = 0;
    cl_uint cudaMajor = 0, cudaMinor = 0, cudaPatchLevel;

    err = clGetPlatformInfo(platform, CL_PLATFORM_NAME,
                            sizeof(name), name, &readSize);
    if (err != CL_SUCCESS || readSize >= sizeof(name))
    {
        return CL_FALSE;
    }

    err = clGetPlatformInfo(platform, CL_PLATFORM_VERSION,
                            sizeof(version), version, &readSize);
    if (err != CL_SUCCESS || readSize >= sizeof(version))
    {
        return CL_FALSE;
    }

    if (strcmp(name, "NVIDIA CUDA"))
    {
        return CL_FALSE;
    }

    /* Inline PTX was a CUDA 4 feature. Should probably test for this in a better way */
    if (  strcmp(name, "NVIDIA CUDA")
        || (sscanf(version,
                   "OpenCL %u.%u CUDA %u.%u.%u",
                   &clMinor, &clMajor,
                   &cudaMajor, &cudaMinor, &cudaPatchLevel) != 5))
    {
        return CL_FALSE;
    }

    return (cudaMajor >= 4);
}

