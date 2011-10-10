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

#include "mw_cl_device.h"
#include "mw_cl_util.h"
#include "milkyway_util.h"
#include "mw_cl_show_types.h"

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
    cl_uint vliw;
    cl_uint doubleFrac;
    cl_uint wavefrontSize;
} AMDGPUData;

/* A reasonable guess */
static const AMDGPUData invalidAMDGPUData = { "Invalid", MW_CAL_TARGET_INVALID, 5, 5, 64 };

/* Table of more detailed info we want that OpenCL won't give us */
static const AMDGPUData amdGPUData[] =
{
    /* I found these names from strings on amdocl.dll/libamdocl64.so
       There are a bunch more extra ones in fglrx.ko (e.g. "Hemlock"),
       but I don't think they're relevant.
     */

    /* R600 */
    { "ATI RV670",  MW_CAL_TARGET_670,      5, 0, 64 },
    { "ATI RV630",  MW_CAL_TARGET_630,      5, 0, 32 },
    { "ATI RV610",  MW_CAL_TARGET_610,      5, 0, 32 },
    { "ATI RV600",  MW_CAL_TARGET_600,      5, 0, 16 },

    /* R700 */
    { "ATI RV770",  MW_CAL_TARGET_770,      5, 5, 64 },
    { "ATI RV730",  MW_CAL_TARGET_730,      5, 0, 32 },
    { "ATI RV710",  MW_CAL_TARGET_710,      5, 0, 16 },

    /* Evergreen */
    { "Cypress",    MW_CAL_TARGET_CYPRESS,  5, 5, 64 },
    { "Juniper",    MW_CAL_TARGET_JUNIPER,  5, 0, 64 },
    { "Redwood",    MW_CAL_TARGET_REDWOOD,  5, 0, 64 },
    { "Cedar",      MW_CAL_TARGET_CEDAR,    5, 0, 32 },

    /* Northern Islands */
    { "Cayman",     MW_CAL_TARGET_CAYMAN,   4, 4, 64 },
    { "Barts",      MW_CAL_TARGET_BARTS,    5, 0, 64 },
    { "Turks",      MW_CAL_TARGET_TURKS,    5, 0, 64 },
    { "Caicos",     MW_CAL_TARGET_CAICOS,   5, 0, 64 },


    /* Southern Islands */
    /*
      Not actually out yet, so not really sure
    { "Tahiti",     MW_CAL_TARGET_TAHITI,   4, 4, 64 },
    { "Thames",     MW_CAL_TARGET_THAMES,   5, 0, 64 },
    { "Lombok",     MW_CAL_TARGET_LOMBOK,   5, 0, 64 },
    */

#if 0
    /* These are there, but I don't know about them */
    { "WinterPark",                        , 5, 0, 64 },
    { "BeaverCreek",                       , 5, 0, 64 },
    { "Loveland",                          , 5, 0, 64 },

    { "Lions",                             , 5, 0, 64 },
    { "Tigers",                            , 5, 0, 64 },
    { "Bears",                             , 5, 0, 64 },
#endif

    { NULL,          MW_CAL_TARGET_INVALID, 5, 0, 64 }
};

cl_uint uavIdFromMWCALtargetEnum(MWCALtargetEnum x)
{
    return (x >= MW_CAL_TARGET_CYPRESS) ? 11 : 1;
}

static const AMDGPUData* lookupAMDGPUInfo(const DevInfo* di)
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

cl_double amdEstimateGFLOPs(const DevInfo* di, cl_bool useDouble)
{
    cl_uint vliw = di->vliw;
    cl_uint doubleFrac = di->doubleFrac;
    cl_ulong flops, flopsFloat, flopsDouble;
    cl_double gflops;

    flopsFloat = 2 * (di->maxCompUnits * vliw * 16) * (cl_ulong) di->clockFreq * 1000000;
    flopsDouble = flopsFloat / doubleFrac;

    mw_printf("Estimated AMD GPU (VLIW%d) GFLOP/s: %.0f SP GFLOP/s, %.0f DP FLOP/s\n",
              vliw,
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

cl_bool hasNvidiaCompilerFlags(const DevInfo* di)
{
    return strstr(di->exts, "cl_nv_compiler_options") != NULL;
}

cl_bool deviceVendorIsAMD(const DevInfo* di)
{
    return (strncmp(di->vendor, "Advanced Micro Devices, Inc.", sizeof(di->vendor)) == 0);
}

cl_bool deviceVendorIsNvidia(const DevInfo* di)
{
    return (strncmp(di->vendor, "Nvidia Corporation", sizeof(di->vendor)) == 0);
}

/* True if devices compute capability >= requested version */
cl_bool minComputeCapabilityCheck(const DevInfo* di, cl_uint major, cl_uint minor)
{
    return     di->computeCapabilityMajor > major
           || (   di->computeCapabilityMajor == major
               && di->computeCapabilityMinor >= minor);
}

/* Exact check on compute capability version */
cl_bool computeCapabilityIs(const DevInfo* di, cl_uint major, cl_uint minor)
{
    return di->computeCapabilityMajor == major && di->computeCapabilityMinor == minor;
}

/* approximate ratio of float : double flops */
static cl_uint cudaEstimateDoubleFrac(const DevInfo* di)
{
    /* FIXME: This also differs with generation.
       Is there a better way to find out the generation and if
     */

    if (minComputeCapabilityCheck(di, 2, 0))
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
static cl_uint cudaCoresPerComputeUnit(const DevInfo* di)
{
    if (minComputeCapabilityCheck(di, 2, 0))
        return 32;

    return 8;     /* 1.x is 8 */
}

cl_double cudaEstimateGFLOPs(const DevInfo* di, cl_bool useDouble)
{
    cl_ulong flopsFloat, flopsDouble, flops;
    cl_uint doubleRat = cudaEstimateDoubleFrac(di);
    cl_uint corePerCU = cudaCoresPerComputeUnit(di);
    cl_uint numCUDACores = di->maxCompUnits * corePerCU;
    cl_double gflops;

    flopsFloat = 2 * numCUDACores * (cl_ulong) di->clockFreq * 1000000;
    flopsDouble = flopsFloat / doubleRat;

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

cl_bool isNvidiaGPUDevice(const DevInfo* di)
{
    return (di->vendorID == MW_NVIDIA) && (di->devType == CL_DEVICE_TYPE_GPU);
}

cl_bool isAMDGPUDevice(const DevInfo* di)
{
    /* Not sure if the vendor ID for AMD is the same with their
       CPUs.  Also something else weird was going on with the
       vendor ID, so check the name just in case.
    */

    return (di->vendorID == MW_AMD_ATI || deviceVendorIsAMD(di));
}

cl_double deviceEstimateGFLOPs(const DevInfo* di, cl_bool useDouble)
{
    cl_double gflops = 0.0;

    if (di->devType == CL_DEVICE_TYPE_GPU)
    {
        if (isNvidiaGPUDevice(di))
        {
            gflops = cudaEstimateGFLOPs(di, useDouble);
        }
        else if (isAMDGPUDevice(di))
        {
            gflops = amdEstimateGFLOPs(di, useDouble);
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

    /* TODO: Correct way to test for Tesla type things?
       Is there a way we can find if a GPU is connected to a display in general?
     */
    di->nonOutput = ((di->devType != CL_DEVICE_TYPE_GPU) || (strstr(di->devName, "Tesla") != NULL));


    if (isNvidiaGPUDevice(di))
    {
        di->vliw = 1;
        di->doubleFrac = cudaEstimateDoubleFrac(di);
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
    else if (isAMDGPUDevice(di))
    {
        amdData = lookupAMDGPUInfo(di);

        di->vliw       = amdData->vliw;
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
        mwPerrorCL("Error getting device information", err);
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
              "VLIW:                %u\n"
              "Double extension:    %s\n"
              "Double fraction:     1/%u\n"
              "Extensions:          %s\n",
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
              di->vliw,
              showMWDoubleExts(di->doubleExts),
              di->doubleFrac,
              di->exts
        );
}

void mwPrintDevInfoShort(const DevInfo* di)
{
    mw_printf("Device %s (%s:0x%x) (%s)\n"
              "Driver version:      %s\n"
              "Version:             %s\n"
              "Compute capability:  %u.%u\n"
              "Image support:       %s\n"
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
              showCLBool(di->imgSupport),
              di->maxCompUnits,
              di->clockFreq,
              di->memSize,
              di->localMemSize,
              di->maxConstBufSize,
              showMWDoubleExts(di->doubleExts)
        );
}

static void mwGetPlatformInfo(PlatformInfo* pi, cl_platform_id platform)
{
    cl_int err;
    size_t readSize;

    err = clGetPlatformInfo(platform, CL_PLATFORM_NAME,
                            sizeof(pi->name), pi->name, &readSize);
    if (readSize > sizeof(pi->name))
        mwPerrorCL("Failed to read platform name", err);

    err = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR,
                            sizeof(pi->vendor), pi->vendor, &readSize);
    if (readSize > sizeof(pi->vendor))
        mwPerrorCL("Failed to read platform vendor", err);

    err = clGetPlatformInfo(platform, CL_PLATFORM_VERSION,
                            sizeof(pi->version), pi->version, &readSize);
    if (readSize > sizeof(pi->version))
        mwPerrorCL("Failed to read platform version", err);

    err = clGetPlatformInfo(platform, CL_PLATFORM_EXTENSIONS,
                            sizeof(pi->extensions), pi->extensions, &readSize);
    if (readSize > sizeof(pi->extensions))
        mwPerrorCL("Failed to read platform extensions", err);

    err = clGetPlatformInfo(platform, CL_PLATFORM_PROFILE,
                            sizeof(pi->profile), pi->profile, &readSize);
    if (readSize > sizeof(pi->profile))
        mwPerrorCL("Failed to read platform profile", err);
}

static void mwPrintPlatformInfo(PlatformInfo* pi, cl_uint n)
{
    mw_printf("Platform %u information:\n"
              "  Name:       %s\n"
              "  Version:    %s\n"
              "  Vendor:     %s\n"
              "  Extensions: %s\n"
              "  Profile:    %s\n",
              n,
              pi->name,
              pi->version,
              pi->vendor,
              pi->extensions,
              pi->profile
        );
}

void mwPrintPlatforms(cl_platform_id* platforms, cl_uint n_platforms)
{
    cl_uint i;
    PlatformInfo pi = EMPTY_PLATFORM_INFO;

    for (i = 0; i < n_platforms; ++i)
    {
        mwGetPlatformInfo(&pi, platforms[i]);
        mwPrintPlatformInfo(&pi, i);
    }
}

cl_platform_id* mwGetAllPlatformIDs(cl_uint* n_platforms_out)
{
    cl_uint n_platform = 0;
    cl_platform_id* ids;
    cl_int err;

    err = clGetPlatformIDs(0, NULL, &n_platform);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL("Error getting number of platform", err);
        return NULL;
    }

    if (n_platform == 0)
    {
        mw_printf("No CL platforms found\n");
        return NULL;
    }

    ids = mwMalloc(sizeof(cl_platform_id) * n_platform);
    err = clGetPlatformIDs(n_platform, ids, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL("Error getting platform IDs", err);
        free(ids);
        return NULL;
    }

    mw_printf("Found %u platform(s)\n", n_platform);

    *n_platforms_out = n_platform;
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
        mwPerrorCL("Failed to find number of devices", err);
        return NULL;
    }

    if (numDev == 0)
    {
        mw_printf("Didn't find any CL devices\n");
        return NULL;
    }

    mw_printf("Found %u CL device(s)\n", numDev);

    devs = (cl_device_id*) mwMalloc(sizeof(cl_device_id) * numDev);
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, numDev, devs, &numDev);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL("Failed to get device IDs", err);
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
        mwPerrorCL("Failed to get device type", err);

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
        mwPerrorCL("Error reading platform extensions (readSize = "ZU")\n", err, readSize);
        return CL_FALSE;
    }

    return (strstr(exts, "cl_amd_offline_devices") != NULL);
}

