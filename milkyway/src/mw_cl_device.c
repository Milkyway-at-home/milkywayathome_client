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

#include <stdio.h>

#include "milkyway_cl.h"
#include "milkyway_util.h"
#include "mw_cl_show_types.h"
#include "mw_cl_setup.h"

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

    if (err)
        warn("Error getting device information: %s\n", showCLInt(err));
    else
        di->doubleExts = mwGetDoubleExts(di->exts);

    return err;
}

void mwPrintDevInfo(const DevInfo* di)
{
    warn("Device %s (%s:0x%x)\n"
         "Type:                %s\n"
         "Driver version:      %s\n"
         "Version:             %s\n"
         "Little endian:       %s\n"
         "Error correction:    %s\n"
         "Image support:       %s\n"
         "Address bits:        %u\n"
         "Max compute units:   %u\n"
         "Clock frequency:     %u Mhz\n"
         "Global mem size:     %llu\n"
         "Max mem alloc:       %llu\n"
         "Global mem cache:    %llu\n"
         "Cacheline size:      %u\n"
         "Local mem type:      %s\n"
         "Local mem size:      %llu\n"
         "Max const args:      %u\n"
         "Max const buf size:  %llu\n"
         "Max parameter size:  %zu\n"
         "Max work group size: %zu\n"
         "Max work item dim:   %u\n"
         "Max work item sizes: { %zu, %zu, %zu }\n"
         "Mem base addr align: %u\n"
         "Min type align size: %u\n"
         "Timer resolution:    %zu ns\n"
         "Double extension:    %s\n"
         "Extensions:          %s\n" ,
         di->devName,
         di->vendor,
         di->vendorID,
         showCLDeviceType(di->devType),
         di->driver,
         di->version,
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
         showMWDoubleExts(di->doubleExts),
         di->exts
        );
}

static void mwGetPlatformInfo(PlatformInfo* pi, cl_platform_id platform)
{
	cl_int err;
	size_t readSize;

	err = clGetPlatformInfo(platform, CL_PLATFORM_PROFILE,
                            sizeof(pi->name), pi->name, &readSize);
	if (readSize > sizeof(pi->name))
		warn("Failed to read platform profile: %s\n", showCLInt(err));


	err = clGetPlatformInfo(platform, CL_PLATFORM_VERSION,
                            sizeof(pi->version), pi->version, &readSize);
	if (readSize > sizeof(pi->version))
		warn("Failed to read platform version: %s\n", showCLInt(err));


	err = clGetPlatformInfo(platform, CL_PLATFORM_NAME,
                            sizeof(pi->name), pi->name, &readSize);
	if (readSize > sizeof(pi->name))
		warn("Failed to read platform name: %s\n", showCLInt(err));


	err = clGetPlatformInfo(platform, CL_PLATFORM_EXTENSIONS,
                            sizeof(pi->extensions), pi->extensions, &readSize);
	if (readSize > sizeof(pi->extensions))
		warn("Failed to read platform extensions: %s\n", showCLInt(err));
}

static void mwPrintPlatformInfo(PlatformInfo* pi, cl_uint n)
{
	warn("Platform %u information:\n"
	     "  Platform name:       %s\n"
	     "  Platform version:    %s\n"
	     "  Platform vendor:     %s\n"
	     "  Platform profile:    %s\n"
	     "  Platform extensions: %s\n",
	     n,
	     pi->name,
	     pi->version,
	     pi->vendor,
	     pi->profile,
	     pi->extensions);
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

cl_platform_id* mwGetAllPlatformIDs(CLInfo* ci, cl_uint* n_platforms_out)
{
    cl_uint n_platform = 0;
    cl_platform_id* ids;
    cl_int err;

    err = clGetPlatformIDs(0, NULL, &n_platform);
    if (err != CL_SUCCESS)
    {
        warn("Error getting number of platform: %s\n", showCLInt(err));
        return NULL;
    }

    ids = mallocSafe(sizeof(cl_platform_id) * n_platform);
    err = clGetPlatformIDs(n_platform, ids, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Error getting platform IDs: %s\n", showCLInt(err));
        free(ids);
        return NULL;
    }

    warn("Found %u platforms\n", n_platform);

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
        warn("Failed to find number of devices: %s\n", showCLInt(err));
        return NULL;
    }

    if (numDev == 0)
    {
        warn("Didn't find any CL devices\n");
        return NULL;
    }

    warn("Found %u CL devices\n", numDev);

    devs = (cl_device_id*) mallocSafe(sizeof(cl_device_id) * numDev);
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, numDev, devs, &numDev);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get device IDs: %s\n", showCLInt(err));
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
        warn("Failed to get device type: %s\n", showCLInt(err));

    return err;
}

cl_int mwSelectDevice(CLInfo* ci, const cl_device_id* devs, const CLRequest* clr, const cl_uint nDev)
{
    cl_int err = CL_SUCCESS;

    if (clr->devNum >= nDev)
    {
        warn("Requested device is out of range of number found devices\n");
        return -1;
    }

    ci->dev = devs[clr->devNum];
    err = mwGetDeviceType(ci->dev, &ci->devType);
    if (err != CL_SUCCESS)
        warn("Failed to find type of device %u\n", clr->devNum);

    return err;
}


