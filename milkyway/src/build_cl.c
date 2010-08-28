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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "milkyway_cl.h"
#include "milkyway_util.h"
#include "show_cl_types.h"
#include "build_cl.h"

/* This doesn't seem to exist on OS X, but the callback on ATI on
 * Linux/Windows dies without it */
#ifndef CL_CALLBACK
  #define CL_CALLBACK
#endif

cl_int printCLExtensions(cl_device_id dev)
{
    cl_int err;
    char exts[1024];

    err = clGetDeviceInfo(dev, CL_DEVICE_EXTENSIONS, sizeof(exts), exts, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to print CL extensions for device: %s", exts);
        return err;
    }
    else
        printf("Extensions: %s\n", exts);

    return CL_SUCCESS;
}

cl_int printDeviceInfo(cl_device_id dev, cl_device_type type)
{
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
    cl_int err = CL_SUCCESS;
    cl_uint cachelineSize;
    cl_bool littleEndian;
    cl_bool errCorrect;
    char devName[128];
    char vendor[128];
    char version[128];
    char driver[128];
    //char clCVer[128];

    size_t maxWorkItemSizes[3] = { 0, 0, 0 };

    err |= clGetDeviceInfo(dev, CL_DEVICE_NAME,                     sizeof(devName),  devName, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VENDOR,                   sizeof(vendor),   vendor, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VENDOR_ID,                sizeof(cl_uint),  &vendorID, NULL);
    err |= clGetDeviceInfo(dev, CL_DRIVER_VERSION,                  sizeof(driver),   driver, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VERSION,                  sizeof(version),  version, NULL);
  //err |= clGetDeviceInfo(dev, CL_DEVICE_OPENCL_C_VERSION,         sizeof(clCVer),   clCVer, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_ENDIAN_LITTLE,            sizeof(cl_bool),  &littleEndian, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(cl_bool),  &errCorrect, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_ADDRESS_BITS,             sizeof(cl_uint),  &addrBits, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS,        sizeof(cl_uint),  &maxCompUnits, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_CLOCK_FREQUENCY,      sizeof(cl_uint),  &clockFreq, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_SIZE,          sizeof(cl_ulong), &memSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_MEM_ALLOC_SIZE,       sizeof(cl_ulong), &maxMemAlloc, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,    sizeof(cl_ulong), &gMemCache, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, sizeof(cl_uint), &cachelineSize, NULL);

  //err |= clGetDeviceInfo(dev, CL_DEVICE_HOST_UNIFIED_MEMORY,      sizeof(cl_ulong), &unifiedMem, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(cl_device_local_mem_type), &localMemType, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_LOCAL_MEM_SIZE,           sizeof(cl_ulong), &localMemSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_CONSTANT_ARGS,        sizeof(cl_uint),  &maxConstArgs, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong), &maxConstBufSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_PARAMETER_SIZE, sizeof(size_t), &maxParamSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &maxWorkGroupSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &maxWorkItemDim, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(maxWorkItemSizes), maxWorkItemSizes, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MEM_BASE_ADDR_ALIGN, sizeof(cl_uint), &memBaseAddrAlign, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, sizeof(cl_uint), &minAlignSize, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_PROFILING_TIMER_RESOLUTION, sizeof(size_t), &timerRes, NULL);




    printf("Device %s (%s:0x%x)\n"
           "Type:                %s\n"
           "Driver version:      %s\n"
           "Version:             %s\n"
           "Little endian:       %s\n"
           "Error correction:    %s\n"
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
           ,
           devName,
           vendor,
           vendorID,
           showCLDeviceType(type),
           driver,
           version,
           showCLBool(littleEndian),
           showCLBool(errCorrect),
           addrBits,
           maxCompUnits,
           clockFreq,
           memSize,
           maxMemAlloc,
           gMemCache,
           cachelineSize,
           //showBool(unifiedMem),
           showCLDeviceLocalMemType(localMemType),
           localMemSize,
           maxConstArgs,
           maxConstBufSize,
           maxParamSize,
           maxWorkGroupSize,
           maxWorkItemDim,
           maxWorkItemSizes[0], maxWorkItemSizes[1], maxWorkItemSizes[2],
           memBaseAddrAlign,
           minAlignSize,
           timerRes
        );

    if (err)
        warn("Error getting device information: %s\n", showCLInt(err));
    return err;
}

cl_int destroyCLInfo(CLInfo* ci)
{
    cl_int err = CL_SUCCESS;
    err |= clReleaseCommandQueue(ci->queue);
    err |= clReleaseProgram(ci->prog);
    err |= clReleaseKernel(ci->kern);
    err |= clReleaseContext(ci->clctx);

    /* TODO: or'ing the err and showing = useless */
    if (err)
        warn("Error cleaning up CLInfo: %s\n", showCLInt(err));

    return err;
}

static char* getBuildLog(CLInfo* ci)
{
    size_t logSize, readSize;
    char* buildLog;

    clGetProgramBuildInfo(ci->prog,
                          ci->dev,
                          CL_PROGRAM_BUILD_LOG,
                          0,
                          NULL,
                          &logSize);

    buildLog = callocSafe(sizeof(char), logSize + 1);

    clGetProgramBuildInfo(ci->prog,
                          ci->dev,
                          CL_PROGRAM_BUILD_LOG,
                          logSize,
                          buildLog,
                          &readSize);

    if (readSize != logSize)
        warn("Failed to read complete build log\n");

    return buildLog;
}

static void CL_CALLBACK milkywayBuildCB(cl_program prog, void* user_data)
{
    cl_int infoErr;
    cl_build_status stat;
    CLInfo* ci;
    char* buildLog;

    if (!user_data)
    {
        warn("milkywayBuildCB got null user_data\n");
        return;
    }

    ci = (CLInfo*) user_data;

    infoErr = clGetProgramBuildInfo(ci->prog,
                                    ci->dev,
                                    CL_PROGRAM_BUILD_STATUS,
                                    sizeof(stat),
                                    &stat,
                                    NULL);

    if (infoErr != CL_SUCCESS)
        warn("Get build status failed: %s\n", showCLInt(infoErr));
    else
        warn("Build status: %s\n", showCLBuildStatus(stat));

    buildLog = getBuildLog(ci);

    warn("Build log: \n%s\n", buildLog);
    free(buildLog);
}

static cl_int milkywayBuildProgram(CLInfo* ci, const char** src, cl_uint srcCount, const char* compileDefs)
{
    cl_int err = CL_SUCCESS;

    ci->prog = clCreateProgramWithSource(ci->clctx, srcCount, src, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating program: %s\n", showCLInt(err));
        return err;
    }

    err = clBuildProgram(ci->prog, 1, &ci->dev, compileDefs, milkywayBuildCB, ci);
    if (err != CL_SUCCESS)
        warn("clBuildProgram: Build failure: %s\n", showCLInt(err));

    return err;
}

static cl_int createCtxQueue(CLInfo* ci)
{
    cl_int err = CL_SUCCESS;

    ci->clctx = clCreateContext(NULL, 1, &ci->dev, NULL, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating context: %s\n", showCLInt(err));
        return err;
    }

    ci->queue = clCreateCommandQueue(ci->clctx, ci->dev, 0, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating command queue: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static cl_int getDevInfo(CLInfo* ci, cl_device_type type)
{
    cl_int err;
    cl_platform_id platform;
    cl_uint n_platform;

    err = clGetPlatformIDs(1, &platform, &n_platform);
    if (err != CL_SUCCESS)
    {
        warn("Error getting CL platform IDs: %s\n", showCLInt(err));
        return err;
    }

    warn("Found %u platforms\n", n_platform);

    err = clGetDeviceIDs(platform, type, 1, &ci->dev, &ci->devCount);
    if (err != CL_SUCCESS)
    {
        warn("Error getting device: %s\n", showCLInt(err));
        return err;
    }

    if (ci->devCount == 0)
    {
        warn("Didn't find any %s devices\n", showCLDeviceType(type));
        return -1; /* FIXME: Meaningful error? */
    }

    ci->devType = type;

    err = printDeviceInfo(ci->dev, type);
    if (err != CL_SUCCESS)
    {
        warn("Error getting printing device information: %s\n", showCLInt(err));
        return err;
    }

    err = printCLExtensions(ci->dev);
    if (err != CL_SUCCESS)
    {
        warn("Error getting printing device extensions: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

/* Query one device specified by type, create a context, command
 * queue, and compile the kernel for it
 * Returns: non-zero on failure
 *  TODO: Multiple device support
 *  TODO: Caching of compiled binaries
 */
cl_int getCLInfo(CLInfo* ci,
                 cl_device_type type,
                 const char* kernName,
                 const char** src,
                 const cl_uint srcCount,
                 const char* compileDefs)
{
    cl_int err;

    err = getDevInfo(ci, type);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get information about device\n");
        return err;
    }

    err = createCtxQueue(ci);
    if (err != CL_SUCCESS)
    {
        warn("Error creating CL context and command queue: %s\n", showCLInt(err));
        return err;
    }

    err = milkywayBuildProgram(ci, src, srcCount, compileDefs);
    if (err != CL_SUCCESS)
    {
        warn("Error building program: %s\n", showCLInt(err));
        return err;
    }

    ci->kern = clCreateKernel(ci->prog, kernName, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating kernel '%s': %s\n", kernName, showCLInt(err));
        return err;
    }

    clUnloadCompiler();

    return CL_SUCCESS;
}

