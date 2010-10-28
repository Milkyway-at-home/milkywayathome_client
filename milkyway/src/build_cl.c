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

cl_int mwSetOutOfOrder(CLInfo* ci)
{
    cl_int err;

    err = clSetCommandQueueProperty(ci->queue,
                                    CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE,
                                    CL_TRUE,
                                    NULL);
    if (err != CL_SUCCESS)
        warn("Setting out of order on queue failed: %s\n", showCLInt(err));

    return err;
}

cl_int mwEnableProfiling(CLInfo* ci)
{
    cl_int err;

    err = clSetCommandQueueProperty(ci->queue, CL_QUEUE_PROFILING_ENABLE, CL_TRUE, NULL);
    if (err != CL_SUCCESS)
        warn("Enabling profiling on queue failed: %s\n", showCLInt(err));

    return err;
}

cl_int mwDisableProfiling(CLInfo* ci)
{
    cl_int err;

    err = clSetCommandQueueProperty(ci->queue, CL_QUEUE_PROFILING_ENABLE, CL_FALSE, NULL);
    if (err != CL_SUCCESS)
        warn("Disabling profiling on queue failed: %s\n", showCLInt(err));

    return err;
}

/* Timing in nanoseconds */
cl_ulong mwEventTimeNS(cl_event ev)
{
    cl_int err;
    cl_ulong ts, te;

    err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &ts, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get event start time: %s\n", showCLInt(err));
        return 0;
    }

    err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &te, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get event end time: %s\n", showCLInt(err));
        return 0;
    }

    return te - ts;
}

/* Timing in seconds */
double mwEventTime(cl_event ev)
{
    return (double) mwEventTimeNS(ev) / 1.0e9;
}

/* Wait for an event then release it */
cl_int mwWaitReleaseEvent(cl_event* ev)
{
    cl_int err;

    assert(ev);

    err = clWaitForEvents(1, ev);
    if (err != CL_SUCCESS)
    {
        warn("%s: Failed to wait for event: %s\n", FUNC_NAME, showCLInt(err));
        return err;
    }

    err = clReleaseEvent(*ev);
    if (err != CL_SUCCESS)
    {
        warn("%s: Failed to release event: %s\n", FUNC_NAME, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

#ifdef CL_VERSION_1_1

cl_event mwCreateEvent(CLInfo* ci)
{
    cl_int err;
    cl_event ev;

    ev = clCreateUserEvent(ci->clctx, &err);

    if (err != CL_SUCCESS)
    {
        warn("Failed to create custom event: %s\n", showCLInt(err));
        return NULL;
    }

    return ev;
}

cl_int mwFinishEvent(cl_event ev)
{
    cl_int err;

    err = clSetUserEventStatus(ev, CL_COMPLETE);
    if (err != CL_SUCCESS)
        warn("Failed to mark custom event as completed: %s\n", showCLInt(err));

    return err;
}

#endif /* CL_VERSION_1_1 */

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
        warn("Extensions: %s\n", exts);

    return CL_SUCCESS;
}

cl_int getDevInfo(DevInfo* di, cl_device_id dev)
{
    cl_int err = CL_SUCCESS;

    err |= clGetDeviceInfo(dev, CL_DEVICE_TYPE,                     sizeof(di->devType),  &di->devType, NULL);

    err |= clGetDeviceInfo(dev, CL_DEVICE_NAME,                     sizeof(di->devName),  di->devName, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VENDOR,                   sizeof(di->vendor),   di->vendor, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VENDOR_ID,                sizeof(cl_uint),  &di->vendorID, NULL);
    err |= clGetDeviceInfo(dev, CL_DRIVER_VERSION,                  sizeof(di->driver),   di->driver, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_VERSION,                  sizeof(di->version),  di->version, NULL);
  //err |= clGetDeviceInfo(dev, CL_DEVICE_OPENCL_C_VERSION,         sizeof(di->clCVer),   di->clCVer, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_ENDIAN_LITTLE,            sizeof(cl_bool),  &di->littleEndian, NULL);
    err |= clGetDeviceInfo(dev, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(cl_bool),  &di->errCorrect, NULL);
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

    if (err)
        warn("Error getting device information: %s\n", showCLInt(err));
    return err;
}

void printDevInfo(const DevInfo* di)
{
    warn("Device %s (%s:0x%x)\n"
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
         "Timer resolution:    %zu ns\n" ,
         di->devName,
         di->vendor,
         di->vendorID,
         showCLDeviceType(di->devType),
         di->driver,
         di->version,
         showCLBool(di->littleEndian),
         showCLBool(di->errCorrect),
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
         di->timerRes
        );
}

cl_int destroyCLInfo(CLInfo* ci)
{
    cl_int err = CL_SUCCESS;
  /* Depending on where things fail, some of these will be NULL, and
   * will spew errors when trying to cleanup. */
    if (ci->queue)
        err |= clReleaseCommandQueue(ci->queue);
    if (ci->bufQueue)
        err |= clReleaseCommandQueue(ci->bufQueue);
    if (ci->prog)
        err |= clReleaseProgram(ci->prog);
    if (ci->kern)
        err |= clReleaseKernel(ci->kern);
    if (ci->clctx)
        err |= clReleaseContext(ci->clctx);

    /* TODO: or'ing the err and showing = useless */
    if (err)
        warn("Error cleaning up CLInfo: %s\n", showCLInt(err));

    return err;
}

cl_int getWorkGroupInfo(CLInfo* ci, WGInfo* wgi)
{
    cl_int err = CL_SUCCESS;

    err |= clGetKernelWorkGroupInfo(ci->kern,
                                    ci->dev,
                                    CL_KERNEL_COMPILE_WORK_GROUP_SIZE,
                                    sizeof(wgi->cwgs),
                                    wgi->cwgs,
                                    NULL);

    err |= clGetKernelWorkGroupInfo(ci->kern,
                                    ci->dev,
                                    CL_KERNEL_WORK_GROUP_SIZE,
                                    sizeof(size_t),
                                    &wgi->wgs,
                                    NULL);

    err |= clGetKernelWorkGroupInfo(ci->kern,
                                    ci->dev,
                                    CL_KERNEL_LOCAL_MEM_SIZE,
                                    sizeof(cl_ulong),
                                    &wgi->lms,
                                    NULL);
    if (err != CL_SUCCESS)
        warn("Failed to get kernel work group info: %s\n", showCLInt(err));

    return err;
}

void printWorkGroupInfo(const WGInfo* wgi)
{
    warn("Kernel work group info:\n"
         "  Work group size = %zu\n"
         "  Kernel local mem size = %llu\n"
         "  Compile work group size = { %zu, %zu, %zu }\n",
         wgi->wgs,
         wgi->lms,
         wgi->cwgs[0], wgi->cwgs[1], wgi->cwgs[2]);
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

/* Build program and create kernel */
static cl_int mwBuildProgram(CLInfo* ci, const char* options, const char* kernName)
{
    cl_int err = CL_SUCCESS;

    err = clBuildProgram(ci->prog, 1, &ci->dev, options, milkywayBuildCB, ci);
    if (err != CL_SUCCESS)
        warn("clBuildProgram: Build failure: %s\n", showCLInt(err));

    ci->kern = clCreateKernel(ci->prog, kernName, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating kernel '%s': %s\n", kernName, showCLInt(err));
        return err;
    }

    return err;
}

unsigned char* mwGetProgramBinary(CLInfo* ci, size_t* binSizeOut)
{
    cl_int err;
    size_t binSize;
    size_t readSize;
    unsigned char* bin = NULL;

    err = clGetProgramInfo(ci->prog, CL_PROGRAM_BINARY_SIZES, sizeof(binSize), &binSize, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get program binary size: %s\n", showCLInt(err));
        return NULL;
    }

    bin = (unsigned char*) mallocSafe(binSize);
    err = clGetProgramInfo(ci->prog, CL_PROGRAM_BINARIES, binSize, &bin, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Error getting program binary: %s\n", showCLInt(err));
        free(bin);
        bin = NULL;
        binSize = 0;
    }

    if (binSizeOut)
        *binSizeOut = binSize;
    return bin;
}

static cl_int createCtxQueue(CLInfo* ci, cl_bool useBufQueue)
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

    if (useBufQueue)
    {
        ci->bufQueue = clCreateCommandQueue(ci->clctx, ci->dev, 0, &err);
        if (err != CL_SUCCESS)
        {
            warn("Error creating buffer command queue: %s\n", showCLInt(err));
            return err;
        }
    }

    return CL_SUCCESS;
}

cl_int printDevInfoExts(const CLInfo* ci, const DevInfo* di)
{
    cl_int err;

    printDevInfo(di);

    err = printCLExtensions(ci->dev);
    if (err != CL_SUCCESS)
    {
        warn("Error getting printing device extensions: %s\n", showCLInt(err));
        return err;
    }
    return CL_SUCCESS;
}

typedef struct
{
	char name[128];
	char vendor[128];
	char version[128];
	char profile[128];
	char extensions[512];
} PlatformInfo;

#define EMPTY_PLATFORM_INFO { "", "", "", "", "" }

static void getPlatformInfo(PlatformInfo* pi, cl_platform_id platform)
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

static void printPlatformInfo(PlatformInfo* pi, cl_uint n)
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

static void printPlatforms(cl_platform_id* platforms, cl_uint n_platforms)
{
	cl_uint i;
	PlatformInfo pi = EMPTY_PLATFORM_INFO;

	for (i = 0; i < n_platforms; ++i)
	{
		getPlatformInfo(&pi, platforms[i]);
		printPlatformInfo(&pi, i);
	}
}

static cl_platform_id* getAllPlatformIDs(CLInfo* ci, cl_uint* n_platforms_out)
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

static cl_device_id* getAllDevices(cl_platform_id platform, cl_uint* numDevOut)
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

static cl_int getDeviceType(cl_device_id dev, cl_device_type* devType)
{
    cl_int err = CL_SUCCESS;

    err = clGetDeviceInfo(dev, CL_DEVICE_TYPE, sizeof(cl_device_type), devType, NULL);
    if (err != CL_SUCCESS)
        warn("Failed to get device type: %s\n", showCLInt(err));

    return err;
}

static cl_int selectDevice(CLInfo* ci, const cl_device_id* devs, const CLRequest* clr, const cl_uint nDev)
{
    cl_int err = CL_SUCCESS;

    if (clr->devNum >= nDev)
    {
        warn("Requested device is out of range of number found devices\n");
        return -1;
    }

    ci->dev = devs[clr->devNum];
    err = getDeviceType(ci->dev, &ci->devType);
    if (err != CL_SUCCESS)
        warn("Failed to find type of device %u\n", clr->devNum);

    return err;
}

static cl_int getCLInfo(CLInfo* ci, const CLRequest* clr)
{
    cl_int err = CL_SUCCESS;
    cl_uint n_platform;
    DevInfo di;
    cl_platform_id* ids;
    cl_device_id* devs;
    cl_uint nDev;

    ids = getAllPlatformIDs(ci, &n_platform);
    if (!ids)
    {
        warn("Failed to get any platforms\n");
        return -1;
    }

    printPlatforms(ids, n_platform);
    warn("Using device %u on platform %u\n", clr->devNum, clr->platform);

    devs = getAllDevices(ids[clr->platform], &nDev);
    if (!devs)
    {
        warn("Error getting devices\n");
        free(ids);
        return -1;
    }

    err = selectDevice(ci, devs, clr, nDev);
    if (err != CL_SUCCESS)
    {
        warn("Failed to select a device: %s\n", showCLInt(err));
        err = -1;
    }

    free(ids);
    free(devs);

    return err;
}

cl_int mwSetProgramFromBin(CLInfo* ci, const char* kernName, const unsigned char* bin, size_t binSize)
{
    cl_int err;
    cl_int binStatus;

    ci->prog = clCreateProgramWithBinary(ci->clctx, 1, &ci->dev, &binSize, &bin, &binStatus, &err);
    warn("Binary status: %s\n", showCLInt(err));
    if (err != CL_SUCCESS)
    {
        warn("Failed to create program from binary: %s\n", showCLInt(err));
        return err;
    }

    if (binStatus != CL_SUCCESS)
    {
        warn("Reading binary failed: %s\n", showCLInt(err));
        return binStatus;
    }

    err = mwBuildProgram(ci, NULL, kernName);
    if (err != CL_SUCCESS)
    {
        warn("Error building program from binary: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

cl_int mwSetProgramFromSrc(CLInfo* ci,
                           const char* kernName,
                           const char** src,
                           const cl_uint srcCount,
                           const char* compileDefs)
{
    cl_int err;

    ci->prog = clCreateProgramWithSource(ci->clctx, srcCount, src, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating program: %s\n", showCLInt(err));
        return err;
    }

    err = mwBuildProgram(ci, compileDefs, kernName);
    if (err != CL_SUCCESS)
    {
        warn("Error building program from source: %s\n", showCLInt(err));
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


/* Query one device specified by type, create a context and command
 * queue, as well as retrieve device information
 */
cl_int mwSetupCL(CLInfo* ci,
                 DevInfo* di,   /* get detailed device information for selected device */
                 const CLRequest* clr)
{
    cl_int err;

    err = getCLInfo(ci, clr);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get information about device\n");
        return err;
    }

    err = getDevInfo(di, ci->dev);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get device info\n");
        return err;
    }

    err = printDevInfoExts(ci, di);
    if (err != CL_SUCCESS)
    {
        warn("Failed to print device extensions\n");
        return err;
    }

    err = createCtxQueue(ci, CL_FALSE);
    if (err != CL_SUCCESS)
    {
        warn("Error creating CL context and command queue: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

