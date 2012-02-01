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

#include "milkyway_util.h"
#include "milkyway_cl_show_types.h"

const char* showCLDeviceType(const cl_device_type x)
{
    switch (x)
    {
        case CL_DEVICE_TYPE_CPU:
            return "CL_DEVICE_TYPE_CPU";
        case CL_DEVICE_TYPE_GPU:
            return "CL_DEVICE_TYPE_GPU";
        case CL_DEVICE_TYPE_ACCELERATOR:
            return "CL_DEVICE_TYPE_ACCELERATOR";
        case CL_DEVICE_TYPE_DEFAULT:
            return "CL_DEVICE_TYPE_DEFAULT";
        case CL_DEVICE_TYPE_ALL:
            return "CL_DEVICE_TYPE_ALL";
        default:
            mw_printf("Trying to show unknown cl_device_type %d\n", (int) x);
            return "Unhandled cl_device_type";
    }
}

const char* showCLBuildStatus(const cl_build_status x)
{
    switch (x)
    {
        case CL_BUILD_NONE:
            return "CL_BUILD_NONE";
        case CL_BUILD_ERROR:
            return "CL_BUILD_ERROR";
        case CL_BUILD_SUCCESS:
            return "CL_BUILD_SUCCESS";
        case CL_BUILD_IN_PROGRESS:
            return "CL_BUILD_IN_PROGRESS";
        default:
            mw_printf("Trying to show unknown cl_build_status %d\n", x);
            return "Unknown cl_build_status";
    }
}

/* ErrorCode -> String for easier debugging  */
const char* showCLInt(const cl_int x)
{
    switch (x)
    {
        case MW_CL_ERROR:  /* Custom error */
            return "MW_CL_ERROR";

        case CL_SUCCESS:
            return "CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND:
            return "CL_DEVICE_NOT_FOUND";
        case CL_DEVICE_NOT_AVAILABLE:
            return "CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE:
            return "CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES:
            return "CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY:
            return "CL_OUT_OF_HOST_MEMORY";
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP:
            return "CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH:
            return "CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE:
            return "CL_BUILD_PROGRAM_FAILURE";
        case CL_MAP_FAILURE:
            return "CL_MAP_FAILURE";

        case CL_INVALID_VALUE:
            return "CL_INVALID_VALUE";
        case CL_INVALID_DEVICE_TYPE:
            return "CL_INVALID_DEVICE_TYPE";
        case CL_INVALID_PLATFORM:
            return "CL_INVALID_PLATFORM";
        case CL_INVALID_DEVICE:
            return "CL_INVALID_DEVICE";
        case CL_INVALID_CONTEXT:
            return "CL_INVALID_CONTEXT";
        case CL_INVALID_QUEUE_PROPERTIES:
            return "CL_INVALID_QUEUE_PROPERTIES";
        case CL_INVALID_COMMAND_QUEUE:
            return "CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_HOST_PTR:
            return "CL_INVALID_HOST_PTR";
        case CL_INVALID_MEM_OBJECT:
            return "CL_INVALID_MEM_OBJECT";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case CL_INVALID_IMAGE_SIZE:
            return "CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_SAMPLER:
            return "CL_INVALID_SAMPLER";
        case CL_INVALID_BINARY:
            return "CL_INVALID_BINARY";
        case CL_INVALID_BUILD_OPTIONS:
            return "CL_INVALID_BUILD_OPTIONS";
        case CL_INVALID_PROGRAM:
            return "CL_INVALID_PROGRAM";
        case CL_INVALID_PROGRAM_EXECUTABLE:
            return "CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_KERNEL_NAME:
            return "CL_INVALID_KERNEL_NAME";
        case CL_INVALID_KERNEL_DEFINITION:
            return "CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL:
            return "CL_INVALID_KERNEL";
        case CL_INVALID_ARG_INDEX:
            return "CL_INVALID_ARG_INDEX";
        case CL_INVALID_ARG_VALUE:
            return "CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE:
            return "CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS:
            return "CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_WORK_DIMENSION:
            return "CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE:
            return "CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE:
            return "CL_INVALID_WORK_ITEM_SIZE";
        case CL_INVALID_GLOBAL_OFFSET:
            return "CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_EVENT_WAIT_LIST:
            return "CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_EVENT:
            return "CL_INVALID_EVENT";
        case CL_INVALID_OPERATION:
            return "CL_INVALID_OPERATION";
        case CL_INVALID_GL_OBJECT:
            return "CL_INVALID_GL_OBJECT";
        case CL_INVALID_BUFFER_SIZE:
            return "CL_INVALID_BUFFER_SIZE";
        case CL_INVALID_MIP_LEVEL:
            return "CL_INVALID_MIP_LEVEL";
        case CL_INVALID_GLOBAL_WORK_SIZE:
            return "CL_INVALID_GLOBAL_WORK_SIZE";

     #ifdef CL_VERSION_1_1
        case CL_PLATFORM_NOT_FOUND_KHR:
            return "CL_PLATFORM_NOT_FOUND_KHR";
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
            return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        case CL_MISALIGNED_SUB_BUFFER_OFFSET:
            return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
            /*
        case CL_INVALID_D3D10_RESOURCE_KHR:
            return "CL_INVALID_D3D10_RESOURCE_KHR";
            */
     #endif /* CL_VERSION_1_1 */

        default:
            mw_printf("Trying to show unknown cl_int %d\n", x);
            return "Unknown cl_int";
    }
}

const char* showCLMemFlags(const cl_mem_flags x)
{
    switch (x)
    {
        case CL_MEM_READ_WRITE:
            return "CL_MEM_READ_WRITE";
        case CL_MEM_WRITE_ONLY:
            return "CL_MEM_WRITE_ONLY";
        case CL_MEM_READ_ONLY:
            return "CL_MEM_READ_ONLY";
        case CL_MEM_USE_HOST_PTR:
            return "CL_MEM_USE_HOST_PTR";
        case CL_MEM_ALLOC_HOST_PTR:
            return "CL_ALLOC_USE_HOST_PTR";
        default:
            mw_printf("Trying to show unknown cl_mem_flags %d\n", (int) x);
            return "Unhandled cl_mem_flags";
    }
}

const char* showCLDeviceFPConfig(const cl_device_fp_config x)
{
    switch (x)
    {
        case CL_FP_DENORM:
            return "CL_FP_DENORM";
        case CL_FP_INF_NAN:
            return "CL_FP_INF_NAN";
        case CL_FP_ROUND_TO_NEAREST:
            return "CL_FP_ROUND_TO_NEAREST";
        case CL_FP_ROUND_TO_ZERO:
            return "CL_FP_ROUND_TO_ZERO";
        case CL_FP_ROUND_TO_INF:
            return "CL_FP_ROUND_TO_INF";
            /* FIXME:: Typeo in opencl docs? */
        case CL_FP_FMA:
            return "CL_FP_FMA";

            /*
        case CL_FP_SOFT_FLOAT:
            return "CL_FP_SOFT_FLOAT";
            */
        default:
            mw_printf("Trying to show unknown cl_device_local_mem_type %d\n", (int) x);
            return "Unhandled cl_local_mem_type";
    }
}

const char* showCLDeviceLocalMemType(const cl_device_local_mem_type x)
{
    switch (x)
    {
        case CL_LOCAL:
            return "CL_LOCAL";
        case CL_GLOBAL:
            return "CL_GLOBAL";
        default:
            mw_printf("Trying to show unknown cl_device_local_mem_type %d\n", (int) x);
            return "Unhandled cl_local_mem_type";
    }
}

const char* showCLDeviceExecCapabilities(const cl_device_exec_capabilities x)
{
    switch (x)
    {
        case CL_EXEC_KERNEL:
            return "CL_EXEC_KERNEL";
        case CL_EXEC_NATIVE_KERNEL:
            return "CL_EXEC_NATIVE_KERNEL";
        default:
            mw_printf("Trying to show unknown cl_device_exec_capabilities %d\n", (int) x);
            return "Unhandled cl_device_exec_capabilities";
    }
}

const char* showCLCommandQueueProperties(const cl_command_queue_properties x)
{
    switch (x)
    {
        case CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE:
            return "CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE";
        case CL_QUEUE_PROFILING_ENABLE:
            return "CL_QUEUE_PROFILING_ENABLE";
        default:
            mw_printf("Trying to show unknown cl_command_queue_properties %d\n", (int) x);
            return "Unhandled cl_command_queue_properties";
    }
}

const char* showCLBool(const cl_bool x)
{
    switch (x)
    {
        case CL_TRUE:
            return "CL_TRUE";
        case CL_FALSE:
            return "CL_FALSE";
        default:
            mw_printf("Trying to show unknown cl_bool %d\n", (int) x);
            return "Unhandled cl_bool";
    }
}

const char* showCLDeviceMemCacheType(const cl_device_mem_cache_type x)
{
    switch (x)
    {
        case CL_NONE:
            return "CL_NONE";
        case CL_READ_ONLY_CACHE:
            return "CL_READ_ONLY_CACHE";
        case CL_READ_WRITE_CACHE:
            return "CL_READ_WRITE_CACHE";
        default:
            mw_printf("Trying to show unknown cl_device_mem_cache_type %d\n", (int) x);
            return "Unhandled cl_device_mem_cache_type";
    }
}

const char* showCLKernelInfo(const cl_kernel_info x)
{
    switch (x)
    {
        case CL_KERNEL_FUNCTION_NAME:
            return "CL_KERNEL_FUNCTION_NAME";
        case CL_KERNEL_NUM_ARGS:
            return "CL_KERNEL_NUM_ARGS";
        case CL_KERNEL_REFERENCE_COUNT:
            return "CL_KERNEL_REFERENCE_COUNT";
        case CL_KERNEL_CONTEXT:
            return "CL_KERNEL_CONTEXT";
        case CL_KERNEL_PROGRAM:
            return "CL_KERNEL_PROGRAM";
        default:
            mw_printf("Trying to show unknown cl_kernel_info %d\n", (int) x);
            return "Unhandled cl_kernel_info";
    }
}

const char* showMWDoubleExts(const MWDoubleExts x)
{

    if (x & MW_CL_KHR_FP64)  /* Preferred extension */
    {
        return "cl_khr_fp64";
    }
    else if (x & MW_CL_AMD_FP64)
    {
        return "cl_amd_fp64";
    }
    else
    {
        return "(none)";
    }

    mw_unreachable();
    return NULL;
}

const char* showMWCALtargetEnum(const MWCALtargetEnum x)
{
    switch (x)
    {
        case MW_CAL_TARGET_UNKNOWN:
            return "MW_CAL_TARGET_UNKNOWN";
        case MW_CAL_TARGET_INVALID:
            return "MW_CAL_TARGET_INVALID";
        case MW_CAL_TARGET_600:
            return "MW_CAL_TARGET_600";
        case MW_CAL_TARGET_610:
            return "MW_CAL_TARGET_610";
        case MW_CAL_TARGET_630:
            return "MW_CAL_TARGET_630";
        case MW_CAL_TARGET_670:
            return "MW_CAL_TARGET_670";
        case MW_CAL_TARGET_7XX:
            return "MW_CAL_TARGET_7XX";
        case MW_CAL_TARGET_770:
            return "MW_CAL_TARGET_770";
        case MW_CAL_TARGET_710:
            return "MW_CAL_TARGET_710";
        case MW_CAL_TARGET_730:
            return "MW_CAL_TARGET_730";
        case MW_CAL_TARGET_CYPRESS:
            return "MW_CAL_TARGET_CYPRESS";
        case MW_CAL_TARGET_JUNIPER:
            return "MW_CAL_TARGET_JUNIPER";
        case MW_CAL_TARGET_REDWOOD:
            return "MW_CAL_TARGET_REDWOOD";
        case MW_CAL_TARGET_CEDAR:
            return "MW_CAL_TARGET_CEDAR";

        case MW_CAL_TARGET_SUMO:
            return "MW_CAL_TARGET_SUMO";
        case MW_CAL_TARGET_SUPERSUMO:
            return "MW_CAL_TARGET_SUPERSUMO";

        case MW_CAL_TARGET_WRESTLER:
            return "MW_CAL_TARGET_WRESTLER";
        case MW_CAL_TARGET_CAYMAN:
            return "MW_CAL_TARGET_CAYMAN";
        case MW_CAL_TARGET_RESERVED2:
            return "MW_CAL_TARGET_RESERVED2";
        case MW_CAL_TARGET_BARTS:
            return "MW_CAL_TARGET_BARTS";

        case MW_CAL_TARGET_TURKS:
            return "MW_CAL_TARGET_TURKS";
        case MW_CAL_TARGET_CAICOS:
            return "MW_CAL_TARGET_CAICOS";

        case MW_CAL_TARGET_TAHITI:
            return "MW_CAL_TARGET_TAHITI";
        case MW_CAL_TARGET_THAMES:
            return "MW_CAL_TARGET_THAMES";
        case MW_CAL_TARGET_LOMBOK:
            return "MW_CAL_TARGET_LOMBOK";

        default:
            return "Invalid MWCALtargetEnum";
    }
}

