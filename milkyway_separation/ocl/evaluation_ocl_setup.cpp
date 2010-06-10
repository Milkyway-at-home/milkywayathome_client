/*
Copyright 2010 Anthony Waters, Travis Desell,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#include "evaluation_ocl.h"
#include "evaluation_ocl_priv.h"

char* get_platform_info(cl_platform_id platform,
                        cl_platform_info info)
{
    size_t size;
    check_error(clGetPlatformInfo(platform,
                                  info,
                                  0, 0, &size));
    char* profile = new char[size];
    check_error(clGetPlatformInfo(platform,
                                  info,
                                  size, profile, 0));
    return profile;
}

cl_platform_id* get_platforms()
{
    cl_uint num_platforms;
    check_error(clGetPlatformIDs(0, 0, &num_platforms));
    printf("Found %d platforms\n", num_platforms);
    if (num_platforms > 0)
    {
        cl_platform_id* platforms = new cl_platform_id[num_platforms];
        check_error(clGetPlatformIDs(num_platforms,
                                     platforms, 0));
        for (unsigned int i = 0; i < num_platforms; ++i)
        {
            char* profile = get_platform_info(platforms[i],
                                              CL_PLATFORM_PROFILE);
            printf("Platform profile: %s\n", profile);
            delete [] profile;
            char* version = get_platform_info(platforms[i],
                                              CL_PLATFORM_VERSION);
            printf("Platform version: %s\n", version);
            delete [] version;
            char* name = get_platform_info(platforms[i],
                                           CL_PLATFORM_NAME);
            printf("Platform name: %s\n", name);
            delete [] name;
            char* vendor = get_platform_info(platforms[i],
                                             CL_PLATFORM_VENDOR);
            printf("Platform vendor: %s\n", vendor);
            delete [] vendor;
            char* extensions = get_platform_info(platforms[i],
                                                 CL_PLATFORM_EXTENSIONS);
            printf("Platform extensions: %s\n", extensions);
            delete [] extensions;
        }
        return platforms;
    }
    else
    {
        exit(1);
    }
}

cl_context get_context(cl_device_id device)
{
    cl_int err;
    cl_context context = clCreateContext(0, 1,
                                         &device, 0, 0,
                                         &err);
    check_error(err);
    return context;
}

char* get_device_info(cl_device_id device,
                      cl_device_info info)
{
    size_t size;
    check_error(clGetDeviceInfo(device,
                                info,
                                0, 0, &size));
    char* profile = new char[size];
    check_error(clGetDeviceInfo(device,
                                info,
                                size, profile, 0));
    return profile;
}

cl_device_id* get_devices(cl_platform_id platform)
{
    cl_uint num_devices;
    cl_uint device_flag = CL_DEVICE_TYPE_CPU;
    cl_int err = clGetDeviceIDs(platform,
                                device_flag,
                                0, 0, &num_devices);
    if (err == CL_DEVICE_NOT_FOUND)
    {
        printf("Error finding an avaliable GPU, trying to search for a CPU\n");
        device_flag = CL_DEVICE_TYPE_ALL;
        check_error(clGetDeviceIDs(platform,
                                   device_flag,
                                   0, 0, &num_devices));
    }
    printf("Found %d devices\n", num_devices);
    if (num_devices > 0)
    {
        cl_device_id* devices = new cl_device_id[num_devices];
        check_error(clGetDeviceIDs(platform,
                                   device_flag,
                                   num_devices, devices, 0));
        for (unsigned int i = 0; i < num_devices; ++i)
        {
            char* name = get_device_info(devices[i],
                                         CL_DEVICE_NAME);
            printf("Device[%d] name: %s\n", i, name);
            delete [] name;
            char* vendor = get_device_info(devices[i],
                                           CL_DEVICE_VENDOR);
            printf("Device[%d] vendor: %s\n", i, vendor);
            delete [] vendor;
            char* version = get_device_info(devices[i],
                                            CL_DEVICE_VERSION);
            printf("Device[%d] device version: %s\n", i, version);
            delete [] version;
            char* driver_version = get_device_info(devices[i],
                                                   CL_DRIVER_VERSION);
            printf("Device[%d] driver version: %s\n", i, driver_version);
            delete [] driver_version;
            char* extensions = get_device_info(devices[i],
                                               CL_DEVICE_EXTENSIONS);
            printf("Device[%d] extensions: %s\n", i, extensions);
            delete [] extensions;
        }
        return devices;
    }
    else
    {
        exit(1);
    }
}

void setup_command_queue(ocl_mem_t* ocl_mem)
{
    cl_int err;
    cl_command_queue queue = clCreateCommandQueue(ocl_mem->context,
                             ocl_mem->devices[0],
                             CL_QUEUE_PROFILING_ENABLE,
                             &err);
    check_error(err);
    ocl_mem->queue = queue;
}
