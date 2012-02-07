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
#include "milkyway_math.h"
#include "milkyway_cl_show_types.h"
#include "milkyway_cl_setup.h"
#include "milkyway_cl_program.h"
#include "milkyway_cl_device.h"
#include "milkyway_cl_util.h"

static void CL_CALLBACK contextCallback(const char* errInfo,
                                        const void* privateInfo,
                                        size_t cb,
                                        void* userData)
{
    (void) privateInfo, (void) cb, (void) userData;
    mw_printf("CL context error: %s\n", errInfo);
}

#if defined(cl_APPLE_ContextLoggingFunctions) && defined(__APPLE__)
  #define MW_CONTEXT_LOGGER clLogMessagesToStderrAPPLE
#else
  #define MW_CONTEXT_LOGGER contextCallback
#endif /* cl_APPLE_ContextLoggingFunctions */

static cl_int mwCreateCtxQueue(CLInfo* ci, cl_bool useBufQueue, cl_bool enableProfiling)
{
    cl_int err = CL_SUCCESS;
    cl_command_queue_properties props = enableProfiling ? CL_QUEUE_PROFILING_ENABLE : 0;

    ci->clctx = clCreateContext(NULL, 1, &ci->dev, MW_CONTEXT_LOGGER, NULL, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating context");
        return err;
    }

    ci->queue = clCreateCommandQueue(ci->clctx, ci->dev, props, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating command queue");
        return err;
    }

    if (useBufQueue)
    {
        ci->bufQueue = clCreateCommandQueue(ci->clctx, ci->dev, 0, &err);
        if (err != CL_SUCCESS)
        {
            mwPerrorCL(err, "Error creating buffer command queue");
            return err;
        }
    }

    return CL_SUCCESS;
}

/* Return CL_UINT_MAX if it doesn't find one */
static cl_uint choosePlatform(const char* prefVendor, const cl_platform_id* platforms, cl_uint nPlatform)
{
    cl_uint i;
    char platVendor[256];
    char prefBuf[256];
    cl_int err;

    if (!platforms || nPlatform == 0)
        return CL_UINT_MAX;

    /* No strnstr() on Windows, also be paranoid */
    memset(prefBuf, 0, sizeof(prefBuf));
    strncpy(prefBuf, prefVendor, sizeof(prefBuf));
    prefBuf[sizeof(prefBuf) - 1] = '\0';

    /* Out of the available platforms, see if one has a matching vendor */
    for (i = 0; i < nPlatform; ++i)
    {
        err = clGetPlatformInfo(platforms[i],
                                CL_PLATFORM_VENDOR,
                                sizeof(platVendor),
                                platVendor,
                                NULL);
        if (err != CL_SUCCESS)
        {
            mwPerrorCL(err, "Error getting platform vendor");
            return CL_UINT_MAX;
        }

        if (strstr(platVendor, prefBuf))
        {
            return i;
        }
    }

    return CL_UINT_MAX;
}

static cl_int mwGetCLInfo(CLInfo* ci, const CLRequest* clr)
{
    cl_int err = CL_SUCCESS;
    cl_uint nPlatform = 0;
    cl_uint nDev = 0;
    cl_platform_id* ids;
    cl_device_id* devs;
    cl_uint platformChoice = 0;

    memset(ci, 0, sizeof(*ci));

    ids = mwGetAllPlatformIDs(&nPlatform);
    if (!ids)
        return MW_CL_ERROR;

    if (mwIsFirstRun())
    {
        mwPrintPlatforms(ids, nPlatform);
    }

    /* We have this set by default to UINT_MAX, so if it's in a
     * legitimate range, it was specified */
    if (clr->platform < nPlatform)
    {
        platformChoice = clr->platform;
    }
    else if (clr->preferredPlatformVendor && strcmp(clr->preferredPlatformVendor, "")) /* Didn't specify platform by index, try picking one by name */
    {
        platformChoice = choosePlatform(clr->preferredPlatformVendor, ids, nPlatform);
    }
    else
    {
        platformChoice = 0;
    }

    if (platformChoice >= nPlatform)
    {
        mw_printf("Didn't find preferred platform\n");
        platformChoice = 0;
    }

    mw_printf("Using device %u on platform %u\n", clr->devNum, platformChoice);

    ci->plat = ids[platformChoice];
    devs = mwGetAllDevices(ci->plat, &nDev);
    if (!devs)
    {
        free(ids);
        return MW_CL_ERROR;
    }

    err = mwSelectDevice(ci, devs, clr, nDev);
    free(ids);
    free(devs);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to select a device");
        return err;
    }

    return CL_SUCCESS;
}


/* Query one device specified by type, create a context and command
 * queue, as well as retrieve device information
 */
cl_int mwSetupCL(CLInfo* ci, const CLRequest* clr)
{
    cl_int err;

    err = mwGetCLInfo(ci, clr);
    if (err != CL_SUCCESS)
    {
        mw_printf("Failed to get information about device\n");
        return err;
    }

    err = mwGetDevInfo(&ci->di, ci->dev);
    if (err != CL_SUCCESS)
    {
        mw_printf("Failed to get device info\n");
        return err;
    }

    if (mwIsFirstRun())
    {
        if (clr->verbose)
        {
            mwPrintDevInfo(&ci->di);
        }
        else
        {
            mwPrintDevInfoShort(&ci->di);
        }
    }

    if (clr->pollingMode <= MW_POLL_WORKAROUND_CL_WAIT_FOR_EVENTS)
    {
        /* With the default, we will try to use clWaitForEvents()
         * unless we know the driver is bad and busy waits, in which
         * case we will do it manually.
         */

        if (mwDriverHasHighCPUWaitIssue(ci))
        {
            ci->pollingMode = MW_POLL_SLEEP_CL_WAIT_FOR_EVENTS;
        }
        else
        {
            ci->pollingMode = MW_POLL_CL_WAIT_FOR_EVENTS;
        }
    }
    else
    {
        ci->pollingMode = clr->pollingMode;
    }

    return mwCreateCtxQueue(ci, CL_FALSE, clr->enableProfiling);
}

cl_int mwDestroyCLInfo(CLInfo* ci)
{
    cl_int err = CL_SUCCESS;
  /* Depending on where things fail, some of these will be NULL, and
   * will spew errors when trying to cleanup. */
    if (ci->queue)
        err |= clReleaseCommandQueue(ci->queue);
    if (ci->bufQueue)
        err |= clReleaseCommandQueue(ci->bufQueue);
    if (ci->clctx)
        err |= clReleaseContext(ci->clctx);

    /* TODO: or'ing the err and showing = useless */
    if (err)
        mwPerrorCL(err, "Error cleaning up CLInfo");

    return err;
}

