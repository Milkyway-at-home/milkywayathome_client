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

#include "milkyway_util.h"
#include "milkyway_math.h"
#include "mw_cl_show_types.h"
#include "mw_cl_setup.h"
#include "mw_cl_program.h"
#include "mw_cl_device.h"
#include "mw_cl_util.h"

static cl_int mwCreateCtxQueue(CLInfo* ci, cl_bool useBufQueue)
{
    cl_int err = CL_SUCCESS;

    ci->clctx = clCreateContext(NULL, 1, &ci->dev, NULL, NULL, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating context", err);
        return err;
    }

    ci->queue = clCreateCommandQueue(ci->clctx, ci->dev, 0, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating command queue", err);
        return err;
    }

    if (useBufQueue)
    {
        ci->bufQueue = clCreateCommandQueue(ci->clctx, ci->dev, 0, &err);
        if (err != CL_SUCCESS)
        {
            mwCLWarn("Error creating buffer command queue", err);
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
            mwCLWarn("Error getting platform vendor", err);
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
    cl_uint platformChoice;

    ids = mwGetAllPlatformIDs(ci, &nPlatform);
    if (!ids)
        return MW_CL_ERROR;

    mwPrintPlatforms(ids, nPlatform);

    /* We have this set by default to UINT_MAX, so if it's in a
     * legitimate range, it was specified */
    if (clr->platform < nPlatform)
    {
        platformChoice = clr->platform;
    }
    else if (strcmp(clr->preferredPlatformVendor, "")) /* Didn't specify platform by index, try picking one by name */
    {
        platformChoice = choosePlatform(clr->preferredPlatformVendor, ids, nPlatform);
    }
    else
    {
        platformChoice = 0;
    }

    if (platformChoice >= nPlatform)
    {
        warn("Didn't find preferred platform\n");
        platformChoice = 0;
    }

    warn("Using device %u on platform %u\n", clr->devNum, platformChoice);

    devs = mwGetAllDevices(ids[platformChoice], &nDev);
    if (!devs)
    {
        warn("Error getting devices\n");
        free(ids);
        return MW_CL_ERROR;
    }

    err = mwSelectDevice(ci, devs, clr, nDev);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to select a device", err);
        err = -1;
    }

    free(ids);
    free(devs);

    return err;
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
        warn("Failed to get information about device\n");
        return err;
    }

    err = mwCreateCtxQueue(ci, CL_FALSE);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating CL context and command queue", err);
        return err;
    }

    return CL_SUCCESS;
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
    if (ci->prog)
        err |= clReleaseProgram(ci->prog);
    if (ci->kern)
        err |= clReleaseKernel(ci->kern);
    if (ci->clctx)
        err |= clReleaseContext(ci->clctx);

    /* TODO: or'ing the err and showing = useless */
    if (err)
        mwCLWarn("Error cleaning up CLInfo", err);

    return err;
}

