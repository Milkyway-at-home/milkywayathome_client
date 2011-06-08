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

#include "milkyway_util.h"
#include "mw_cl_show_types.h"
#include "mw_cl_util.h"
#include "mw_cl_setup.h"

/* Timing in nanoseconds */
cl_ulong mwEventTimeNS(cl_event ev)
{
    cl_int err;
    cl_ulong ts, te;

    err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &ts, NULL);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to get event start time", err);
        return 0;
    }

    err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &te, NULL);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to get event end time", err);
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
        mwCLWarn("%s: Failed to wait for event", err, FUNC_NAME);
        return err;
    }

    err = clReleaseEvent(*ev);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("%s: Failed to release event", err, FUNC_NAME);
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
        mwCLWarn("Failed to create custom event", err);
        return NULL;
    }

    return ev;
}

cl_int mwFinishEvent(cl_event ev)
{
    cl_int err;

    err = clSetUserEventStatus(ev, CL_COMPLETE);
    if (err != CL_SUCCESS)
        mwCLWarn("Failed to mark custom event as completed", err);

    return err;
}

#endif /* CL_VERSION_1_1 */

cl_int mwGetWorkGroupInfo(const CLInfo* ci, WGInfo* wgi)
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
        mwCLWarn("Failed to get kernel work group info", err);

    return err;
}

void mwPrintWorkGroupInfo(const WGInfo* wgi)
{
    warn("Kernel work group info:\n"
         "  Work group size = "ZU"\n"
         "  Kernel local mem size = "LLU"\n"
         "  Compile work group size = { "ZU", "ZU", "ZU" }\n",
         wgi->wgs,
         wgi->lms,
         wgi->cwgs[0], wgi->cwgs[1], wgi->cwgs[2]);
}

