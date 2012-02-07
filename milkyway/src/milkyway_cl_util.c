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

#include <stdarg.h>

#include "milkyway_util.h"
#include "milkyway_cl_show_types.h"
#include "milkyway_cl_util.h"
#include "milkyway_cl_device.h"
#include "milkyway_cl_setup.h"

void mwPerrorCL(cl_int err, const char* fmt, ...)
{
    va_list argPtr;

    va_start(argPtr, fmt);
    vfprintf(stderr, fmt, argPtr);
    va_end(argPtr);

    fprintf(stderr, " (%d): %s\n", err, showCLInt(err));
}


/* Timing in nanoseconds */
cl_ulong mwEventTimeNS(cl_event ev)
{
    cl_int err;
    cl_ulong ts, te;

    err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &ts, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to get event start time");
        return 0;
    }

    err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &te, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to get event end time");
        return 0;
    }

    assert(te >= ts);

    return te - ts;
}

/* Timing in seconds */
double mwEventTime(cl_event ev)
{
    return (double) mwEventTimeNS(ev) * 1.0e-9;
}

double mwEventTimeMS(cl_event ev)
{
    return (double) mwEventTimeNS(ev) * 1.0e-6;
}

double mwReleaseEventWithTiming(cl_event ev)
{
    double t;

    t = mwEventTime(ev);
    clReleaseEvent(ev);

    return t;
}

/* Wait for an event then release it */
cl_int mwWaitReleaseEvent(cl_event* ev)
{
    cl_int err;

    assert(ev);

    err = clWaitForEvents(1, ev);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to wait for event");
        return err;
    }

    err = clReleaseEvent(*ev);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to release event");
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
        mwPerrorCL(err, "Failed to create custom event");
        return NULL;
    }

    return ev;
}

cl_int mwFinishEvent(cl_event ev)
{
    cl_int err;

    err = clSetUserEventStatus(ev, CL_COMPLETE);
    if (err != CL_SUCCESS)
        mwPerrorCL(err, "Failed to mark custom event as completed");

    return err;
}

#endif /* CL_VERSION_1_1 */

cl_int mwGetWorkGroupInfo(cl_kernel kern, const CLInfo* ci, WGInfo* wgi)
{
    cl_int err = CL_SUCCESS;

    err |= clGetKernelWorkGroupInfo(kern,
                                    ci->dev,
                                    CL_KERNEL_COMPILE_WORK_GROUP_SIZE,
                                    sizeof(wgi->cwgs),
                                    wgi->cwgs,
                                    NULL);

    err |= clGetKernelWorkGroupInfo(kern,
                                    ci->dev,
                                    CL_KERNEL_WORK_GROUP_SIZE,
                                    sizeof(size_t),
                                    &wgi->wgs,
                                    NULL);

    err |= clGetKernelWorkGroupInfo(kern,
                                    ci->dev,
                                    CL_KERNEL_LOCAL_MEM_SIZE,
                                    sizeof(cl_ulong),
                                    &wgi->lms,
                                    NULL);
    if (err != CL_SUCCESS)
        mwPerrorCL(err, "Failed to get kernel work group info");

    return err;
}

void mwPrintWorkGroupInfo(const WGInfo* wgi)
{
    mw_printf("Kernel work group info:\n"
              "  Work group size = "ZU"\n"
              "  Kernel local mem size = "LLU"\n"
              "  Compile work group size = { "ZU", "ZU", "ZU" }\n",
              wgi->wgs,
              wgi->lms,
              wgi->cwgs[0], wgi->cwgs[1], wgi->cwgs[2]);
}

cl_mem mwCreateZeroReadWriteBuffer(CLInfo* ci, size_t size)
{
    void* p;
    cl_mem mem = NULL;
    cl_int err = CL_SUCCESS;

    mem = clCreateBuffer(ci->clctx, CL_MEM_READ_WRITE, size, NULL, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to create zero buffer of size "ZU, size);
        goto fail;
    }

    p = clEnqueueMapBuffer(ci->queue, mem, CL_TRUE, CL_MAP_WRITE,
                           0, size, 0, NULL, NULL, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error mapping zero buffer");
        goto fail;
    }

    memset(p, 0, size);

    err = clEnqueueUnmapMemObject(ci->queue, mem, p, 0, NULL, NULL);
    if (err != CL_SUCCESS)
        mwPerrorCL(err, "Failed to unmap zero buffer");
fail:
    return mem;
}


cl_mem mwDuplicateBuffer(CLInfo* ci, cl_mem buf)
{
    cl_mem bufCopy;
    size_t size;
    cl_mem_flags flags;
    cl_int err;
    cl_event ev;

    if (!buf)
    {
        return NULL;
    }

    err = clGetMemObjectInfo(buf, CL_MEM_FLAGS, sizeof(flags), &flags, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to get memory flags for buffer duplication");
        return NULL;
    }

    err = clGetMemObjectInfo(buf, CL_MEM_SIZE, sizeof(size), &size, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to get memory size for buffer duplication");
        return NULL;
    }

    /* We may have initialized that one from a host pointer, but not this one */
    flags ^= CL_MEM_COPY_HOST_PTR;

    /* Create a copy of the same size */
    bufCopy = clCreateBuffer(ci->clctx, flags, size, NULL, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to create copy buffer of size "ZU, size);
        return NULL;
    }

    err = clEnqueueCopyBuffer(ci->queue, buf, bufCopy, 0, 0, size, 0, NULL, &ev);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to enqueue buffer copy of size"ZU, size);
        clReleaseMemObject(bufCopy);
        return NULL;
    }

    err = mwWaitReleaseEvent(&ev);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to wait for buffer copy");
        clReleaseMemObject(bufCopy);
        return NULL;
    }

    return bufCopy;
}

/* Optionally manually poll for event completion.  This is because some
   versions of both the Nvidia and AMD drivers developed an issue
   where clWaitForEvents/clFinish would cause 100% CPU usage.
*/
cl_int mwCLWaitForEvent(CLInfo* ci, cl_event ev, cl_uint initialWait)
{
    cl_int err;
    cl_int pollingMode = ci->pollingMode;

    if (pollingMode <= MW_POLL_CL_WAIT_FOR_EVENTS)
    {
        return clWaitForEvents(1, &ev);
    }
    else if (pollingMode == MW_POLL_SLEEP_CL_WAIT_FOR_EVENTS)
    {
        err = clFlush(ci->queue);
        if (err != CL_SUCCESS)
            return err;

        mwMilliSleep(initialWait);
        return clWaitForEvents(1, &ev);
    }
    else  /* Manually poll for pollingMode milliseconds */
    {
        cl_int execStatus;

        err = clFlush(ci->queue);  /* Make sure the task is submitted before we wait for it */
        if (err != CL_SUCCESS)
            return err;

        mwMilliSleep(initialWait);

        do
        {
            err = clGetEventInfo(ev,
                                 CL_EVENT_COMMAND_EXECUTION_STATUS,
                                 sizeof(cl_int),
                                 &execStatus,
                                 NULL);
            if (err != CL_SUCCESS)
                return err;

            mwMilliSleep(pollingMode);
        }
        while (execStatus != CL_COMPLETE);

        return CL_SUCCESS;
    }
}

cl_bool mwDriverHasHighCPUWaitIssue(CLInfo* ci)
{
    const DevInfo* di = &ci->di;
    int major = 0, minor = 0, patchLevel = 0;

    if (mwIsNvidiaGPUDevice(di))
    {
        if (sscanf(di->driver, "%u.%u", &major, &minor) != 2)
        {
            return CL_FALSE;
        }

        /* Issue started around 270.xx and hasn't been fixed yet. */
        return (major >= 270);
    }
    else if (mwIsAMDGPUDevice(di))
    {
        /* Sometimes it has other stuff like (VM) after. Not sure what that means. */
        if (sscanf(di->driver, "CAL %d.%d.%d", &major, &minor, &patchLevel) != 3)
        {
            return CL_FALSE;
        }

        /* I think it happened in 11.7 and 11.8 */
        return (major == 1 && minor == 4 && patchLevel >= 1457 && patchLevel < 1546);
    }
    else
    {
        return CL_FALSE;
    }
}

