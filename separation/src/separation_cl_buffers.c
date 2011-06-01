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
#include "mw_cl.h"
#include "separation_cl_buffers.h"
#include "r_points.h"
#include "calculated_constants.h"

static cl_mem createReadWriteBuffer(cl_context clctx, size_t size, cl_int* err)
{
    return clCreateBuffer(clctx, CL_MEM_READ_WRITE, size, NULL, err);
}

static cl_mem createZeroReadWriteBuffer(CLInfo* ci, size_t size, cl_int* errOut)
{
    void* p;
    cl_mem mem = NULL;
    cl_int err = CL_SUCCESS;

    mem = clCreateBuffer(ci->clctx, CL_MEM_READ_WRITE, size, NULL, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Failed to create zero buffer", err);
        goto fail;
    }

    p = clEnqueueMapBuffer(ci->queue, mem, CL_TRUE, CL_MAP_WRITE,
                           0, size, 0, NULL, NULL, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error mapping zero buffer", err);
        goto fail;
    }

    memset(p, 0, size);

    err = clEnqueueUnmapMemObject(ci->queue, mem, p, 0, NULL, NULL);
    if (err != CL_SUCCESS)
        mwCLWarn("Failed to unmap zero buffer", err);

fail:
    *errOut = err;
    return mem;
}

static cl_int createOutBgBuffer(CLInfo* ci,
                                SeparationCLMem* cm,
                                const SeparationSizes* sizes)
{
    cl_int err;

    cm->outBg = createZeroReadWriteBuffer(ci, sizes->outBg, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating out bg buffer of size "ZU, err, sizes->outBg);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int createOutStreamsBuffer(CLInfo* ci, SeparationCLMem* cm, const SeparationSizes* sizes)
{
    cl_int err;

    cm->outStreams = createZeroReadWriteBuffer(ci, sizes->outStreams, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating out probs buffer of size "ZU, err, sizes->outStreams);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int createSCBuffer(CLInfo* ci,
                             SeparationCLMem* cm,
                             const StreamConstants* sc,
                             const SeparationSizes* sizes,
                             const cl_mem_flags constBufFlags)
{
    cl_int err;

    cm->sc = clCreateBuffer(ci->clctx, constBufFlags, sizes->sc, (void*) sc, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating stream constants buffer of size "ZU, err, sizes->sc);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int createAPBuffer(CLInfo* ci,
                             SeparationCLMem* cm,
                             const AstronomyParameters* ap,
                             const SeparationSizes* sizes,
                             const cl_mem_flags constBufFlags)
{
    cl_int err;
    cm->ap = clCreateBuffer(ci->clctx, constBufFlags, sizes->ap, (void*) ap, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating astronomy parameters buffer of size "ZU, err, sizes->ap);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int createIABuffer(CLInfo* ci,
                             SeparationCLMem* cm,
                             const IntegralArea* ia,
                             const SeparationSizes* sizes,
                             const cl_mem_flags constBufFlags)
{
    cl_int err;

    cm->ia = clCreateBuffer(ci->clctx, constBufFlags, sizes->ia, (void*) ia, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating integral area buffer of size "ZU, err, sizes->ia);
        return err;
    }

    return CL_SUCCESS;
}

static cl_int createRBuffers(CLInfo* ci,
                             SeparationCLMem* cm,
                             const AstronomyParameters* ap,
                             const IntegralArea* ia,
                             const StreamGauss sg,
                             const SeparationSizes* sizes,
                             cl_mem_flags constBufFlags,
                             cl_bool useImages)
{
    cl_int err;
    RPoints* r_pts;
    RConsts* rc;
    cl_image_format format = { CL_RGBA, CL_UNSIGNED_INT32 };

    r_pts = precalculateRPts(ap, ia, sg, &rc, useImages);

    if (useImages)
    {
        cm->rPts = clCreateImage2D(ci->clctx, constBufFlags, &format,
                                   ia->r_steps, ap->convolve,
                                   0, r_pts, &err);
    }
    else
    {
        cm->rPts = clCreateBuffer(ci->clctx, constBufFlags, sizes->rPts, r_pts, &err);
    }

    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating stream r points buffer of size "ZU, err, sizes->rPts);
        return err;
    }

    cm->rc = clCreateBuffer(ci->clctx, constBufFlags, sizes->rc, rc, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating stream r consts buffer of size "ZU, err, sizes->rc);
        return err;
    }

    cm->sg_dx = clCreateBuffer(ci->clctx, constBufFlags, sizes->sg_dx, sg.dx, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating stream sg_dx buffer of size "ZU, err, sizes->sg_dx);
        return err;
    }

    mwFreeA(r_pts);
    mwFreeA(rc);

    return CL_SUCCESS;
}

static cl_int createLBTrigBuffer(CLInfo* ci,
                                 SeparationCLMem* cm,
                                 const AstronomyParameters* ap,
                                 const IntegralArea* ia,
                                 const SeparationSizes* sizes,
                                 const cl_mem_flags constBufFlags)
{
    cl_int err;
    LBTrig* lbts;

    lbts = precalculateLBTrig(ap, ia, CL_FALSE);
    cm->lbts = clCreateBuffer(ci->clctx, constBufFlags, sizes->lbts, lbts, &err);
    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating lb_trig buffer of size "ZU, err, sizes->lbts);
        return err;
    }

    mwFreeA(lbts);

    return CL_SUCCESS;
}

void calculateSizes(SeparationSizes* sizes, const AstronomyParameters* ap, const IntegralArea* ia)
{
    /* globals */
    sizes->outBg = sizeof(real) * ia->mu_steps * ia->r_steps;
    sizes->outStreams = sizeof(real) * ia->mu_steps * ia->r_steps * ap->number_streams;

    sizes->rPts = sizeof(RPoints) * ap->convolve * ia->r_steps;
    sizes->lbts = sizeof(LBTrig) * ia->mu_steps * ia->nu_steps;

    /* Constant buffer things */
    sizes->ap = sizeof(AstronomyParameters);
    sizes->ia = sizeof(IntegralArea);
    sizes->sc = sizeof(StreamConstants) * ap->number_streams;
    sizes->rc = sizeof(RConsts) * ia->r_steps;
    sizes->sg_dx = sizeof(real) * ap->convolve;
}

cl_int createSeparationBuffers(CLInfo* ci,
                               SeparationCLMem* cm,
                               const AstronomyParameters* ap,
                               const IntegralArea* ia,
                               const StreamConstants* sc,
                               const StreamGauss sg,
                               const SeparationSizes* sizes,
                               cl_bool useImages)  /* Use images for some buffers if wanted / available. */
{
    cl_int err = CL_SUCCESS;
    cl_mem_flags constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    err |= createOutBgBuffer(ci, cm, sizes);
    err |= createOutStreamsBuffer(ci, cm, sizes);
    err |= createAPBuffer(ci, cm, ap, sizes, constBufFlags);
    err |= createIABuffer(ci, cm, ia, sizes, constBufFlags);
    err |= createSCBuffer(ci, cm, sc, sizes, constBufFlags);
    err |= createRBuffers(ci, cm, ap, ia, sg, sizes, constBufFlags, useImages);
    err |= createLBTrigBuffer(ci, cm, ap, ia, sizes, constBufFlags);

    return err;
}

void releaseSeparationBuffers(SeparationCLMem* cm)
{
    clReleaseMemObject(cm->outStreams);
    clReleaseMemObject(cm->outBg);

    clReleaseMemObject(cm->ap);
    clReleaseMemObject(cm->ia);
    clReleaseMemObject(cm->sc);
    clReleaseMemObject(cm->rPts);
    clReleaseMemObject(cm->rc);
    clReleaseMemObject(cm->sg_dx);
    clReleaseMemObject(cm->lbts);
}

real* mapIntegralResults(CLInfo* ci, SeparationCLMem* cm, size_t resultsSize)
{
    cl_int err;
    real* mapOutBg;

    mapOutBg = (real*) clEnqueueMapBuffer(ci->queue,
                                          cm->outBg,
                                          CL_TRUE, CL_MAP_READ,
                                          0, resultsSize,
                                          0, NULL,
                                          NULL,
                                          &err);
    if (err != CL_SUCCESS)
        mwCLWarn("Error mapping integral result buffer", err);

    return mapOutBg;
}

real* mapStreamsResults(CLInfo* ci, SeparationCLMem* cm, size_t streamsResultsSize)
{
    cl_int err;
    real* mapOutStreams;

    mapOutStreams = (real*) clEnqueueMapBuffer(ci->queue,
                                               cm->outStreams,
                                               CL_TRUE, CL_MAP_READ,
                                               0, streamsResultsSize,
                                               0, NULL,
                                               NULL,
                                               &err);
    if (err != CL_SUCCESS)
        mwCLWarn("Error mapping stream result buffer", err);

    return mapOutStreams;
}

