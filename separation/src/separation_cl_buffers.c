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

static inline cl_mem createReadWriteBuffer(cl_context clctx, size_t size, cl_int* err)
{
    return clCreateBuffer(clctx, CL_MEM_READ_WRITE, size, NULL, err);
}

static inline cl_mem createZeroReadWriteBuffer(CLInfo* ci, size_t size, cl_int* errOut)
{
    void* p;
    cl_mem mem = NULL;
    cl_int err = CL_SUCCESS;

    mem = clCreateBuffer(ci->clctx, CL_MEM_READ_WRITE, size, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Failed to create zero buffer: %s\n", showCLInt(err));
        goto fail;
    }

    p = clEnqueueMapBuffer(ci->queue, mem, CL_TRUE, CL_MAP_WRITE,
                           0, size, 0, NULL, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error mapping zero buffer: %s\n", showCLInt(err));
        goto fail;
    }

    memset(p, 0, size);

    err = clEnqueueUnmapMemObject(ci->queue, mem, p, 0, NULL, NULL);
    if (err != CL_SUCCESS)
        warn("Failed to unmap zero buffer: %s\n", showCLInt(err));

fail:
    *errOut = err;
    return mem;
}

static inline cl_int createOutMuBuffer(CLInfo* ci,
                                       SeparationCLMem* cm,
                                       const SeparationSizes* sizes)
{
    cl_int err;

    cm->outMu = createZeroReadWriteBuffer(ci, sizes->outMu, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out mu buffer of size %llu: %s\n", (cl_ulong) sizes->outMu, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createOutProbsBuffer(CLInfo* ci,
                                          SeparationCLMem* cm,
                                          const SeparationSizes* sizes)
{
    cl_int err;

    cm->outProbs = createZeroReadWriteBuffer(ci, sizes->outProbs, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out probs buffer of size %llu: %s\n", (cl_ulong) sizes->outProbs, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createSCBuffer(CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const StreamConstants* sc,
                                    const SeparationSizes* sizes,
                                    const cl_mem_flags constBufFlags)
{
    cl_int err;

    cm->sc = clCreateBuffer(ci->clctx,
                            constBufFlags,
                            sizes->sc,
                            (void*) sc,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream constants buffer of size %llu: %s\n", (cl_ulong) sizes->sc, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createAPBuffer(CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const AstronomyParameters* ap,
                                    const SeparationSizes* sizes,
                                    const cl_mem_flags constBufFlags)
{
    cl_int err;
    cm->ap = clCreateBuffer(ci->clctx, constBufFlags, sizes->ap, (void*) ap, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating astronomy parameters buffer of size %llu: %s\n", (cl_ulong) sizes->ap, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createIABuffer(CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const IntegralArea* ia,
                                    const SeparationSizes* sizes,
                                    const cl_mem_flags constBufFlags)
{
    cl_int err;

    cm->ia = clCreateBuffer(ci->clctx, constBufFlags, sizes->ia, (void*) ia, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating integral area buffer of size %llu: %s\n", (cl_ulong) sizes->ia, showCLInt(err));
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
                             cl_bool useImages)
{
    cl_int err;
    RPoints* r_pts;
    RConsts* rc;
    const cl_mem_flags constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
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
        warn("Error creating stream r points buffer of size %llu: %s\n", (cl_ulong) sizes->rPts, showCLInt(err));
        return err;
    }

    cm->rc = clCreateBuffer(ci->clctx, constBufFlags, sizes->rc, rc, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream r consts buffer of size %llu: %s\n", (cl_ulong) sizes->rc, showCLInt(err));
        return err;
    }

    cm->sg_dx = clCreateBuffer(ci->clctx, constBufFlags, sizes->sg_dx, sg.dx, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream sg_dx buffer of size %llu: %s\n", (cl_ulong) sizes->sg_dx, showCLInt(err));
        return err;
    }

    mwAlignedFree(r_pts);
    mwAlignedFree(rc);

    return CL_SUCCESS;
}

static cl_int createLBTrigBuffer(CLInfo* ci,
                                 SeparationCLMem* cm,
                                 const AstronomyParameters* ap,
                                 const IntegralArea* ia,
                                 const SeparationSizes* sizes)
{
    cl_int err;
    LBTrig* lbts;
    const cl_mem_flags constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    lbts = precalculateLBTrig(ap, ia);
    cm->lbts = clCreateBuffer(ci->clctx, constBufFlags, sizes->lbts, lbts, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating lb_trig buffer of size %llu: %s\n", (cl_ulong) sizes->lbts, showCLInt(err));
        return err;
    }

    mwAlignedFree(lbts);

    return CL_SUCCESS;
}

void calculateSizes(SeparationSizes* sizes, const AstronomyParameters* ap, const IntegralArea* ia)
{
    sizes->outMu = sizeof(real) * ia->mu_steps * ia->r_steps;
    sizes->outProbs = sizeof(real) * ia->mu_steps * ia->r_steps * ap->number_streams;
    sizes->ap = sizeof(AstronomyParameters);
    sizes->ia = sizeof(IntegralArea);
    sizes->sc = sizeof(StreamConstants) * ap->number_streams;
    sizes->rPts = sizeof(RPoints) * ap->convolve * ia->r_steps;
    sizes->rc = sizeof(RConsts) * ia->r_steps;
    sizes->sg_dx = sizeof(real) * ap->convolve;
    sizes->lbts = sizeof(LBTrig) * ia->mu_steps * ia->nu_steps;
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
    cl_mem_flags constBufFlags;

    if (ci->devType == CL_DEVICE_TYPE_CPU)
        constBufFlags = CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR;
    else
        constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    err |= createOutMuBuffer(ci, cm, sizes);
    err |= createOutProbsBuffer(ci, cm, sizes);
    err |= createAPBuffer(ci, cm, ap, sizes, constBufFlags);
    err |= createIABuffer(ci, cm, ia, sizes, constBufFlags);
    err |= createSCBuffer(ci, cm, sc, sizes, constBufFlags);
    err |= createRBuffers(ci, cm, ap, ia, sg, sizes, useImages);
    err |= createLBTrigBuffer(ci, cm, ap, ia, sizes);

    return err;
}

void releaseSeparationBuffers(SeparationCLMem* cm)
{
    clReleaseMemObject(cm->outProbs);
    clReleaseMemObject(cm->outMu);

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
    real* mapOutMu;

    mapOutMu = (real*) clEnqueueMapBuffer(ci->queue,
                                          cm->outMu,
                                          CL_TRUE, CL_MAP_READ,
                                          0, resultsSize,
                                          0, NULL,
                                          NULL,
                                          &err);
    if (err != CL_SUCCESS)
        warn("Error mapping integral result buffer: %s\n", showCLInt(err));

    return mapOutMu;
}

real* mapProbsResults(CLInfo* ci, SeparationCLMem* cm, size_t probsResultsSize)
{
    cl_int err;
    real* mapOutProbs;

    mapOutProbs = (real*) clEnqueueMapBuffer(ci->queue,
                                             cm->outProbs,
                                             CL_TRUE, CL_MAP_READ,
                                             0, probsResultsSize,
                                             0, NULL,
                                             NULL,
                                             &err);
    if (err != CL_SUCCESS)
        warn("Error mapping probs result buffer: %s\n", showCLInt(err));

    return mapOutProbs;
}

