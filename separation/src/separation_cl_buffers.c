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
        warn("Error creating out mu buffer of size %zu: %s\n", sizes->outMu, showCLInt(err));
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
        warn("Error creating out probs buffer of size %zu: %s\n", sizes->outProbs, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createSCBuffer(CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const STREAM_CONSTANTS* sc,
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
        warn("Error creating stream constants buffer of size %zu: %s\n", sizes->sc, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createAPBuffer(CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const ASTRONOMY_PARAMETERS* ap,
                                    const SeparationSizes* sizes,
                                    const cl_mem_flags constBufFlags)
{
    cl_int err;
    cm->ap = clCreateBuffer(ci->clctx, constBufFlags, sizes->ap, (void*) ap, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating astronomy parameters buffer of size %zu: %s\n", sizes->ap, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createIABuffer(CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const INTEGRAL_AREA* ia,
                                    const SeparationSizes* sizes,
                                    const cl_mem_flags constBufFlags)
{
    cl_int err;

    cm->ia = clCreateBuffer(ci->clctx, constBufFlags, sizes->ia, (void*) ia, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating integral area buffer of size %zu: %s\n", sizes->ia, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createRBuffers(CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const ASTRONOMY_PARAMETERS* ap,
                                    const INTEGRAL_AREA* ia,
                                    const STREAM_GAUSS sg,
                                    const SeparationSizes* sizes)
{
    cl_int err;

    R_POINTS* r_pts;
    R_CONSTS* rc;
    const cl_mem_flags constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    r_pts = precalculate_r_pts(ap, ia, sg, &rc);

    cm->rPts = clCreateBuffer(ci->clctx, constBufFlags, sizes->rPts, r_pts, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream r points buffer of size %zu: %s\n", sizes->rPts, showCLInt(err));
        return err;
    }

    cm->rc = clCreateBuffer(ci->clctx, constBufFlags, sizes->rc, rc, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream r consts buffer of size %zu: %s\n", sizes->rc, showCLInt(err));
        return err;
    }

    cm->sg_dx = clCreateBuffer(ci->clctx, constBufFlags, sizes->sg_dx, sg.dx, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream sg_dx buffer of size %zu: %s\n", sizes->sg_dx, showCLInt(err));
        return err;
    }

    mwAlignedFree(r_pts);
    mwAlignedFree(rc);

    return CL_SUCCESS;
}

static cl_int createLBTrigBuffer(CLInfo* ci,
                                 SeparationCLMem* cm,
                                 const ASTRONOMY_PARAMETERS* ap,
                                 const INTEGRAL_AREA* ia,
                                 const SeparationSizes* sizes)
{
    cl_int err;
    LB_TRIG* lbts;
    const cl_mem_flags constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    lbts = precalculateLBTrig(ap, ia);
    cm->lbts = clCreateBuffer(ci->clctx, constBufFlags, sizes->lbts, lbts, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating lb_trig buffer of size %zu: %s\n", sizes->lbts, showCLInt(err));
        return err;
    }

    mwAlignedFree(lbts);

    return CL_SUCCESS;
}

void calculateSizes(SeparationSizes* sizes, const ASTRONOMY_PARAMETERS* ap, const INTEGRAL_AREA* ia)
{
    sizes->outMu = sizeof(real) * ia->mu_steps * ia->r_steps;
    sizes->outProbs = sizeof(real) * ia->mu_steps * ia->r_steps * ap->number_streams;
    sizes->ap = sizeof(ASTRONOMY_PARAMETERS);
    sizes->ia = sizeof(INTEGRAL_AREA);
    sizes->sc = sizeof(STREAM_CONSTANTS) * ap->number_streams;
    sizes->rPts = sizeof(R_POINTS) * ap->convolve * ia->r_steps;
    sizes->rc = sizeof(R_CONSTS) * ia->r_steps;
    sizes->sg_dx = sizeof(real) * ap->convolve;
    sizes->lbts = sizeof(LB_TRIG) * ia->mu_steps * ia->nu_steps;
}

cl_int createSeparationBuffers(CLInfo* ci,
                               SeparationCLMem* cm,
                               const ASTRONOMY_PARAMETERS* ap,
                               const INTEGRAL_AREA* ia,
                               const STREAM_CONSTANTS* sc,
                               const STREAM_GAUSS sg,
                               const SeparationSizes* sizes)
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
    err |= createRBuffers(ci, cm, ap, ia, sg, sizes);
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

