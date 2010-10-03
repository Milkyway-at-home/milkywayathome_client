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
#include "show_cl_types.h"
#include "separation_cl_buffers.h"
#include "r_points.h"

static inline cl_mem createWriteBuffer(cl_context clctx, size_t size, cl_int* err)
{
    return  clCreateBuffer(clctx, CL_MEM_WRITE_ONLY, size, NULL, err);
}

static inline cl_int createOutMuBuffer(const unsigned int mu_steps,
                                       const unsigned int r_steps,
                                       CLInfo* ci,
                                       SeparationCLMem* cm)
{
    cl_int err;
    size_t size = sizeof(real) * mu_steps * r_steps;

    cm->outMu = createWriteBuffer(ci->clctx, size, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out mu buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    cm->outMu_tmp = createWriteBuffer(ci->clctx, size, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out mu temp buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createOutProbsBuffer(const unsigned int mu_steps,
                                          const unsigned int r_steps,
                                          const unsigned int number_streams,
                                          CLInfo* ci,
                                          SeparationCLMem* cm)
{
    cl_int err;
    size_t size = sizeof(real) * mu_steps * r_steps * number_streams;

    cm->outProbs = createWriteBuffer(ci->clctx, size, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out probs buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    cm->outProbs_tmp = createWriteBuffer(ci->clctx, size, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out probs temp buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createSCBuffer(const STREAM_CONSTANTS* sc,
                                    const unsigned int number_streams,
                                    CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const cl_mem_flags constBufFlags)
{
    cl_int err;
    size_t size = sizeof(STREAM_CONSTANTS) * number_streams;

    cm->sc = clCreateBuffer(ci->clctx,
                            constBufFlags,
                            size,
                            (void*) sc,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream constants buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createAPBuffer(const ASTRONOMY_PARAMETERS* ap,
                                    CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const cl_mem_flags constBufFlags)
{
    cl_int err;
    cm->ap = clCreateBuffer(ci->clctx,
                            constBufFlags,
                            sizeof(ASTRONOMY_PARAMETERS),
                            (void*) ap,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating astronomy parameters buffer of size %zu: %s\n", sizeof(ASTRONOMY_PARAMETERS), showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createIABuffer(const INTEGRAL_AREA* ia,
                                    CLInfo* ci,
                                    SeparationCLMem* cm,
                                    const cl_mem_flags constBufFlags)
{
    cl_int err;
    cm->ia = clCreateBuffer(ci->clctx, constBufFlags, sizeof(INTEGRAL_AREA), (void*) ia, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating integral area buffer of size %zu: %s\n", sizeof(INTEGRAL_AREA), showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createRBuffers(const ASTRONOMY_PARAMETERS* ap,
                                    const INTEGRAL_AREA* ia,
                                    const STREAM_GAUSS* sg,
                                    CLInfo* ci,
                                    SeparationCLMem* cm)
{
    cl_int err;

    R_POINTS* r_pts;
    R_CONSTS* rc;
    size_t rPtsSize = sizeof(R_POINTS) * ap->convolve * ia->r_steps;
    size_t rcSize = sizeof(R_CONSTS) * ia->r_steps;

    r_pts = precalculate_r_pts(ap, ia, sg, &rc);

    cm->rPts = clCreateBuffer(ci->clctx,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                              rPtsSize,
                              r_pts,
                              &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream r points buffer of size %zu: %s\n", rPtsSize, showCLInt(err));
        return err;
    }

    cm->rc = clCreateBuffer(ci->clctx,
                            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            rPtsSize,
                            rc,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream r consts buffer of size %zu: %s\n", rcSize, showCLInt(err));
        return err;
    }

    mwAlignedFree(r_pts);
    mwAlignedFree(rc);

    return CL_SUCCESS;
}


cl_int createSeparationBuffers(const ASTRONOMY_PARAMETERS* ap,
                               const INTEGRAL_AREA* ia,
                               const STREAM_CONSTANTS* sc,
                               const STREAM_GAUSS* sg,
                               CLInfo* ci,
                               SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;
    cl_mem_flags constBufFlags;

    if (ci->devType == CL_DEVICE_TYPE_CPU)
        constBufFlags = CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR;
    else
        constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    err |= createOutMuBuffer(ia->mu_steps, ia->r_steps, ci, cm);
    err |= createOutProbsBuffer(ia->mu_steps, ia->r_steps, ap->number_streams, ci, cm);
    err |= createAPBuffer(ap, ci, cm, constBufFlags);
    err |= createIABuffer(ia, ci, cm, constBufFlags);
    err |= createSCBuffer(sc, ap->number_streams, ci, cm, constBufFlags);
    err |= createRBuffers(ap, ia, sg, ci, cm);

    return err;
}

void releaseSeparationBuffers(SeparationCLMem* cm)
{
    clReleaseMemObject(cm->outProbs);
    clReleaseMemObject(cm->outProbs_tmp);

    clReleaseMemObject(cm->outMu);
    clReleaseMemObject(cm->outMu_tmp);

    clReleaseMemObject(cm->ap);
    clReleaseMemObject(cm->ia);
    clReleaseMemObject(cm->sc);
    clReleaseMemObject(cm->rPts);
    clReleaseMemObject(cm->rc);
}


static inline void swapOutputBuffers(SeparationCLMem* cm)
{
    cl_mem tmp;

    tmp           = cm->outMu_tmp;
    cm->outMu_tmp = cm->outMu;
    cm->outMu     = tmp;

    tmp              = cm->outProbs_tmp;
    cm->outProbs_tmp = cm->outProbs;
    cm->outProbs     = tmp;
}

cl_int separationSwapOutputBuffers(const CLInfo* ci, SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;

    printf("Swapping buffers\n");

    swapOutputBuffers(cm);

    /* Output buffers */
    err |= clSetKernelArg(ci->kern, 0, sizeof(cl_mem), &cm->outMu_tmp);
    err |= clSetKernelArg(ci->kern, 1, sizeof(cl_mem), &cm->outProbs_tmp);

    if (err != CL_SUCCESS)
        warn("Failed to set output buffer arguments: %s\n", showCLInt(err));

    return err;
}

