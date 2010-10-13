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
    return clCreateBuffer(clctx, CL_MEM_WRITE_ONLY, size, NULL, err);
}

static inline cl_int createOutMuBuffer(CLInfo* ci,
                                       SeparationCLMem* cm,
                                       const SeparationSizes* sizes)
{
    cl_int err;

    cm->outMu = createWriteBuffer(ci->clctx, sizes->outMu, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out mu buffer of size %zu: %s\n", sizes->outMu, showCLInt(err));
        return err;
    }

    cm->outMu_tmp = createWriteBuffer(ci->clctx, sizes->outMu, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out mu temp buffer of size %zu: %s\n", sizes->outMu, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createOutProbsBuffer(CLInfo* ci,
                                          SeparationCLMem* cm,
                                          const SeparationSizes* sizes)
{
    cl_int err;

    cm->outProbs = createWriteBuffer(ci->clctx, sizes->outProbs, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out probs buffer of size %zu: %s\n", sizes->outProbs, showCLInt(err));
        return err;
    }

    cm->outProbs_tmp = createWriteBuffer(ci->clctx, sizes->outProbs, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out probs temp buffer of size %zu: %s\n", sizes->outProbs, showCLInt(err));
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
                                    const STREAM_GAUSS* sg,
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

    mwAlignedFree(r_pts);
    mwAlignedFree(rc);

    return CL_SUCCESS;
}

static void calculateSizes(SeparationSizes* sizes, const ASTRONOMY_PARAMETERS* ap, const INTEGRAL_AREA* ia)
{
    sizes->outMu = sizeof(real) * ia->mu_steps * ia->r_steps;
    sizes->outProbs = sizeof(real) * ia->mu_steps * ia->r_steps * ap->number_streams;
    sizes->ap = sizeof(ASTRONOMY_PARAMETERS);
    sizes->ia = sizeof(INTEGRAL_AREA);
    sizes->sc = sizeof(STREAM_CONSTANTS) * ap->number_streams;
    sizes->rPts = sizeof(R_POINTS) * ap->convolve * ia->r_steps;
    sizes->rc = sizeof(R_CONSTS) * ia->r_steps;
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

    SeparationSizes sizes;

    calculateSizes(&sizes, ap, ia);

    if (ci->devType == CL_DEVICE_TYPE_CPU)
        constBufFlags = CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR;
    else
        constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    err |= createOutMuBuffer(ci, cm, &sizes);
    err |= createOutProbsBuffer(ci, cm, &sizes);
    err |= createAPBuffer(ci, cm, ap, &sizes, constBufFlags);
    err |= createIABuffer(ci, cm, ia, &sizes, constBufFlags);
    err |= createSCBuffer(ci, cm, sc, &sizes, constBufFlags);
    err |= createRBuffers(ci, cm, ap, ia, sg, &sizes);

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


void swapOutputBuffers(SeparationCLMem* cm)
{
    cl_mem tmp;

    tmp           = cm->outMu_tmp;
    cm->outMu_tmp = cm->outMu;
    cm->outMu     = tmp;

    tmp              = cm->outProbs_tmp;
    cm->outProbs_tmp = cm->outProbs;
    cm->outProbs     = tmp;
}

/* Set kernel arguments to the temporary output buffers */
cl_int separationSetOutputBuffers(CLInfo* ci, SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;

    /* Set output buffer arguments */
    err |= clSetKernelArg(ci->kern, 0, sizeof(cl_mem), &cm->outMu_tmp);
    err |= clSetKernelArg(ci->kern, 1, sizeof(cl_mem), &cm->outProbs_tmp);

    if (err != CL_SUCCESS)
        warn("Failed to set output buffer arguments: %s\n", showCLInt(err));

    return err;
}

