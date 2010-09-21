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

static inline cl_int createOutMuBuffer(const unsigned int r_steps,
                                       const unsigned int nu_steps,
                                       CLInfo* ci,
                                       SeparationCLMem* cm)
{
    cl_int err;
    size_t size = sizeof(BG_PROB) * r_steps * nu_steps;
    cm->outMu = clCreateBuffer(ci->clctx,
                               CL_MEM_WRITE_ONLY,
                               size,
                               NULL,
                               &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out mu buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline cl_int createOutProbsBuffer(const unsigned int r_steps,
                                          const unsigned int nu_steps,
                                          const unsigned int number_streams,
                                          CLInfo* ci,
                                          SeparationCLMem* cm)
{
    cl_int err;
    cm->outProbs = clCreateBuffer(ci->clctx,
                                  CL_MEM_WRITE_ONLY,
                                  sizeof(ST_PROBS) * r_steps * nu_steps * number_streams,
                                  NULL,
                                  &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out probs buffer: %s\n", showCLInt(err));
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

static inline cl_int createRPtsBuffer(const R_POINTS* r_pts_all,
                                      const unsigned int nconvolve,
                                      const unsigned int r_steps,
                                      CLInfo* ci,
                                      SeparationCLMem* cm,
                                      const cl_mem_flags constBufFlags)
{
    cl_int err;
    size_t size = sizeof(R_POINTS) * nconvolve * r_steps;
    cm->rPts = clCreateBuffer(ci->clctx,
                              constBufFlags,
                              size,
                              (void*) r_pts_all,
                              &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating r_pts buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

cl_int createSeparationBuffers(const ASTRONOMY_PARAMETERS* ap,
                               const INTEGRAL_AREA* ia,
                               const STREAM_CONSTANTS* sc,
                               const R_POINTS* r_pts_all,
                               CLInfo* ci,
                               SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;
    cl_mem_flags constBufFlags;

    if (ci->devType == CL_DEVICE_TYPE_CPU)
        constBufFlags = CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR;
    else
        constBufFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    err |= createOutMuBuffer(ia->r_steps, ia->nu_steps, ci, cm);
    err |= createOutProbsBuffer(ia->r_steps, ia->nu_steps, ap->number_streams, ci, cm);
    err |= createAPBuffer(ap, ci, cm, constBufFlags);
    err |= createIABuffer(ia, ci, cm, constBufFlags);
    err |= createSCBuffer(sc, ap->number_streams, ci, cm, constBufFlags);
    err |= createRPtsBuffer(r_pts_all, ap->convolve, ia->r_steps, ci, cm, constBufFlags);

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
}

