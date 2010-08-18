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
#include "separation.h"
#include "show_cl_types.h"
#include "separation_cl_buffers.h"

inline static cl_int createOutNuBuffer(const unsigned int r_steps,
                                       CLInfo* ci,
                                       SeparationCLMem* cm)
{
    cl_int err;
    cm->outNu = clCreateBuffer(ci->clctx,
                               CL_MEM_WRITE_ONLY,
                               sizeof(BG_PROB) * r_steps,
                               NULL,
                               &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out nu buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

inline static cl_int createSCBuffer(const STREAM_CONSTANTS* sc,
                                    const unsigned int number_streams,
                                    CLInfo* ci,
                                    SeparationCLMem* cm)
{
    cl_int err;
    size_t size = sizeof(STREAM_CONSTANTS) * number_streams;
    cm->sc = clCreateBuffer(ci->clctx,
                            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            size,
                            sc,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream constants buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

inline static cl_int createNuConstsBuffer(const NU_CONSTANTS* nu_consts,
                                          const unsigned int nu_steps,
                                          CLInfo* ci,
                                          SeparationCLMem* cm)
{
    cl_int err;
    size_t size = sizeof(NU_CONSTANTS) * nu_steps;
    cm->nuConsts = clCreateBuffer(ci->clctx,
                                  CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  size,
                                  nu_consts,
                                  &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating nu constants buffer of size %zu: %s\n", size, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

cl_int createSeparationBuffers(const ASTRONOMY_PARAMETERS* ap,
                               const INTEGRAL_AREA* ia,
                               const STREAM_CONSTANTS* sc,
                               const NU_CONSTANTS* nu_st,
                               CLInfo* ci,
                               SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;

    err |= createOutNuBuffer(ia->r_steps, ci, cm);
    err |= createSCBuffer(sc, ap->number_streams, ci, cm);
    err |= createNuConstsBuffer(nu_st, ia->nu_steps, ci, cm);

    return err;
}

void releaseSeparationBuffers(SeparationCLMem* cm)
{
    clReleaseMemObject(cm->outNu);
    clReleaseMemObject(cm->sc);
    clReleaseMemObject(cm->nuConsts);
}

