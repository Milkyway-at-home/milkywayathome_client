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

inline static cl_int createOutNuBuffer(const unsigned int nu_steps,
                                       SeparationCLInfo* ci,
                                       SeparationCLMem* cm)
{
    cl_int err;
    cm->ia = clCreateBuffer(ci->clctx,
                            CL_MEM_WRITE_ONLY,
                            sizeof(BG_PROB) * nu_steps,
                            NULL,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating out nu buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}


inline static cl_int createAPBuffer(const ASTRONOMY_PARAMETERS* ap,
                                    SeparationCLInfo* ci,
                                    SeparationCLMem* cm)
{
    cl_int err;
    cm->ap = clCreateBuffer(ci->clctx,
                            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            sizeof(ASTRONOMY_PARAMETERS),
                            ap,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating astronomy parameters buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

inline static cl_int createSCBuffer(const STREAM_CONSTANTS* sc,
                                    const unsigned int number_streams,
                                    SeparationCLInfo* ci,
                                    SeparationCLMem* cm)
{
    cl_int err;
    cm->sc = clCreateBuffer(ci->clctx,
                            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            sizeof(STREAM_CONSTANTS) * number_streams,
                            sc,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating stream constants buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

inline static cl_int createRConstsBuffer(const R_CONSTANTS* r_consts,
                                         const unsigned int r_steps,
                                         SeparationCLInfo* ci,
                                         SeparationCLMem* cm)
{
    cl_int err;
    cm->rConsts = clCreateBuffer(ci->clctx,
                                 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                 sizeof(R_CONSTANTS) * r_steps,
                                 r_consts,
                                 &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating r constants buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

inline static cl_int createRPointsBuffer(const R_POINTS* r_points,
                                         const unsigned int r_steps,
                                         const unsigned int convolve,
                                         SeparationCLInfo* ci,
                                         SeparationCLMem* cm)
{
    cl_int err;
    cm->rPoints = clCreateBuffer(ci->clctx,
                                 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                 sizeof(R_POINTS) * r_steps * convolve,
                                 r_points,
                                 &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating r points buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

inline static cl_int createNuConstsBuffer(const NU_CONSTANTS* nu_consts,
                                         const unsigned int nu_steps,
                                         SeparationCLInfo* ci,
                                         SeparationCLMem* cm)
{
    cl_int err;
    cm->rConsts = clCreateBuffer(ci->clctx,
                                 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                 sizeof(NU_CONSTANTS) * nu_steps,
                                 nu_consts,
                                 &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating nu constants buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

inline static cl_int createIABuffer(const INTEGRAL_AREA* ia,
                                    SeparationCLInfo* ci,
                                    SeparationCLMem* cm)
{
    cl_int err;
    cm->ia = clCreateBuffer(ci->clctx,
                            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            sizeof(INTEGRAL_AREA),
                            ia,
                            &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating integral area buffer: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

cl_int createSeparationBuffers(const ASTRONOMY_PARAMETERS* ap,
                               const STREAM_CONSTANTS* sc,
                               const R_CONSTANTS* r_consts,
                               const R_POINTS* r_points,
                               const NU_CONSTANTS* nu_st,
                               const INTEGRAL_AREA* ia,
                               SeparationCLInfo* ci,
                               SeparationCLMem* cm)
{
    cl_int err = CL_SUCCESS;

    err |= createOutNuBuffer(ia->nu_steps, ci, cm);
    err |= createAPBuffer(ap, ci, cm);
    err |= createSCBuffer(sc, ap->number_streams, ci, cm);
    err |= createRConstsBuffer(r_consts, ia->r_steps, ci, cm);
    err |= createRPointsBuffer(r_points, ia->r_steps, ap->convolve, ci, cm);
    err |= createNuConstsBuffer(nu_st, ia->nu_steps, ci, cm);
    err |= createIABuffer(ia, ci, cm);

    return err;
}

