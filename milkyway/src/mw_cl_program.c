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

#include <stdio.h>

#include "milkyway_cl.h"
#include "milkyway_util.h"
#include "mw_cl_show_types.h"
#include "mw_cl_program.h"

/* This doesn't seem to exist on OS X, but the callback on ATI on
 * Linux/Windows dies without it */
#ifndef CL_CALLBACK
  #define CL_CALLBACK
#endif

static char* mwGetBuildLog(CLInfo* ci)
{
    size_t logSize, readSize;
    char* buildLog;

    clGetProgramBuildInfo(ci->prog,
                          ci->dev,
                          CL_PROGRAM_BUILD_LOG,
                          0,
                          NULL,
                          &logSize);

    buildLog = callocSafe(sizeof(char), logSize + 1);

    clGetProgramBuildInfo(ci->prog,
                          ci->dev,
                          CL_PROGRAM_BUILD_LOG,
                          logSize,
                          buildLog,
                          &readSize);

    if (readSize != logSize)
        warn("Failed to read complete build log\n");

    return buildLog;
}

static void CL_CALLBACK milkywayBuildCB(cl_program prog, void* user_data)
{
    cl_int infoErr;
    cl_build_status stat;
    CLInfo* ci;
    char* buildLog;

    if (!user_data)
    {
        warn("milkywayBuildCB got null user_data\n");
        return;
    }

    ci = (CLInfo*) user_data;

    infoErr = clGetProgramBuildInfo(ci->prog,
                                    ci->dev,
                                    CL_PROGRAM_BUILD_STATUS,
                                    sizeof(stat),
                                    &stat,
                                    NULL);

    if (infoErr != CL_SUCCESS)
        warn("Get build status failed: %s\n", showCLInt(infoErr));
    else
        warn("Build status: %s\n", showCLBuildStatus(stat));

    buildLog = mwGetBuildLog(ci);

    warn("Build log: \n%s\n", buildLog);
    free(buildLog);
}

/* Build program and create kernel */
cl_int mwBuildProgram(CLInfo* ci, const char* options, const char* kernName)
{
    cl_int err = CL_SUCCESS;

    err = clBuildProgram(ci->prog, 1, &ci->dev, options, milkywayBuildCB, ci);
    if (err != CL_SUCCESS)
        warn("clBuildProgram: Build failure: %s\n", showCLInt(err));

    ci->kern = clCreateKernel(ci->prog, kernName, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating kernel '%s': %s\n", kernName, showCLInt(err));
        return err;
    }

    return err;
}

unsigned char* mwGetProgramBinary(CLInfo* ci, size_t* binSizeOut)
{
    cl_int err;
    size_t binSize;
    size_t readSize;
    unsigned char* bin = NULL;

    err = clGetProgramInfo(ci->prog, CL_PROGRAM_BINARY_SIZES, sizeof(binSize), &binSize, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to get program binary size: %s\n", showCLInt(err));
        return NULL;
    }

    bin = (unsigned char*) mallocSafe(binSize);
    err = clGetProgramInfo(ci->prog, CL_PROGRAM_BINARIES, binSize, &bin, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Error getting program binary: %s\n", showCLInt(err));
        free(bin);
        bin = NULL;
        binSize = 0;
    }

    if (binSizeOut)
        *binSizeOut = binSize;
    return bin;
}


cl_int mwSetProgramFromBin(CLInfo* ci, const char* kernName, const unsigned char* bin, size_t binSize)
{
    cl_int err;
    cl_int binStatus;

    ci->prog = clCreateProgramWithBinary(ci->clctx, 1, &ci->dev, &binSize, &bin, &binStatus, &err);
    warn("Binary status: %s\n", showCLInt(err));
    if (err != CL_SUCCESS)
    {
        warn("Failed to create program from binary: %s\n", showCLInt(err));
        return err;
    }

    if (binStatus != CL_SUCCESS)
    {
        warn("Reading binary failed: %s\n", showCLInt(err));
        return binStatus;
    }

    err = mwBuildProgram(ci, NULL, kernName);
    if (err != CL_SUCCESS)
    {
        warn("Error building program from binary: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

cl_int mwSetProgramFromSrc(CLInfo* ci,
                           const char* kernName,
                           const char** src,
                           const cl_uint srcCount,
                           const char* compileDefs)
{
    cl_int err;

    ci->prog = clCreateProgramWithSource(ci->clctx, srcCount, src, NULL, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating program: %s\n", showCLInt(err));
        return err;
    }

    err = mwBuildProgram(ci, compileDefs, kernName);
    if (err != CL_SUCCESS)
    {
        warn("Error building program from source: %s\n", showCLInt(err));
        return err;
    }

    ci->kern = clCreateKernel(ci->prog, kernName, &err);
    if (err != CL_SUCCESS)
    {
        warn("Error creating kernel '%s': %s\n", kernName, showCLInt(err));
        return err;
    }

    clUnloadCompiler();

    return CL_SUCCESS;
}


