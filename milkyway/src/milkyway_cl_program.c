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

#include <errno.h>

#include "milkyway_util.h"
#include "milkyway_cl.h"
#include "milkyway_cl_program.h"

static char* mwGetBuildLog(CLInfo* ci)
{
    cl_int err;
    size_t logSize, readSize;
    char* buildLog;

    err = clGetProgramBuildInfo(ci->prog,
                                ci->dev,
                                CL_PROGRAM_BUILD_LOG,
                                0,
                                NULL,
                                &logSize);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to get build log size");
        return NULL;
    }

    buildLog = mwCalloc(sizeof(char), logSize + 1);

    err = clGetProgramBuildInfo(ci->prog,
                                ci->dev,
                                CL_PROGRAM_BUILD_LOG,
                                logSize,
                                buildLog,
                                &readSize);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to read program build log");
        free(buildLog);
        return NULL;
    }

    if (readSize != logSize)
        mw_printf("Failed to read complete build log\n");

    return buildLog;
}

static void CL_CALLBACK milkywayBuildCB(cl_program prog, void* user_data)
{
    cl_int infoErr;
    cl_build_status status;
    CLInfo* ci;
    char* buildLog;

    if (!user_data)
    {
        mw_printf("milkywayBuildCB got null user_data\n");
        return;
    }

    ci = (CLInfo*) user_data;

    assert(ci->prog == prog);

    infoErr = clGetProgramBuildInfo(prog,
                                    ci->dev,
                                    CL_PROGRAM_BUILD_STATUS,
                                    sizeof(status),
                                    &status,
                                    NULL);
    if (infoErr != CL_SUCCESS)
    {
        mwPerrorCL(infoErr, "Get build status failed");
    }
    else
    {
        mw_printf("Build status: %s\n", showCLBuildStatus(status));
    }

    buildLog = mwGetBuildLog(ci);

    mw_printf("Build log: \n%s\n", buildLog);
    free(buildLog);
}

/* Build program and create kernel */
cl_int mwBuildProgram(CLInfo* ci, const char* options)
{
    cl_int err = CL_SUCCESS;

    err = clBuildProgram(ci->prog, 1, &ci->dev, options, milkywayBuildCB, ci);
    if (err != CL_SUCCESS)
        mwPerrorCL(err, "clBuildProgram: Build failure");

    return err;
}

unsigned char* mwGetProgramBinary(CLInfo* ci, size_t* binSizeOut)
{
    cl_int err;
    size_t binSize = 0;
    unsigned char* bin = NULL;

    err = clGetProgramInfo(ci->prog, CL_PROGRAM_BINARY_SIZES, sizeof(binSize), &binSize, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to get program binary size");
        return NULL;
    }

    if (binSize == 0)
    {
        mw_printf("Program binary size == 0\n");
        return NULL;
    }

    bin = (unsigned char*) mwCalloc(binSize + 1, 1);
    err = clGetProgramInfo(ci->prog, CL_PROGRAM_BINARIES, sizeof(bin), &bin, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error getting program binary");
        free(bin);
        bin = NULL;
        binSize = 0;
    }

    if (binSizeOut)
        *binSizeOut = binSize;
    return bin;
}

int mwSaveProgramBinaryToFile(CLInfo* ci, const char* filename)
{
    size_t binSize;
    unsigned char* bin;
    FILE* f;
    int rc = 0;

    bin = mwGetProgramBinary(ci, &binSize);
    f = mw_fopen(filename, "wb");
    if (!f)
    {
        mwPerror("Error opening CL binary output '%s'", filename);
        return errno;
    }

    if (fwrite(bin, binSize, 1, f) != 1)
    {
        mwPerror("Error writing program binary to file '%s'\n", filename);
        rc = errno;
    }

    if (fclose(f) < 0)
    {
        mwPerror("Error closing program binary file '%s'", filename);
        rc = errno;
    }
    free(bin);

    return rc;
}



cl_int mwSetProgramFromBin(CLInfo* ci,const unsigned char* bin, size_t binSize)
{
    cl_int err;
    cl_int binStatus;

    ci->prog = clCreateProgramWithBinary(ci->clctx, 1, &ci->dev, &binSize, &bin, &binStatus, &err);
    mwPerrorCL(err, "Binary status");
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to create program from binary");
        return err;
    }

    if (binStatus != CL_SUCCESS)
    {
        mwPerrorCL(err, "Reading binary failed");
        return binStatus;
    }

    err = mwBuildProgram(ci, NULL);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error building program from binary");
        return err;
    }

    return CL_SUCCESS;
}

cl_int mwSetProgramFromSrc(CLInfo* ci,
                           cl_uint srcCount,
                           const char** src,
                           const size_t* lengths,
                           const char* compileDefs)
{
    cl_int err;

    ci->prog = clCreateProgramWithSource(ci->clctx, srcCount, src, lengths, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error creating program");
        return err;
    }

    err = mwBuildProgram(ci, compileDefs);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Error building program from source");
        return err;
    }

    clUnloadCompiler();

    return CL_SUCCESS;
}


cl_int mwCreateKernel(cl_kernel* kern, CLInfo* ci, const char* name)
{
    cl_int err = CL_SUCCESS;

    *kern = clCreateKernel(ci->prog, name, &err);
    if (err != CL_SUCCESS)
    {
        mwPerrorCL(err, "Failed to create kernel '%s'", name);
        return err;
    }

    return err;
}


