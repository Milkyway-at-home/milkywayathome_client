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

/* Quick test for comparing basic math library function output on
 * different platforms */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if TEST_OPENCL
  #ifdef __APPLE__
    #include <OpenCL/cl.h>
  #else
    #include <CL/cl.h>
  #endif /* __APPLE__ */
#endif /* TEST_OPENCL */

typedef struct
{
    float rnd1;
    float rnd2;

    float sqrtr;
    float cbrtr;
    float log1pr;
    float expm1r;
    float logr;
    float sinr;
    float cosr;
    float tanr;
    float powr;
    float sqrr;
    float expr;
    float invr;
    float fmar;
} ResultSet;

#define BUFSIZE 10240

const char* structSrc =
"typedef struct\n"
"{\n"
"    float rnd1;\n"
"    float rnd2;\n"
"\n"
"    float sqrtr;\n"
"    float cbrtr;\n"
"    float log1pr;\n"
"    float expm1r;\n"
"    float logr;\n"
"    float sinr;\n"
"    float cosr;\n"
"    float tanr;\n"
"    float powr;\n"
"    float sqrr;\n"
"    float expr;\n"
"    float invr;\n"
"    float fmar;\n"
    "} ResultSet;\n";


const char* precisionTestSrc =
"__kernel\n"
"void precisionTest(__global ResultSet* results,\n"
"                   __global float* randoms,\n"
"                   unsigned int idx)\n"
"{\n"
"    float rnd1, rnd2;\n"
"    __global ResultSet* res;\n"
"    const size_t id = get_global_id(0);\n"
"    rnd2 = randoms[id];\n"
"    res  = &results[id];\n"
"    rnd1 = rnd2;\n"
"    res->rnd1   = rnd1;\n"
"    res->rnd2   = rnd2;\n"
"    res->sqrtr  = sqrt(rnd1);\n"
"    res->cbrtr  = cbrt(rnd1);\n"
"    res->log1pr = log1p(rnd1);\n"
"    res->expm1r = expm1(rnd1);\n"
"    res->expr   = exp(rnd1);\n"
"    res->logr   = log(rnd1);\n"
"    res->sinr   = sin(rnd1);\n"
"    res->cosr   = cos(rnd1);\n"
"    res->tanr   = tan(rnd1);\n"
"    res->powr   = pow(rnd1, (float)1.5);\n"
"    res->sqrr   = rnd1 * rnd1;\n"
"    res->invr   = 1.0 / rnd1;\n"
"    res->fmar   = fma((float)rnd1, (float)rnd1, (float)rnd1);\n"
"}\n"
"";


void precisionTest(ResultSet* results,
                   float* randoms,
                   unsigned int idx)
{
    float rnd1, rnd2;
    ResultSet* res;

    rnd2 = randoms[idx];
    res  = &results[idx];

    rnd1 = rnd2;
    res->rnd1   = rnd1;
    res->rnd2   = rnd2;
    res->sqrtr  = sqrtf(rnd1);
    res->cbrtr  = cbrtf(rnd1);
    res->log1pr = log1pf(rnd1);
    res->expm1r = expm1f(rnd1);
    res->expr   = expf(rnd1);
    res->logr   = logf(rnd1);
    res->sinr   = sinf(rnd1);
    res->cosr   = cosf(rnd1);
    res->tanr   = tanf(rnd1);
    res->powr   = powf(rnd1, 1.5);
    res->sqrr   = rnd1 * rnd1;
    res->invr   = 1.0 / rnd1;
    res->fmar   = fma(rnd1, rnd1, rnd1);
}

void printResult(FILE* f, ResultSet* r)
{
    fprintf(f,
            "\n--------------------\n"
            "rnd1  = %g\n"
            "rnd2  = %g\n"
            "sqrt  = %g\n"
            "expm1 = %g\n"
            "exp   = %g\n"
            "log   = %g\n"
            "sin   = %g\n"
            "cos   = %g\n"
            "pow   = %g\n"
            "sqr   = %g\n"
            "\n--------------------\n" ,
            r->rnd1,
            r->rnd2,
            r->sqrtr,
            r->expm1r,
            r->expr,
            r->logr,
            r->sinr,
            r->cosr,
            r->powr,
            r->sqrr);
}

float* fillRandoms(unsigned int n)
{
    unsigned int i;
    float* arr = malloc(sizeof(float) * n);

    for (i = 0; i < n; ++i)
        arr[i] = drand48();

    return arr;
}

ResultSet* runTests(float* randoms, unsigned int n)
{
    unsigned int i;
    ResultSet* results = malloc(sizeof(ResultSet) * n);

    for (i = 0; i < n; ++i)
        precisionTest(results, randoms, i);
    return results;
}

void printResults(FILE* f, ResultSet* res, unsigned int n)
{
    unsigned int i;
    for (i = 0; i < n; ++i)
        printResult(f, &res[i]);
}

#if TEST_OPENCL

static ResultSet* runTestsCL(cl_device_type type,
                             float* randoms,
                             unsigned int n)
{
    cl_uint maxComputeUnits, clockFreq;
    cl_ulong memSize;
    char* compileDefinitions;
    unsigned int devCount;

    cl_device_id dev;
    cl_program prog;
    cl_kernel kern;

    cl_context ctx;
    cl_int err;
    cl_command_queue queue;

    ResultSet* res;

    err = clGetDeviceIDs(NULL, type, 1, &dev, &devCount);
    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error getting device: %d\n", err);
        return NULL;
    }

    if (devCount == 0)
    {
        fprintf(stderr, "Didn't find any devices\n");
        return NULL;
    }

    ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err);
    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error creating context: %d\n", err);
        return NULL;
    }

    queue = clCreateCommandQueue(ctx, dev, 0, &err);
    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error creating command Queue: %d\n", err);
        return NULL;
    }

    /* Print some device information */
    clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS,   sizeof(cl_uint),  &maxComputeUnits, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint),  &clockFreq, NULL);
    clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_SIZE,     sizeof(cl_ulong), &memSize, NULL);

    printf("arst device %d: %u %u %lu\n",
           (int) type, maxComputeUnits, clockFreq, (unsigned long) memSize);

    const char* allSrc[] = { structSrc, precisionTestSrc };

    prog = clCreateProgramWithSource(ctx, 2, allSrc, NULL, &err);
    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error creating program: %d\n", err);
        return NULL;
    }

    asprintf(&compileDefinitions, "-TEST_OPENCL=%d ", 0);

    err = clBuildProgram(prog,
                         1,
                         &dev,
                         NULL,
                         /* compileDefinitions */
                         NULL,
                         NULL);

    free(compileDefinitions);
    if (err != CL_SUCCESS)
    {
        char buildLog[BUFSIZE] = "";
        size_t failSize;

        clGetProgramBuildInfo(prog,
                              dev,
                              CL_PROGRAM_BUILD_LOG,
                              sizeof(buildLog),
                              buildLog,
                              &failSize);

        if (failSize > BUFSIZE)
        {
            char* bigBuf = calloc(sizeof(char), failSize + 1);

            clGetProgramBuildInfo(prog,
                                  dev,
                                  CL_PROGRAM_BUILD_LOG,
                                  failSize,
                                  bigBuf,
                                  NULL);

            printf("Large build message: \n%s\n", bigBuf);
            free(bigBuf);
        }

        fprintf(stderr, "Build failure: %d: log = %s\n", err, buildLog);
        return NULL;
    }

    kern = clCreateKernel(prog, "precisionTest", &err);
    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error creating kernel: %d\n", err);
        return NULL;
    }

    clUnloadCompiler();


    cl_mem resultBuf;
    cl_mem randBuf;


    randBuf = clCreateBuffer(ctx,
                             CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             sizeof(float) * n,
                             randoms,
                             &err);
    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error creating rand buffer: %d\n", err);
        return NULL;
    }

    resultBuf = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, sizeof(ResultSet) * n, NULL, &err);
    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error creating result buffer: %d\n", err);
        return NULL;
    }


    err |= clSetKernelArg(kern, 0, sizeof(cl_mem), &resultBuf);
    err |= clSetKernelArg(kern, 1, sizeof(cl_mem), &randBuf);
    err |= clSetKernelArg(kern, 2, sizeof(unsigned int), &n);
    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error setting kernel arguments: %d\n", err);
        return NULL;
    }

    const size_t global[] = { n };
    err = clEnqueueNDRangeKernel(queue,
                                 kern,
                                 1,
                                 NULL, global, NULL,
                                 0, NULL, NULL);
    if (err)
    {
        fprintf(stderr, "Failed to enqueue kernel\n");
        return NULL;
    }



    res = malloc(sizeof(ResultSet) * n);
    err = clEnqueueReadBuffer(queue,
                              resultBuf,
                              CL_TRUE,
                              0, sizeof(ResultSet) * n, res,
                              0, NULL, NULL);

    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error reading CL result buffer\n");
        free(res);
        clReleaseMemObject(resultBuf);
        clReleaseMemObject(randBuf);
        clReleaseCommandQueue(queue);
        clReleaseProgram(prog) ;
        clReleaseKernel(kern);
        clReleaseContext(ctx);
        return NULL;
    }

    clReleaseMemObject(resultBuf);
    clReleaseMemObject(randBuf);
    clReleaseCommandQueue(queue);
    clReleaseProgram(prog) ;
    clReleaseKernel(kern);
    clReleaseContext(ctx);

    return res;
}

#endif /* TEST_OPENCL */

#if !TEST_OPENCL
  #define cl_device_type int
  #define CL_DEVICE_TYPE_CPU 0
  #define CL_DEVICE_TYPE_GPU 1
#endif

void runPrecisionTest(cl_device_type device, const long seed, const unsigned int n)
{
    float* randoms;
    ResultSet* results;

    srand48(seed);

    randoms = fillRandoms(n);

  #if TEST_OPENCL
    printf("Running OpenCL test\n");
    results = runTestsCL(device, randoms, n);
  #else
    printf("Running normal test\n");
    results = runTests(randoms, n);
  #endif

    if (!results)
    {
        free(randoms);
        fprintf(stderr, "Failed to get results\n");
        return;
    }

    printResults(stdout, results, n);

    free(randoms);
    free(results);
}

int main(int argc, char** argv)
{
    unsigned int n = 10000;
    cl_device_type device = CL_DEVICE_TYPE_CPU;
    long seed = 0;

    if (argc >= 3)
    {
        n = strtod(argv[1], NULL);
        device = (cl_device_type) strtol(argv[2], NULL, 10);
    }

    if (argc >= 4)
        seed = (long) strtol(argv[3], NULL, 10);

    runPrecisionTest(device, seed, n);

    return 0;
}

