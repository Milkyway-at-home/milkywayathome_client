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


#include <stdio.h>
#include <stdlib.h>

/* Windows doesn't have drand48, so do stupid things */
#define RANDOM_DOUBLE (rand() / RAND_MAX)


#if USE_FDLIBM
  #include <fdlibm.h>
#else
  #include <math.h>
#endif

#if TEST_OPENCL
  #ifdef __APPLE__
    #include <OpenCL/cl.h>
  #else
    #include <CL/cl.h>
  #endif /* __APPLE__ */
#endif /* TEST_OPENCL */

typedef struct
{
    double rnd1;
    double rnd2;

    double sqrtr;
    double cbrtr;
    double log1pr;
    double expm1r;
    double logr;
    double sinr;
    double cosr;
    double tanr;
    double powr;
    double sqrr;
    double expr;
    double invr;
    /* double fmar; */
} ResultSet;

#define BUFSIZE 10240

const char* structSrc =
"typedef struct\n"
"{\n"
"    double rnd1;\n"
"    double rnd2;\n"
"\n"
"    double sqrtr;\n"
"    double cbrtr;\n"
"    double log1pr;\n"
"    double expm1r;\n"
"    double logr;\n"
"    double sinr;\n"
"    double cosr;\n"
"    double tanr;\n"
"    double powr;\n"
"    double sqrr;\n"
"    double expr;\n"
"    double invr;\n"
    "} ResultSet;\n";


const char* precisionTestSrc =
"__kernel\n"
"void precisionTest(__global ResultSet* results,\n"
"                   __global double* randoms,\n"
"                   unsigned int idx)\n"
"{\n"
"    double rnd1, rnd2;\n"
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
"    res->powr   = pow(rnd1, (double)1.5);\n"
"    res->sqrr   = rnd1 * rnd1;\n"
"    res->invr   = 1.0 / rnd1;\n"
"}\n"
"";


void precisionTest(ResultSet* results,
                   double* randoms,
                   unsigned int idx)
{
    double rnd1, rnd2;
    ResultSet* res;

    rnd2 = randoms[idx];
    res  = &results[idx];

    rnd1 = rnd2;
    res->rnd1   = rnd1;
    res->rnd2   = rnd2;
    res->sqrtr  = sqrt(rnd1);
    res->cbrtr  = cbrt(rnd1);
    res->log1pr = log1p(rnd1);
    res->expm1r = expm1(rnd1);
    res->expr   = exp(rnd1);
    res->logr   = log(rnd1);
    res->sinr   = sin(rnd1);
    res->cosr   = cos(rnd1);
    res->tanr   = tan(rnd1);
    res->powr   = pow(rnd1, 1.5);
    res->sqrr   = rnd1 * rnd1;
    res->invr   = 1.0 / rnd1;
    /* res->fmar   = fma(rnd1, rnd1, rnd1); */
}

void printResult(FILE* f, ResultSet* r)
{
    fprintf(f,
            "\n--------------------\n"
            "rnd1  = %.30g\n"
            "rnd2  = %.30g\n"
            "sqrt  = %.30g\n"
            "cbrt  = %.30g\n"
            "expm1 = %.30g\n"
            "exp   = %.30g\n"
            "log   = %.30g\n"
            "log1p = %.30g\n"
            "sin   = %.30g\n"
            "cos   = %.30g\n"
            "tan   = %.30g\n"
            "pow   = %.30g\n"
            "sqr   = %.30g\n"
            "invr  = %.30g\n"
            /* "fmar  = %.30g\n" */
            "\n--------------------\n" ,
            r->rnd1,
            r->rnd2,
            r->sqrtr,
            r->cbrtr,
            r->expm1r,
            r->expr,
            r->logr,
            r->log1pr,
            r->sinr,
            r->cosr,
            r->tanr,
            r->powr,
            r->sqrr,
            r->invr
            /* r->fmar */

        );
}

double* fillRandoms(unsigned int n)
{
    unsigned int i;
    double* arr = malloc(sizeof(double) * n);

    for (i = 0; i < n; ++i)
        arr[i] = RANDOM_DOUBLE;
    /* arr[i] = drand48(); */

    return arr;
}

ResultSet* runTests(double* randoms, unsigned int n)
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
                             double* randoms,
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
                             sizeof(double) * n,
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
    double* randoms;
    ResultSet* results;

    /* srand48(seed); */
    srand(seed);

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
    cl_device_type device = CL_DEVICE_TYPE_GPU;
    long seed = 0;

    if (argc >= 2)
        n = strtod(argv[1], NULL);

    if (argc >= 3)
    {
        if ( strtol(argv[2], NULL, 10) == 0)
            device = CL_DEVICE_TYPE_CPU;
        else
            device = CL_DEVICE_TYPE_GPU;
    }

    if (argc >= 4)
        seed = (long) strtol(argv[3], NULL, 10);

    runPrecisionTest(device, seed, n);

    return 0;
}

