/*
  Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
  Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
  and Rensselaer Polytechnic Institute.

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

extern "C++" {
#include "../milkyway/parameters.h"
#include "../milkyway/star_points.h"
#include "../milkyway/evaluation_optimized.h"
#include "coords.h"
#include "cpu_coords.h"
#include "r_constants.h"
#include "pi_constants.h"
#include "gauss_legendre.h"
#include "boinc_api.h"
}

#include <cuda.h>
#include <cuda_runtime.h>
#include <cutil_inline.h>

#define MAX_CONVOLVE 120
#define R_INCREMENT 1
#define NUM_THREADS 64 //per block

#ifdef SINGLE_PRECISION
#define GPU_PRECISION float
#endif
#ifdef DOUBLE_PRECISION
#define GPU_PRECISION double
#endif

GPU_PRECISION** device_V;       //V -- float[nu][r]
GPU_PRECISION** device_sinb;
GPU_PRECISION** device_sinl;
GPU_PRECISION** device_cosb;
GPU_PRECISION** device_cosl;

GPU_PRECISION** host_bg_int;
GPU_PRECISION** host_st_int;
GPU_PRECISION** device_bg_int;
GPU_PRECISION** device_st_int;
GPU_PRECISION** device_bg_int_c;
GPU_PRECISION** device_st_int_c;
GPU_PRECISION** host_bg_int_c;
GPU_PRECISION** host_st_int_c;

__device__ __constant__ GPU_PRECISION constant_inverse_fstream_sigma_sq2[8];

__device__ __constant__ GPU_PRECISION constant_dx[MAX_CONVOLVE];
__device__ __constant__ GPU_PRECISION constant_qgaus_W[MAX_CONVOLVE];

__device__ __constant__ GPU_PRECISION constant_background_weight[1];
__device__ __constant__ GPU_PRECISION constant_exp_sum[1];
__device__ __constant__ GPU_PRECISION constant_stream_weight[4];

#define MAX_STREAMS (4)
__device__ __constant__ GPU_PRECISION constant_fstream_c[3*MAX_STREAMS];
__device__ __constant__ GPU_PRECISION constant_fstream_a[3*MAX_STREAMS];
#ifdef GLOBAL_MEMORY
GPU_PRECISION* device_fstream_c;
GPU_PRECISION* device_fstream_a;
#endif

GPU_PRECISION*   device_stars;

GPU_PRECISION*   device_probability;
GPU_PRECISION*   device_bg_only;
GPU_PRECISION*   device_st_only;
GPU_PRECISION*   host_probability;
GPU_PRECISION*   host_bg_only;
GPU_PRECISION*   host_st_only;

//int max_blocks = 128;
int max_blocks = 50000;

#define kernel3__mu_step    (mu_offset + blockIdx.x)
#define kernel3__mu_steps   (mu_steps)
//TODO check if blockIdx.y is correct
#define kernel3__r_step     (in_step + blockIdx.y)
#define kernel3__r_steps    in_steps
#define kernel3__nu_step    threadIdx.x
#define kernel3__nu_steps   blockDim.x

extern __shared__ GPU_PRECISION shared_mem[];

#ifdef SINGLE_PRECISION

#ifdef CONST_MEMORY
#include "evaluation_gpu6_float_const.cu"
#else
#ifdef GLOBAL_MEMORY
#ifdef SHARED_MEMORY
#include "evaluation_gpu6_float_shared.cu"
#else
#include "evaluation_gpu6_float_global.cu"
#endif
#else
#include "evaluation_gpu6_float.cu"
#endif
#endif


#else
#ifdef CONST_MEMORY
#include "evaluation_gpu6_double_const.cu"
#else
#ifdef GLOBAL_MEMORY
#ifdef SHARED_MEMORY
#include "evaluation_gpu6_double_shared.cu"
#else
#ifdef REGISTERS
#include "evaluation_gpu6_double_registers.cu"
#else
#include "evaluation_gpu6_double_global.cu"
#endif
#endif
#else
#include "evaluation_gpu6_double.cu"
#endif
#endif
#endif
#include "evaluation_gpu6_likelihood.cu"
#include "evaluation_gpu6_macros.cu"

void parse_prefs(char* project_prefs)
{
    printf("project_prefs: %s\n", project_prefs);
    if (project_prefs)
    {
        char* prefs_start = strstr(project_prefs, "<nvidia_block_amount>");
        char* prefs_end = strstr(project_prefs, "</nvidia_block_amount>");
        if (prefs_start && prefs_end)
        {
            int pref_length = prefs_end - (prefs_start + 21);
            printf("pref_length:%d\n", pref_length);
            char* pref_num = new char[pref_length+1];
            strncpy(pref_num, prefs_start + 21, pref_length);
            pref_num[pref_length] = '\0';
            printf("pref_num:%s\n", pref_num);
            max_blocks = atoi(pref_num);
            if (max_blocks < 0 || max_blocks > 50000)
                max_blocks = 128; //default
        }
    }
    fprintf(stderr, "Using %d concurrent blocks\n", max_blocks);
}

bool boinc_setup_gpu(int device)
{
    //from BOINC
    //http://boinc.berkeley.edu/trac/wiki/CudaApps
    CUdevice  hcuDevice;
    CUcontext hcuContext;
    CU_SAFE_CALL_NO_SYNC(cuDeviceGet( &hcuDevice, device));
    CU_SAFE_CALL_NO_SYNC(cuCtxCreate( &hcuContext, 0x4, hcuDevice ));
    return true;
}

int choose_cuda_13()
{
    //check for and find a CUDA 1.3 (double precision)
    //capable card
    int device_count;
    cutilSafeCall(cudaGetDeviceCount(&device_count));
    fprintf(stderr, "Found %d CUDA cards\n", device_count);
    if (device_count < 1)
    {
        fprintf(stderr, "No CUDA cards found, you cannot run the GPU version\n");
        exit(1);
    }
    int* eligable_devices = (int*)malloc(sizeof(int) * device_count);
    int eligable_device_idx = 0;
    int max_gflops = 0;
    int device;
    char* chosen_device = 0;
    for (int idx = 0; idx < device_count; ++idx)
    {
        cudaDeviceProp deviceProp;
        cutilSafeCall(cudaGetDeviceProperties(&deviceProp, idx));
        fprintf(stderr, "Found a %s\n", deviceProp.name);
        if ((deviceProp.major == 1 && deviceProp.minor >= 3) ||
                (deviceProp.major > 1))
        {
            eligable_devices[eligable_device_idx++] = idx;
            fprintf(stderr, "Device can be used it has compute capability >=1.3 support\n");
            //check how many gflops it has
            int gflops = deviceProp.multiProcessorCount * deviceProp.clockRate;
            if (gflops >= max_gflops)
            {
                max_gflops = gflops;
                device = idx;
                if (chosen_device)
                    free(chosen_device);
                chosen_device = (char*) malloc(sizeof(char) * strlen(deviceProp.name) + 1);
                strncpy(chosen_device, deviceProp.name, strlen(deviceProp.name));
                chosen_device[strlen(deviceProp.name)] = '\0';
            }
        }
        else
        {
            fprintf(stderr, "Device cannot be used, it does not have compute capability 1.3 support\n");
        }
    }
    free(eligable_devices);
    if (eligable_device_idx < 1)
    {
        fprintf(stderr, "No compute capability 1.3 cards have been found, exiting...\n");
        free(chosen_device);
        return -1;
    }
    else
    {
        fprintf(stderr, "Chose device %s\n", chosen_device);
        cutilSafeCall(cudaSetDevice(device));
        free(chosen_device);
        return device;
    }
}

int choose_gpu(int argc, char** argv)
{
    CUresult status = cuInit(0);
    if (status != CUDA_SUCCESS)
        return -1;
    //check for the --device command line argument
    int device_arg = -1; //invalid device
    for (unsigned int idx = 1; idx < argc; ++idx)
    {
        if (argv[idx])
        {
            if (strncmp("--device", argv[idx], 8) == 0)
            {
                //next arg is the number of the device
                if (idx + 1 < argc)
                {
                    device_arg = atoi(argv[idx+1]);
                    break;
                }
            }
        }
    }
    if (device_arg >= 0)
        fprintf(stderr, "Device index specified on the command line was %d\n", device_arg);
    else
        fprintf(stderr, "No device was specified on the command line\n");
#ifdef DOUBLE_PRECISION
    fprintf(stderr, "Looking for a Double Precision capable NVIDIA GPU\n");
    //check if the device from the command line has CUDA 1.3 support
    if (device_arg != -1)
    {
        int device_count;
        cutilSafeCall(cudaGetDeviceCount(&device_count));
        if (device_arg >= device_count)
        {
            fprintf(stderr, "Device Index on the command line is greater than the number of devices (%d devices %d specified)\n", device_count, device_arg);
            device_arg = choose_cuda_13();
        }
        else
        {
            cudaDeviceProp deviceProp;
            cutilSafeCall(cudaGetDeviceProperties(&deviceProp, device_arg));
            if ((deviceProp.major == 1 && deviceProp.minor >= 3) ||
                    (deviceProp.major > 1))
            {
                cutilSafeCall(cudaSetDevice(device_arg));
                fprintf(stderr, "The device %s specified on the command line can be used\n", deviceProp.name);
            }
            else
            {
                fprintf(stderr, "The device %s from the command line cannot be used because a device supporting compute capability 1.3 (Double Precision) is required\n",
                        deviceProp.name);
                device_arg = choose_cuda_13();
            }
        }
    }
    else
    {
        device_arg = choose_cuda_13();
    }
    if (device_arg == -1)
    {
        return -1;
    }
    else
    {
        if (!boinc_setup_gpu(device_arg))
        {
            fprintf(stderr, "Unable to setup CUDA sync context (will waste CPU cycles)\n");
        }
    }
    return 0;
#endif
#ifdef SINGLE_PRECISION
    cudaDeviceProp deviceProp;
    if (device_arg == -1)
    {
        //just choose the device with the most gflops
        device_arg  = cutGetMaxGflopsDeviceId();
        cutilSafeCall(cudaGetDeviceProperties(&deviceProp, device_arg));
        cutilSafeCall(cudaSetDevice(device_arg));
    }
    else
    {
        cutilSafeCall(cudaGetDeviceProperties(&deviceProp, device_arg));
        cutilSafeCall(cudaSetDevice(device_arg));
    }
    fprintf(stderr, "Using %s\n", deviceProp.name);
    if (!boinc_setup_gpu(device_arg))
    {
        fprintf(stderr, "Unable to setup cuda sync context (will waste CPU cycles)\n");
    }
    return 0;
#endif
}

/**
   desired_size is the number of threads needed
*/
int get_total_threads(int desired_size)
{
    //desired size needs to be divisible by NUM_THREADS
    if (desired_size % NUM_THREADS != 0)
    {
        return NUM_THREADS * ceil(desired_size /
                                  (double) NUM_THREADS);
    }
    else
    {
        return desired_size;
    }
}

double cuda_evaluator(double*, ASTRONOMY_PARAMETERS* ap, EVALUATION_STATE* es, STAR_POINTS* sp)
{
    int i;
    double likelihood;

    size_t num_integral_block = ap.number_integrals * sizeof(double);
    int* r_steps = (int*)malloc(ap.number_integrals * sizeof(int));
    int* mu_steps = (int*)malloc(ap.number_integrals * sizeof(int));
    int* nu_steps = (int*)malloc(ap.number_integrals * sizeof(int));
    double* r_min = (double*)malloc(num_integral_block);
    double* mu_min = (double*)malloc(num_integral_block);
    double* nu_min = (double*)malloc(num_integral_block);
    double* r_step_size = (double*)malloc(num_integral_block);
    double* mu_step_size = (double*)malloc(num_integral_block);
    double* nu_step_size = (double*)malloc(num_integral_block);
    for (i = 0; i < ap->number_integrals; ++i)
    {
        r_steps[i] = ap->integral[i]->r_steps;
        mu_steps[i] = ap->integral[i]->mu_steps;
        nu_steps[i] = ap->integral[i]->nu_steps;
        r_min[i] = ap->integral[i]->r_min;
        mu_min[i] = ap->integral[i]->mu_min;
        nu_min[i] = ap->integral[i]->nu_min;
        r_step_size[i] = ap->integral[i]->r_step_size;
        mu_step_size[i] = ap->integral[i]->mu_step_size;
        nu_step_size[i] = ap->integral[i]->nu_step_size;
    }

    gpu__initialize();

    likelihood = gpu__likelihood(double* parameters);

    gpu__free_constants(ap);

    free(r_steps);
    free(mu_steps);
    free(nu_steps);
    free(r_min);
    free(mu_min);
    free(nu_min);
    free(r_step_size);
    free(mu_step_size);
    free(nu_step_size);

    return likelihood;
}

void gpu__initialize(ASTRONOMY_PARAMETERS* ap)
{
    device_bg_int = new GPU_PRECISION*[ap->number_integrals];
    device_st_int = new GPU_PRECISION*[ap->number_integrals];
    host_bg_int = new GPU_PRECISION*[ap->number_integrals];
    host_st_int = new GPU_PRECISION*[ap->number_integrals];
    device_bg_int_c = new GPU_PRECISION*[ap->number_integrals];
    device_st_int_c = new GPU_PRECISION*[ap->number_integrals];
    host_bg_int_c = new GPU_PRECISION*[ap->number_integrals];
    host_st_int_c = new GPU_PRECISION*[ap->number_integrals];

    device_V = new GPU_PRECISION*[ap->number_integrals];
    device_sinb = new GPU_PRECISION*[ap->number_integrals];
    device_sinl = new GPU_PRECISION*[ap->number_integrals];
    device_cosb = new GPU_PRECISION*[ap->number_integrals];
    device_cosl = new GPU_PRECISION*[ap->number_integrals];

    printf("%d integrals %d streams\n", ap->number_integrals, ap->number_streams);

    allocate_cu_arrays(ap->number_integrals);
    for (int i = 0; i < ap->number_integrals; i++)
    {
        int int_size = ap->integral[i]->nu_steps *
                       ap->integral[i]->mu_steps;
        int gpu_int_size = get_total_threads(int_size);

        double* sinb = new double[gpu_int_size];
        double* sinl = new double[gpu_int_size];
        double* cosb = new double[gpu_int_size];
        double* cosl = new double[gpu_int_size];

        populate_lb(ap->sgr_coordinates, ap->wedge,
                    ap->integral[i]->mu_steps,
                    ap->integral[i]->mu_min,
                    ap->integral[i]->mu_step_size,
                    ap->integral[i]->nu_steps,
                    ap->integral[i]->nu_min,
                    ap->integral[i]->nu_step_size,
                    sinb, sinl, cosb, cosl);

        GPU_PRECISION* host_sinb = new GPU_PRECISION[gpu_int_size];
        GPU_PRECISION* host_sinl = new GPU_PRECISION[gpu_int_size];
        GPU_PRECISION* host_cosb = new GPU_PRECISION[gpu_int_size];
        GPU_PRECISION* host_cosl = new GPU_PRECISION[gpu_int_size];

        for (int j = 0; j < gpu_int_size; ++j)
        {
            host_sinb[j] = (GPU_PRECISION) sinb[j];
            host_sinl[j] = (GPU_PRECISION) sinl[j];
            host_cosb[j] = (GPU_PRECISION) cosb[j];
            host_cosl[j] = (GPU_PRECISION) cosl[j];
        }

        cutilSafeCall(cudaMalloc((void**) &(device_sinb[i]),
                                 gpu_int_size * sizeof(GPU_PRECISION) ));
        cutilSafeCall(cudaMalloc((void**) &(device_sinl[i]),
                                 gpu_int_size * sizeof(GPU_PRECISION)));
        cutilSafeCall(cudaMalloc((void**) &(device_cosb[i]),
                                 gpu_int_size * sizeof(GPU_PRECISION)));
        cutilSafeCall(cudaMalloc((void**) &(device_cosl[i]),
                                 gpu_int_size * sizeof(GPU_PRECISION)));
        cutilSafeCall(cudaMemcpy(device_sinb[i], host_sinb,
                                 gpu_int_size * sizeof(GPU_PRECISION),
                                 cudaMemcpyHostToDevice));
        cutilSafeCall(cudaMemcpy(device_sinl[i], host_sinl,
                                 gpu_int_size * sizeof(GPU_PRECISION),
                                 cudaMemcpyHostToDevice));
        cutilSafeCall(cudaMemcpy(device_cosb[i], host_cosb,
                                 gpu_int_size * sizeof(GPU_PRECISION),
                                 cudaMemcpyHostToDevice));
        cutilSafeCall(cudaMemcpy(device_cosl[i], host_cosl,
                                 gpu_int_size * sizeof(GPU_PRECISION),
                                 cudaMemcpyHostToDevice));

        delete [] host_sinb;
        delete [] host_sinl;
        delete [] host_cosb;
        delete [] host_cosl;

        double* irv = new double[ap->integral[i]->r_steps];
        double* reff_xr_rp3 = new double[ap->integral[i]->r_steps];
        double** qw_r3_N = new double*[ap->integral[i]->r_steps];
        double** r_point = new double*[ap->integral[i]->r_steps];
        double* ids = new double[ap->integral[i]->nu_steps];
        double* nus = new double[ap->integral[i]->nu_steps];
        double** r_in_mag = new double*[ap->integral[i]->r_steps];
        double** r_in_mag2 = new double*[ap->integral[i]->r_steps];
        GPU_PRECISION* host_V = new GPU_PRECISION[ap->integral[i]->nu_steps *
                ap->integral[i]->r_steps];

        cpu__r_constants(ap->convolve,
                         ap->integral[i]->r_steps,
                         ap->integral[i]->r_min,
                         ap->integral[i]->r_step_size,
                         ap->integral[i]->mu_steps,
                         ap->integral[i]->mu_min,
                         ap->integral[i]->mu_step_size,
                         ap->integral[i]->nu_steps,
                         ap->integral[i]->nu_min,
                         ap->integral[i]->nu_step_size,
                         irv, r_point,
                         r_in_mag, r_in_mag2,
                         qw_r3_N,
                         reff_xr_rp3, nus, ids);
        setup_r_point_texture(ap->integral[i]->r_steps,
                              ap->convolve, i, r_point);
        setup_qw_r3_N_texture(ap->integral[i]->r_steps,
                              ap->convolve, i, qw_r3_N);
        setup_r_in_mag_texture2(ap->integral[i]->r_steps,
                                ap->convolve, i, r_in_mag, r_in_mag2);


        for (int k = 0; k < ap->integral[i]->nu_steps; k++)
        {
            for (int j = 0; j < ap->integral[i]->r_steps; j++)
            {
                host_V[(j * ap->integral[i]->nu_steps) + k] =
                    reff_xr_rp3[j] * irv[j] * ids[k];
            }
        }
        delete [] irv;
        delete [] reff_xr_rp3;
        delete [] qw_r3_N;
        delete [] r_point;
        delete [] ids;
        delete [] nus;
        delete [] r_in_mag;
        delete [] r_in_mag2;

        cutilSafeCall(cudaMalloc((void**) &(device_V[i]),
                                 ap->integral[i]->nu_steps *
                                 ap->integral[i]->r_steps *
                                 sizeof(GPU_PRECISION)));
        cutilSafeCall(cudaMemcpy(device_V[i], host_V, ap->integral[i]->nu_steps *
                                 ap->integral[i]->r_steps * sizeof(GPU_PRECISION),
                                 cudaMemcpyHostToDevice) );

        cutilSafeCall(cudaMalloc((void**) &device_bg_int[i],
                                 gpu_int_size * sizeof(GPU_PRECISION)) );
        cutilSafeCall(cudaMalloc((void**) &device_st_int[i], ap->number_streams *
                                 gpu_int_size * sizeof(GPU_PRECISION)) );
        cutilSafeCall(cudaMalloc((void**) &device_bg_int_c[i], gpu_int_size *
                                 sizeof(GPU_PRECISION)) );
        cutilSafeCall(cudaMalloc((void**) &device_st_int_c[i],
                                 ap->number_streams *
                                 gpu_int_size * sizeof(GPU_PRECISION)) );
        host_bg_int[i] = new GPU_PRECISION[gpu_int_size];
        host_st_int[i] = new GPU_PRECISION[gpu_int_size * ap->number_streams];
        host_bg_int_c[i] = new GPU_PRECISION[gpu_int_size];
        host_st_int_c[i] = new GPU_PRECISION[gpu_int_size * ap->number_streams];
    }


    //printf("initializing constants for %d stars\n", number_stars);
    int gpu_num_stars = get_total_threads(sp->number_stars);
    GPU_PRECISION* host_stars = new GPU_PRECISION[gpu_num_stars * 5];
    for (int i = 0; i < sp->number_stars; i++)
    {
        host_stars[i + gpu_num_stars*0] = (GPU_PRECISION)sin(sp->stars[i][1]
                                          * D_DEG2RAD);
        host_stars[i + gpu_num_stars*1] = (GPU_PRECISION)sin(sp->stars[i][0]
                                          * D_DEG2RAD);
        host_stars[i + gpu_num_stars*2] = (GPU_PRECISION)cos(sp->stars[i][1]
                                          * D_DEG2RAD);
        host_stars[i + gpu_num_stars*3] = (GPU_PRECISION)cos(sp->stars[i][0]
                                          * D_DEG2RAD);
        host_stars[i + gpu_num_stars*4] = (GPU_PRECISION)sp->stars[i][2];
    }
    cutilSafeCall(cudaMalloc((void**) &device_stars, gpu_num_stars *
                             5 * sizeof(GPU_PRECISION)) );
    cutilSafeCall(cudaMemcpy(device_stars, host_stars, gpu_num_stars *
                             5 * sizeof(GPU_PRECISION), cudaMemcpyHostToDevice));
    delete [] host_stars;


    double* d_qgaus_W = new double[ap->convolve];
    double* d_qgaus_X = new double[ap->convolve];

    d_gauss_legendre(-1.0, 1.0, d_qgaus_X, d_qgaus_W, ap->convolve);
    GPU_PRECISION* host_dx = new GPU_PRECISION[ap->convolve];
    GPU_PRECISION* host_qgaus_W = new GPU_PRECISION[ap->convolve];
    for (int i = 0; i < ap->convolve; i++)
    {
        host_dx[i] = (GPU_PRECISION)(3.0 * d_stdev * d_qgaus_X[i]);
        host_qgaus_W[i] = (GPU_PRECISION)d_qgaus_W[i];
    }
    delete [] d_qgaus_W;
    delete [] d_qgaus_X;

    cutilSafeCall(cudaMemcpyToSymbol(constant_dx, host_dx, ap->convolve *
                                     sizeof(GPU_PRECISION), 0,
                                     cudaMemcpyHostToDevice) );
    cutilSafeCall(cudaMemcpyToSymbol(constant_qgaus_W, host_qgaus_W,
                                     ap->convolve * sizeof(GPU_PRECISION), 0,
                                     cudaMemcpyHostToDevice) );
    delete [] host_dx;
    delete [] host_qgaus_W;

    printf("mallocing %d probability\n", gpu_num_stars);
    cutilSafeCall(cudaMalloc((void**) &device_probability,
                             gpu_num_stars * sizeof(GPU_PRECISION)));
    cutilSafeCall(cudaMalloc((void**) &device_bg_only,
                             gpu_num_stars * sizeof(GPU_PRECISION)));
    cutilSafeCall(cudaMalloc((void**) &device_st_only,
                             ap->number_streams * gpu_num_stars * sizeof(GPU_PRECISION)));
    host_probability = new GPU_PRECISION[gpu_num_stars];
    host_bg_only = new GPU_PRECISION[gpu_num_stars];
    host_st_only = new GPU_PRECISION[gpu_num_stars * ap->number_streams];

    unsigned int total, free;
    cuMemGetInfo(&free, &total);
    fprintf(stderr, "Used %d/%d memory, %d remaining\n",
            (total - free) / 1024, total / 1024, free / 1024);
}

void gpu__free_constants(ASTRONOMY_PARAMETERS* ap)
{
    for (int i = 0; i < ap->number_integrals; i++)
    {
        cutilSafeCall(cudaFree(device_V[i]) );
        cutilSafeCall(cudaFree(device_bg_int[i]));
        cutilSafeCall(cudaFree(device_st_int[i]));
        cutilSafeCall(cudaFree(device_bg_int_c[i]));
        cutilSafeCall(cudaFree(device_st_int_c[i]));
        cutilSafeCall(cudaFree(device_sinb[i]));
        cutilSafeCall(cudaFree(device_sinl[i]));
        cutilSafeCall(cudaFree(device_cosb[i]));
        cutilSafeCall(cudaFree(device_cosl[i]));
        delete [] host_bg_int[i];
        delete [] host_st_int[i];
        delete [] host_bg_int_c[i];
        delete [] host_st_int_c[i];
    }

    delete [] host_bg_int;
    delete [] host_st_int;
    delete [] host_bg_int_c;
    delete [] host_st_int_c;
    delete [] device_V;
    delete [] device_sinb;
    delete [] device_sinl;
    delete [] device_cosb;
    delete [] device_cosl;
    delete [] device_bg_int;
    delete [] device_st_int;
    delete [] device_bg_int_c;
    delete [] device_st_int_c;

    cutilSafeCall(cudaFree(device_stars));
    cutilSafeCall(cudaFree(device_probability));
    cutilSafeCall(cudaFree(device_bg_only));
    cutilSafeCall(cudaFree(device_st_only));
    delete [] host_probability;
    delete [] host_bg_only;
    delete [] host_st_only;
}

template <unsigned int number_streams>
__global__ void
gpu__zero_integrals(int mu_offset,
                    int mu_steps,
                    int nu_steps,
                    int total_threads,
                    GPU_PRECISION* background_integrals,
                    GPU_PRECISION* stream_integrals)
{
    int pos = threadIdx.x + (blockDim.x * blockIdx.x) + mu_offset;
    background_integrals[pos] = 0;
    for (int i = 0; i < number_streams; i++)
    {
        stream_integrals[pos] = 0;
        pos += (total_threads);
    }
}

void cpu__sum_integrals(int iteration, double* background_integral,
                        double* stream_integrals)
{
    int int_size = ap->integral[iteration]->nu_steps *
                   ap->integral[iteration]->mu_steps;
    int gpu_int_size = get_total_threads(int_size);
    cutilSafeCall(cudaMemcpy(host_bg_int[iteration],
                             device_bg_int[iteration],
                             gpu_int_size * sizeof(GPU_PRECISION),
                             cudaMemcpyDeviceToHost) );
    cutilSafeCall(cudaMemcpy(host_bg_int_c[iteration],
                             device_bg_int_c[iteration],
                             gpu_int_size * sizeof(GPU_PRECISION),
                             cudaMemcpyDeviceToHost) );

    double sum = 0.0;
    double correction_sum = 0.0;
    for (int i = 0; i < int_size; i++)
    {
        sum += (double)(host_bg_int[iteration][i]);
        correction_sum += (host_bg_int_c[iteration][i]);
        //printf("background_integral[%d/%d]: %.15f\n", i,
        //int_size, host_bg_int[iteration][i]);q
        //printf("(sum(background_integral[%d/%d]: %.15lf\n", i,
        //integral_size[iteration], sum);
    }
    printf("Applying a bg correction of %.25f\n", correction_sum);
    sum += correction_sum;
    if (iteration == 0)
        *background_integral = sum;
    else
        *background_integral -= sum;

    cutilSafeCall(cudaMemcpy(host_st_int[iteration],
                             device_st_int[iteration],
                             ap->number_streams * gpu_int_size *
                             sizeof(GPU_PRECISION), cudaMemcpyDeviceToHost));
    cutilSafeCall(cudaMemcpy(host_st_int_c[iteration],
                             device_st_int_c[iteration],
                             ap->number_streams * gpu_int_size *
                             sizeof(GPU_PRECISION), cudaMemcpyDeviceToHost));
    for (int i = 0; i < ap->number_streams; i++)
    {
        sum = 0.0;
        double correction_sum = 0.0;
        for (int j = 0; j < int_size; j++)
        {
            sum += (double)(host_st_int[iteration]
                            [j + (i * gpu_int_size)]);
            correction_sum += (double)(host_st_int_c[iteration]
                                       [j + (i * gpu_int_size)]);
            // printf("[%d] stream_integral: %.15f\n", j,
            //     host_st_int[iteration][j + (i * gpu_int_size)]);
            // printf("sum: %.15lf\n", sum);
        }
        printf("Applyting a correction of %.25f\n", correction_sum);
        sum += correction_sum;
        if (iteration == 0)
            stream_integrals[i] = sum;
        else
            stream_integrals[i] -= sum;
    }
}

/********
 *  Likelihood calculation
 ********/
void cpu__sum_likelihood(int num_streams, double* probability,
                         double* bg_only, double* st_only,
                         int gpu_num_stars)
{
    cutilSafeCall(cudaMemcpy(host_probability, device_probability,
                             gpu_num_stars * sizeof(GPU_PRECISION),
                             cudaMemcpyDeviceToHost));
    cutilSafeCall(cudaMemcpy(host_bg_only, device_bg_only,
                             gpu_num_stars * sizeof(GPU_PRECISION),
                             cudaMemcpyDeviceToHost));
    cutilSafeCall(cudaMemcpy(host_st_only, device_st_only,
                             num_streams * gpu_num_stars * sizeof(GPU_PRECISION),
                             cudaMemcpyDeviceToHost));
    printf("number_stars:%d\n", sp->number_stars);
    *probability = 0.0;
    *bg_only = 0.0;
    double correction1 = 0.0;
    double correction2 = 0.0;
    double correction3[num_streams];
    for (int i = 0; i < num_streams; ++i)
    {
        correction3[i] = 0.0;
        st_only[i] = 0.0;
    }
    for (int i = 0; i < sp->number_stars; i++)
    {
        //printf("host_probability[%d]=%.15lf\n", i, host_probability[i]);
        //printf("host_bg[%d]=%.15lf\n", i, host_bg_only[i]);
        double temp = *probability;
        *probability += host_probability[i];
        correction1 += (host_probability[i] - (*probability - temp));
        temp = *bg_only;
        *bg_only += host_bg_only[i];
        correction2 += (host_bg_only[i] - (*bg_only - temp));
        for (int j = 0; j < num_streams; ++j)
        {
            temp = st_only[j];
            st_only[j] += host_st_only[j*gpu_num_stars+i];
            correction3[j] += (host_st_only[j*gpu_num_stars+i] - (st_only[j] - temp));
            //printf("host_st[%d,%d]=%.45lf correction=%.15lf\n",
            //     i, j, host_st_only[i+j * gpu_num_stars],
            //     correction3[j]);
        }
    }
    *probability += correction1;
    printf("Applyting a correction of %.25f for probability\n", correction1);
    *probability /= sp->number_stars;
    *bg_only += correction2;
    *bg_only /= sp->number_stars;
    printf("Applyting a correction of %.25f for bg_only\n", correction2);
    for (int i = 0; i < num_streams; ++i)
    {
        st_only[i] += correction3[i];
        st_only[i] /= sp->number_stars;
        printf("Applyting a correction of %.25f for st_only[%d]\n", correction3[i], i);
    }
}



/********
 *  Run the GPU kernels and get the probability
 ********/

#define stream_parameters(x, y) parameters[(x * 6) + y + 3]
#define stream_weights(x) parameters[(x * 6) + 2]
//#define background_weight parameters[0]
#define background_weight (0.0)
//#define alpha parameters[1]
#define q (parameters[0])
#define r0 (parameters[1])
//#define delta parameters[4]

double gpu__likelihood(double* parameters)
{
    double stream_c[3], lbr[3];
    GPU_PRECISION* fstream_a = new GPU_PRECISION[MAX_STREAMS * 3];
    GPU_PRECISION* fstream_c = new GPU_PRECISION[MAX_STREAMS * 3];
    GPU_PRECISION* fstream_sigma_sq2 = new GPU_PRECISION[MAX_STREAMS];
    for (int i = 0; i < MAX_STREAMS; ++i)
    {
        fstream_sigma_sq2[i] = 0.0;
        fstream_a[i*3 + 0] = 0.0;
        fstream_a[i*3 + 1] = 0.0;
        fstream_a[i*3 + 2] = 0.0;
        fstream_c[i*3 + 0] = 0.0;
        fstream_c[i*3 + 0] = 0.0;
        fstream_c[i*3 + 0] = 0.0;
    }

    for (int i = 0; i < ap->number_streams; i++)
    {
        fstream_sigma_sq2[i] = 1 / (GPU_PRECISION)(2.0 * stream_parameters(i, 4)
                               * stream_parameters(i, 4));
        //printf("fstream_sigma_sq2[%d] = %.15f\n", i, fstream_sigma_sq2[i]);
        fstream_a[(i * 3) + 0] = (GPU_PRECISION)(sin(stream_parameters(i, 2))
                                 * cos(stream_parameters(i, 3)));
        fstream_a[(i * 3) + 1] = (GPU_PRECISION)(sin(stream_parameters(i, 2))
                                 * sin(stream_parameters(i, 3)));
        fstream_a[(i * 3) + 2] = (GPU_PRECISION)(cos(stream_parameters(i, 2)));

        if (ap->sgr_coordinates == 0)
        {
            gc_eq_gal(ap->wedge, stream_parameters(i, 0) * D_DEG2RAD, 0 * D_DEG2RAD,
                      &(lbr[0]), &(lbr[1]));
        }
        else
        {
            gc_sgr_gal(ap->wedge, stream_parameters(i, 0) * D_DEG2RAD, 0 * D_DEG2RAD,
                       &(lbr[0]), &(lbr[1]));
        }
        lbr[2] = stream_parameters(i, 1);
        d_lbr2xyz(lbr, stream_c);

        fstream_c[(i * 3) + 0] = (GPU_PRECISION)stream_c[0];
        fstream_c[(i * 3) + 1] = (GPU_PRECISION)stream_c[1];
        fstream_c[(i * 3) + 2] = (GPU_PRECISION)stream_c[2];
        // printf("stream %d fstream_c[0] = %.15f\n", i, fstream_c[i*3]);
        // printf("stream %d fstream_c[1] = %.15f\n", i, fstream_c[i*3 + 1]);
        // printf("stream %d fstream_c[2] = %.15f\n", i, fstream_c[i*3 + 2]);

        //printf("stream %d fstream_a[0] = %.15f\n", i, fstream_a[i*3]);
        //printf("stream %d fstream_a[1] = %.15f\n", i, fstream_a[i*3 + 1]);
        //printf("stream %d fstream_a[2] = %.15f\n", i, fstream_a[i*3 + 2]);
    }
    cutilSafeCall(cudaMemcpyToSymbol(constant_fstream_c, fstream_c, 3 *
                                     MAX_STREAMS *
                                     sizeof(GPU_PRECISION), 0,
                                     cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpyToSymbol(constant_fstream_a, fstream_a, 3 *
                                     MAX_STREAMS *
                                     sizeof(GPU_PRECISION), 0,
                                     cudaMemcpyHostToDevice));
    setup_constant_textures(fstream_a, fstream_c,
                            fstream_sigma_sq2, MAX_STREAMS);
#ifdef GLOBAL_MEMORY
    cutilSafeCall(cudaMalloc((void**) &device_fstream_c, 3 *
                             ap->number_streams *
                             sizeof(GPU_PRECISION)) );
    cutilSafeCall(cudaMemcpy(device_fstream_c, fstream_c, 3 *
                             ap->number_streams *
                             sizeof(GPU_PRECISION),
                             cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMalloc((void**) &device_fstream_a, 3 *
                             ap->number_streams *
                             sizeof(GPU_PRECISION)) );
    cutilSafeCall(cudaMemcpy(device_fstream_a, fstream_a, 3 *
                             ap->number_streams *
                             sizeof(GPU_PRECISION),
                             cudaMemcpyHostToDevice));
#endif
    cutilSafeCall(cudaMemcpyToSymbol(constant_inverse_fstream_sigma_sq2,
                                     fstream_sigma_sq2,
                                     ap->number_streams * sizeof(GPU_PRECISION),
                                     0, cudaMemcpyHostToDevice));
    double background_integral = 0.0;
    double* stream_integrals = new double[ap->number_streams];
    for (int i = 0; i < ap->number_streams; i++)
        stream_integrals[i] = 0.0;

    double coeff = 1.0 / (d_stdev * sqrt(2.0 * D_PI));

    unsigned int total_work = 0;
    unsigned int current_work = 0;
    for (int  i = 0; i < ap->number_integrals; i++)
    {
        total_work += ap->integral[i]->r_steps;
    }
    unsigned int timer = 0;
    cutCreateTimer(&timer);
    cutResetTimer(timer);
    cutStartTimer(timer);

    double q_squared_inverse = 1 / (q * q);
    double bg_a, bg_b, bg_c;
    if (ap->aux_bg_profile == 0)
    {
        bg_a    = 0;
        bg_b    = 0;
        bg_c    = 0;
    }
    else if (ap->aux_bg_profile == 1)
    {
        bg_a    = ap->background_parameters[4]; //vickej2_bg
        bg_b    = ap->background_parameters[5]; //vickej2_bg
        bg_c    = ap->background_parameters[6]; //vickej2_bg
    }
    printf("bg_a:%f\n", bg_a);
    printf("bg_b:%f\n", bg_b);
    printf("bg_c:%f\n", bg_c);
    for (int i = 0; i < ap->number_integrals; i++)
    {
        bind_texture(i);
        int number_streams = ap->number_streams;
        int mu_steps = ap->integral[i]->mu_steps;
        int nu_steps = ap->integral[i]->nu_steps;
        int r_steps = ap->integral[i]->r_steps;
        int int_size = ap->integral[i]->mu_steps *
                       ap->integral[i]->nu_steps;
        int total_threads = get_total_threads(int_size);
        printf("int_size:%d, total_threads:%d\n", int_size, total_threads);
        int total_blocks = total_threads / NUM_THREADS;
        dim3 dimBlock(NUM_THREADS);
        printf("total_blocks:%d\n", total_blocks);
        for (int mu_step = 0; mu_step < total_blocks; mu_step += max_blocks)
        {
            int offset = mu_step * NUM_THREADS;
            dim3 dimGrid(min(max_blocks, total_blocks - mu_step), R_INCREMENT);
            EXECUTE_ZERO_INTEGRALS(device_bg_int[i],
                                   device_st_int[i]);
            EXECUTE_ZERO_INTEGRALS(device_bg_int_c[i],
                                   device_st_int_c[i]);
        }
#ifdef SHARED_MEMORY
        int shared_mem_size = (ap->number_streams + ap->number_streams * 3 * 2) *
                              sizeof(GPU_PRECISION);
#else
        int shared_mem_size = 0; //shared memory not used in unrolled loops
#endif
        printf("Using %d bytes of shared memory\n", shared_mem_size);
        for (int j = 0; j < r_steps; j += R_INCREMENT)
        {
            boinc_fraction_done(current_work / (double) total_work);
            for (int mu_step = 0; mu_step < total_blocks; mu_step += max_blocks)
            {
                dim3 dimGrid(min(max_blocks, total_blocks - mu_step), R_INCREMENT);
                if (ap->aux_bg_profile == 1)
                {
                    const int aux_bg_profile = 1;
                    EXECUTE_INTEGRAL_KERNEL;
                }
                else
                {
                    const int aux_bg_profile = 0;
                    EXECUTE_INTEGRAL_KERNEL;
                }
                cudaError err = cudaThreadSynchronize();
                if (err != cudaSuccess)
                {
                    fprintf(stderr, "Error executing gpu__integral_kernel3 error message: %s\n",
                            cudaGetErrorString(err));
                    exit(1);
                }
            }
            current_work += R_INCREMENT;
        }
        cpu__sum_integrals(i, &background_integral, stream_integrals);
    }
    fprintf(stderr, "<background_integral> %.20lf </background_integral>\n", background_integral);
    fprintf(stderr, "<stream_integrals>");
    for (int i = 0; i < ap->number_streams; ++i)
    {
        fprintf(stderr, " %.20lf", stream_integrals[i]);
    }
    fprintf(stderr, "</stream_integrals>\n");
    cutStopTimer(timer);
    float t = cutGetTimerValue(timer);
    fprintf(stderr, "gpu__integral_kernel3 took %f ms\n", t);
    double* stream_weight = new double[ap->number_streams];
    double exp_weight = exp(background_weight);
    double sum_exp_weights = exp_weight;
    double bg_weight = exp_weight / background_integral;
    for (int i = 0; i < ap->number_streams; i++)
    {
        exp_weight = exp(stream_weights(i));
        sum_exp_weights += exp_weight;
        stream_weight[i] = exp_weight / stream_integrals[i];
    }
    sum_exp_weights *= 0.001;
    GPU_PRECISION f_background_weight[1];
    GPU_PRECISION* f_stream_weight = new GPU_PRECISION[ap->number_streams];
    f_background_weight[0] = (GPU_PRECISION)( bg_weight / sum_exp_weights );
    //printf("bg_weight = %.15f\n", f_background_weight[0]);
    //printf("sum_exp_weights:%.15f\n", sum_exp_weights);
    for (int i = 0; i < ap->number_streams; i++)
    {
        f_stream_weight[i] = (GPU_PRECISION)(stream_weight[i] / sum_exp_weights);
        //printf("st_weight[%d] = %.15f\n", i, f_stream_weight[i]);
        //printf("f_stream_weight[%d] = %.15f\n", i, f_stream_weight[i]);
    }

    cutilSafeCall(cudaMemcpyToSymbol(constant_background_weight,
                                     f_background_weight,
                                     sizeof(GPU_PRECISION), 0,
                                     cudaMemcpyHostToDevice) );
    cutilSafeCall(cudaMemcpyToSymbol(constant_stream_weight,
                                     f_stream_weight,
                                     ap->number_streams * sizeof(GPU_PRECISION),
                                     0, cudaMemcpyHostToDevice) );

    double likelihood = 0.0;
    cutResetTimer(timer);
    cutStartTimer(timer);
    int gpu_num_stars = get_total_threads(sp->number_stars);
    int total_blocks = gpu_num_stars / NUM_THREADS;
    dim3 dimBlock(NUM_THREADS);
    //    printf("num streams:%u\n", number_streams);
    int number_stars = get_total_threads(sp->number_stars);
    for (int mu_step = 0; mu_step < total_blocks; mu_step += max_blocks)
    {
        dim3 dimGrid(min(max_blocks, total_blocks - mu_step), R_INCREMENT);
        int convolve = ap->convolve;
        int offset = mu_step * NUM_THREADS;
        //printf("mu_step:%d, offset:%d\n", mu_step, offset);
        if (ap->aux_bg_profile == 1)
        {
            const int aux_bg_profile = 1;
            EXECUTE_LIKELIHOOD_KERNEL;
        }
        else
        {
            const int aux_bg_profile = 0;
            EXECUTE_LIKELIHOOD_KERNEL;
        }
        cudaError_t err;
        err = cudaThreadSynchronize();
        if (err != cudaSuccess)
        {
            fprintf(stderr, "Error executing gpu__likelihood_kernel error message: %s\n",
                    cudaGetErrorString(err));
            exit(1);
        }
    }
    cutStopTimer(timer);
    t = cutGetTimerValue(timer);
    fprintf(stderr, "gpu__likelihood_kernel took %f ms\n", t);
    double bg_only;
    double* st_only = new double[ap->number_streams];
    cpu__sum_likelihood(ap->number_streams, &likelihood, &bg_only,
                        st_only, number_stars);
    fprintf(stderr, "<background_only_likelihood> %.20lf </background_only_likelihood>\n",
            bg_only - 3.0);
    fprintf(stderr, "<stream_only_likelihood>");
    for (int i = 0; i < ap->number_streams; i++)
    {
        fprintf(stderr, " %.20lf", st_only[i] - 3.0);
    }
    fprintf(stderr, " </stream_only_likelihood>\n");
    //printf("likelihood: %.30lf\n", likelihood);

    delete [] fstream_a;
    delete [] fstream_c;
    delete [] fstream_sigma_sq2;
    delete [] f_stream_weight;
    delete [] stream_integrals;
    delete [] stream_weight;
    cutDeleteTimer(timer);

    return likelihood;

}
