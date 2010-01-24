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
#include "../astronomy/parameters.h"
#include "../astronomy/star_points.h"
#include "../astronomy/evaluation_optimized.h"
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
//#define NUM_BLOCKS 128
#define NUM_BLOCKS 25000
#define NUM_THREADS 64 //per block

#ifdef SINGLE_PRECISION
#define GPU_PRECISION float
#endif
#ifdef DOUBLE_PRECISION
#define GPU_PRECISION double
#endif

extern ASTRONOMY_PARAMETERS *ap;
extern STAR_POINTS *sp;

GPU_PRECISION **device_V;		//V -- float[nu][r]
GPU_PRECISION **device_sinb;
GPU_PRECISION **device_sinl;
GPU_PRECISION **device_cosb;
GPU_PRECISION **device_cosl;

GPU_PRECISION **host_bg_int;
GPU_PRECISION **host_st_int;
GPU_PRECISION **device_bg_int;
GPU_PRECISION **device_st_int;

#ifdef SINGLE_PRECISION
GPU_PRECISION **device_bg_correction;
GPU_PRECISION **device_st_correction;
#endif

__device__ __constant__ GPU_PRECISION constant_inverse_fstream_sigma_sq2[8];

__device__ __constant__ GPU_PRECISION constant_dx[MAX_CONVOLVE];
__device__ __constant__ GPU_PRECISION constant_qgaus_W[MAX_CONVOLVE];

__device__ __constant__ GPU_PRECISION constant_background_weight[1];
__device__ __constant__ GPU_PRECISION constant_stream_weight[4];


GPU_PRECISION	*device_stars;

GPU_PRECISION	*device_probability;
GPU_PRECISION	*device_probability_correction;
GPU_PRECISION	*host_probability;

#define kernel3__mu_step	(mu_offset + blockIdx.x)
#define kernel3__mu_steps	(mu_steps)
//TODO check if blockIdx.y is correct
#define kernel3__r_step		(in_step + blockIdx.y)
#define kernel3__r_steps	in_steps
#define kernel3__nu_step	threadIdx.x
#define kernel3__nu_steps	blockDim.x

extern __shared__ GPU_PRECISION shared_mem[];

#ifdef SINGLE_PRECISION
#include "evaluation_gpu6_float.cu"
#else
#include "evaluation_gpu6_double.cu"
#endif
#include "evaluation_gpu6_likelihood.cu"
#include "evaluation_gpu6_macros.cu"

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
  int *eligable_devices = (int*)malloc(sizeof(int) * device_count);
  int eligable_device_idx = 0;
  int max_gflops = 0;
  int device;
  char *chosen_device = 0;
  for(int idx = 0;idx<device_count;++idx)
    {
      cudaDeviceProp deviceProp;
      cutilSafeCall(cudaGetDeviceProperties(&deviceProp, idx));
      fprintf(stderr, "Found a %s\n", deviceProp.name);
      if (deviceProp.major == 1 && deviceProp.minor == 3)
	{
	  eligable_devices[eligable_device_idx++] = idx;
	  fprintf(stderr, "Device can be used it has compute capability 1.3 support\n");
	  //check how many gflops it has
	  int gflops = deviceProp.multiProcessorCount * deviceProp.clockRate;
	  if (gflops >= max_gflops)
	    {
	      max_gflops = gflops;
	      device = idx;
	      if (chosen_device)
		free(chosen_device);
	      chosen_device = (char*) malloc(sizeof(char) * strlen(deviceProp.name)+1);
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
  if (eligable_device_idx < 1) {
    fprintf(stderr, "No compute capability 1.3 cards have been found, exiting...\n");
    free(chosen_device);
    return -1;
  } else {
    fprintf(stderr, "Chose device %s\n", chosen_device);
    cutilSafeCall(cudaSetDevice(device));  
    free(chosen_device);
    return device;
  }
}

int choose_gpu(int argc, char **argv) {
  CUresult status = cuInit(0);
  if(status != CUDA_SUCCESS)
    return -1;
  //check for the --device command line argument
  int device_arg = -1; //invalid device
  for(unsigned int idx = 1;idx < argc;++idx) {
    if (strncmp("--device", argv[idx], 8) == 0) {
      //next arg is the number of the device
      if (idx+1 < argc)
	{
	  device_arg = atoi(argv[idx+1]);
	  break;
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
  if (device_arg != -1) {
    cudaDeviceProp deviceProp;
    cutilSafeCall(cudaGetDeviceProperties(&deviceProp, device_arg));
    if (deviceProp.major == 1 && deviceProp.minor == 3) {
      cutilSafeCall(cudaSetDevice(device_arg));
      fprintf(stderr, "The device %s specified on the command line can be used\n", deviceProp.name);
    } else {
      fprintf(stderr, "The device %s from the command line cannot be used because a device supporting compute capability 1.3 (Double Precision) is required\n",
	      deviceProp.name);
      device_arg = choose_cuda_13();
    }
  } else {
    device_arg = choose_cuda_13();
  }
  if (device_arg == -1) {
    return -1;
  } else {
    if (!boinc_setup_gpu(device_arg))
      {
	fprintf(stderr, "Unable to setup CUDA sync context (will waste CPU cycles)\n");
      }
  }
  return 0;
#endif
#ifdef SINGLE_PRECISION
  cudaDeviceProp deviceProp;
  if (device_arg == -1) {
    //just choose the device with the most gflops
    device_arg  = cutGetMaxGflopsDeviceId();
    cutilSafeCall(cudaGetDeviceProperties(&deviceProp, device_arg));
    cutilSafeCall(cudaSetDevice(device_arg));
  } else {
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

void gpu__initialize()
{ 
  device_bg_int = new GPU_PRECISION*[ap->number_integrals];
  device_st_int = new GPU_PRECISION*[ap->number_integrals];
  host_bg_int = new GPU_PRECISION*[ap->number_integrals];
  host_st_int = new GPU_PRECISION*[ap->number_integrals];
#ifdef SINGLE_PRECISION
  device_bg_correction = new GPU_PRECISION*[ap->number_integrals];
  device_st_correction = new GPU_PRECISION*[ap->number_integrals];
#endif

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

      GPU_PRECISION *host_sinb = new GPU_PRECISION[gpu_int_size];
      GPU_PRECISION *host_sinl = new GPU_PRECISION[gpu_int_size];
      GPU_PRECISION *host_cosb = new GPU_PRECISION[gpu_int_size];
      GPU_PRECISION *host_cosl = new GPU_PRECISION[gpu_int_size];

      populate_lb(ap->sgr_coordinates, ap->wedge, 
		  ap->integral[i]->mu_steps,
		  ap->integral[i]->mu_min,
		  ap->integral[i]->mu_step_size, 
		  ap->integral[i]->nu_steps,
		  ap->integral[i]->nu_min,
		  ap->integral[i]->nu_step_size,
		  host_sinb, host_sinl, host_cosb, host_cosl);

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

      GPU_PRECISION *irv = new GPU_PRECISION[ap->integral[i]->r_steps];
      GPU_PRECISION *reff_xr_rp3 = new GPU_PRECISION[ap->integral[i]->r_steps];
      GPU_PRECISION **qw_r3_N = new GPU_PRECISION*[ap->integral[i]->r_steps];
      GPU_PRECISION **r_point = new GPU_PRECISION*[ap->integral[i]->r_steps];
      GPU_PRECISION *ids = new GPU_PRECISION[ap->integral[i]->nu_steps];
      GPU_PRECISION *nus = new GPU_PRECISION[ap->integral[i]->nu_steps];
      GPU_PRECISION *host_V = new GPU_PRECISION[ap->integral[i]->nu_steps * 
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
		       irv, r_point, qw_r3_N, 
		       reff_xr_rp3, nus, ids);
      setup_r_point_texture(ap->integral[i]->r_steps, 
			    ap->convolve, i, r_point);
      setup_qw_r3_N_texture(ap->integral[i]->r_steps, 
			    ap->convolve, i, qw_r3_N);

      for (int k = 0; k < ap->integral[i]->nu_steps; k++) {
	for (int j = 0; j < ap->integral[i]->r_steps; j++) {
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
#ifdef SINGLE_PRECISION
      cutilSafeCall(cudaMalloc((void**) &device_bg_correction[i], gpu_int_size *
			       sizeof(GPU_PRECISION)) );
      cutilSafeCall(cudaMalloc((void**) &device_st_correction[i], 
			       ap->number_streams *
			       gpu_int_size * sizeof(GPU_PRECISION)) );
#endif
      host_bg_int[i] = new GPU_PRECISION[gpu_int_size];
      host_st_int[i] = new GPU_PRECISION[gpu_int_size * ap->number_streams];
    }
	

  //printf("initializing constants for %d stars\n", number_stars);
  int gpu_num_stars = get_total_threads(sp->number_stars);
  GPU_PRECISION *host_stars = new GPU_PRECISION[gpu_num_stars * 5];
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


  double *d_qgaus_W = new double[ap->convolve];
  double *d_qgaus_X = new double[ap->convolve];

  d_gauss_legendre(-1.0, 1.0, d_qgaus_X, d_qgaus_W, ap->convolve);
  GPU_PRECISION *host_dx = new GPU_PRECISION[ap->convolve];
  GPU_PRECISION *host_qgaus_W = new GPU_PRECISION[ap->convolve];
  for (int i = 0; i < ap->convolve; i++) {
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
  cutilSafeCall(cudaMalloc((void**) &device_probability_correction, 
			   gpu_num_stars * sizeof(GPU_PRECISION)));
  host_probability = new GPU_PRECISION[gpu_num_stars];

  unsigned int total, free;
  cuMemGetInfo(&free, &total);
  fprintf(stderr, "Used %d/%d memory, %d remaining\n", 
	  (total-free) / 1024, total/1024, free/1024);
}

void gpu__free_constants() {
  for (int i = 0; i < ap->number_integrals; i++) {
    cutilSafeCall(cudaFree(device_V[i]) );
    cutilSafeCall(cudaFree(device_bg_int[i]));
    cutilSafeCall(cudaFree(device_st_int[i]));
#ifdef SINGLE_PRECISION
    cutilSafeCall(cudaFree(device_bg_correction[i]));
    cutilSafeCall(cudaFree(device_st_correction[i]));
#endif
    cutilSafeCall(cudaFree(device_sinb[i]));
    cutilSafeCall(cudaFree(device_sinl[i]));
    cutilSafeCall(cudaFree(device_cosb[i]));
    cutilSafeCall(cudaFree(device_cosl[i]));
    delete [] host_bg_int[i];
    delete [] host_st_int[i];
  }

  delete [] host_bg_int;
  delete [] host_st_int;
  delete [] device_V;
  delete [] device_sinb;
  delete [] device_sinl;
  delete [] device_cosb;
  delete [] device_cosl;
  delete [] device_bg_int;
  delete [] device_st_int;
#ifdef SINGLE_PRECISION
  delete [] device_bg_correction;
  delete [] device_st_correction;
#endif

  cutilSafeCall(cudaFree(device_stars));
  cutilSafeCall(cudaFree(device_probability));
  cutilSafeCall(cudaFree(device_probability_correction));
  delete [] host_probability;
}

template <unsigned int number_streams>
__global__ void 
gpu__zero_integrals(int mu_offset, 
		    int mu_steps, 
		    int nu_steps, 
		    GPU_PRECISION *background_integrals, 
		    GPU_PRECISION *stream_integrals) 
{
  int pos = threadIdx.x + (blockDim.x * blockIdx.x) + mu_offset;
  background_integrals[pos] = 0;
  for (int i = 0; i < number_streams; i++) 
    {
      stream_integrals[pos] = 0;
      pos += (nu_steps * mu_steps);
    }
}

void cpu__sum_integrals(int iteration, double *background_integral, 
			double *stream_integrals) 
{
  int int_size = ap->integral[iteration]->nu_steps * 
    ap->integral[iteration]->mu_steps;
  int gpu_int_size = get_total_threads(int_size);
  cutilSafeCall(cudaMemcpy(host_bg_int[iteration], 
			   device_bg_int[iteration], 
			   gpu_int_size * sizeof(GPU_PRECISION), 
			   cudaMemcpyDeviceToHost) );

  double sum = 0.0;
  for (int i = 0; i < int_size; i++) 
    {
      sum += (double)(host_bg_int[iteration][i]);
      //printf("background_integral[%d/%d]: %.15f\n", i, 
      //int_size, host_bg_int[iteration][i]);
      //printf("(sum(background_integral[%d/%d]: %.15lf\n", i, 
      //integral_size[iteration], sum);
    }
  if (iteration == 0) 
    *background_integral = sum;
  else 
    *background_integral -= sum;

  cutilSafeCall(cudaMemcpy(host_st_int[iteration], 
			   device_st_int[iteration], 
			   ap->number_streams * gpu_int_size *
			   sizeof(GPU_PRECISION), cudaMemcpyDeviceToHost));
  for (int i = 0; i < ap->number_streams; i++) 
    {
      sum = 0.0;
      for (int j = 0; j < int_size; j++) 
	{
	  sum += (double)(host_st_int[iteration]
			  [j + (i * gpu_int_size)]);
	  //printf("[%d] stream_integral: %.15f\n", j, 
	  //host_st_int[iteration][j + (i * int_size)]);
	}
      if (iteration == 0) 
	stream_integrals[i] = sum;
      else 
	stream_integrals[i] -= sum;
    }
}

/********
 *	Likelihood calculation
 ********/

__global__ void 
gpu__zero_likelihood(int offset, 
		     GPU_PRECISION *probability) 
{
  probability[offset + threadIdx.x + (blockDim.x * blockIdx.x)] = 0.0;
}

void cpu__sum_likelihood(double *probability) {
  for (int i = 0; i < sp->number_stars; i++) 
    {
      host_probability[i] = 0;;
    }
  cutilSafeCall(cudaMemcpy(host_probability, device_probability,
			   sp->number_stars * sizeof(GPU_PRECISION),
			   cudaMemcpyDeviceToHost));
  printf("number_stars:%d\n", sp->number_stars);
  for (int i = 0; i < sp->number_stars; i++) 
    {
      *probability += host_probability[i];
    }
}



/********
 *	Run the GPU kernels and get the probability
 ********/

#define stream_parameters(x, y) parameters[(x * 6) + y + 3]
#define stream_weights(x) parameters[(x * 6) + 2]
//#define background_weight parameters[0]
#define background_weight (0.0)
//#define alpha parameters[1]
#define q (parameters[0])
#define r0 (parameters[1])
#define BG_A (parameters[3])
#define BG_B (parameters[4])
#define BG_C (parameters[5])
//#define delta parameters[4]

double gpu__likelihood(double *parameters) {
  double stream_c[3], lbr[3];
  GPU_PRECISION *fstream_a = new GPU_PRECISION[ap->number_streams * 3]; 
  GPU_PRECISION *fstream_c = new GPU_PRECISION[ap->number_streams * 3]; 
  GPU_PRECISION *fstream_sigma_sq2 = new GPU_PRECISION[ap->number_streams]; 

  for (int i = 0; i < ap->number_streams; i++) {
    fstream_sigma_sq2[i] = 1 / (GPU_PRECISION)(2.0 * stream_parameters(i,4) 
					       * stream_parameters(i,4));
    //printf("fstream_sigma_sq2[%d] = %.15f\n", i, fstream_sigma_sq2[i]);
    fstream_a[(i * 3) + 0] = (GPU_PRECISION)(sin(stream_parameters(i,2)) 
					     * cos(stream_parameters(i,3)));
    fstream_a[(i * 3) + 1] = (GPU_PRECISION)(sin(stream_parameters(i,2)) 
					     * sin(stream_parameters(i,3)));
    fstream_a[(i * 3) + 2] = (GPU_PRECISION)(cos(stream_parameters(i,2)));
    
    if (ap->sgr_coordinates == 0) {
      gc_eq_gal(ap->wedge, stream_parameters(i,0) * D_DEG2RAD, 0 * D_DEG2RAD,
		&(lbr[0]), &(lbr[1]));
    } else {
      gc_sgr_gal(ap->wedge, stream_parameters(i,0) * D_DEG2RAD, 0 * D_DEG2RAD,
		 &(lbr[0]), &(lbr[1]));
    }
    lbr[2] = stream_parameters(i,1);
    d_lbr2xyz(lbr, stream_c);

    fstream_c[(i * 3) + 0] = (GPU_PRECISION)stream_c[0]; 
    fstream_c[(i * 3) + 1] = (GPU_PRECISION)stream_c[1];
    fstream_c[(i * 3) + 2] = (GPU_PRECISION)stream_c[2];
    //printf("stream %d fstream_c[0] = %.15f\n", i, fstream_c[i*3]);
    //printf("stream %d fstream_c[1] = %.15f\n", i, fstream_c[i*3 + 1]);
    //printf("stream %d fstream_c[2] = %.15f\n", i, fstream_c[i*3 + 2]);

    //printf("stream %d fstream_a[0] = %.15f\n", i, fstream_a[i*3]);
    //printf("stream %d fstream_a[1] = %.15f\n", i, fstream_a[i*3 + 1]);
    //printf("stream %d fstream_a[2] = %.15f\n", i, fstream_a[i*3 + 2]);
  }
  setup_constant_textures(fstream_a, fstream_c, 
			  fstream_sigma_sq2, ap->number_streams);
  cutilSafeCall(cudaMemcpyToSymbol(constant_inverse_fstream_sigma_sq2,
				   fstream_sigma_sq2, 
				   ap->number_streams * sizeof(GPU_PRECISION),
				   0, cudaMemcpyHostToDevice)); 
  double background_integral = 0.0;
  double *stream_integrals = new double[ap->number_streams];
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
      int total_blocks = total_threads / NUM_THREADS;
      dim3 dimBlock(NUM_THREADS);
      printf("total_blocks:%d\n", total_blocks);
      for(int mu_step = 0; mu_step < total_blocks; mu_step += NUM_BLOCKS) 
	{
	  int offset = mu_step * NUM_THREADS;
	  dim3 dimGrid(min(NUM_BLOCKS, total_blocks - mu_step), R_INCREMENT);
	  //printf("mu_step:%d offset:%d\n", mu_step, offset);
	  EXECUTE_ZERO_INTEGRALS(device_bg_int[i],
				 device_st_int[i]);
#ifdef SINGLE_PRECISION
	  EXECUTE_ZERO_INTEGRALS(device_bg_correction[i],
				 device_st_correction[i]);
#endif
	}
      int aux_bg_profile = ap->aux_bg_profile;
      int shared_mem_size = 0; //shared memory not used in unrolled loops
      for(int j = 0;j<r_steps;j+= R_INCREMENT)
	{
	  boinc_fraction_done(current_work / (double) total_work);
	  for(int mu_step = 0; mu_step < total_blocks; mu_step += NUM_BLOCKS) {
	    dim3 dimGrid(min(NUM_BLOCKS, total_blocks - mu_step), R_INCREMENT);
	    if (aux_bg_profile == 1 && false)
	      {
		EXECUTE_AUX_INTEGRAL_KERNEL;
	      }
	    else
	      {
		EXECUTE_INTEGRAL_KERNEL;
	      }
	    cudaError err = cudaThreadSynchronize();
	    if(err != cudaSuccess)
	      {
		fprintf(stderr, "Error executing gpu__integral_kernel3 error message: %s\n", 
			cudaGetErrorString(err));
		exit(1);
	      }
	  }
	  current_work += R_INCREMENT;
	}
      cpu__sum_integrals(i, &background_integral, stream_integrals);
      printf("background_integral: %.15lf, stream_integral[0]: %.15lf, stream_integral[1]: %.15lf\n", background_integral, stream_integrals[0], stream_integrals[1]);
    }
  cutStopTimer(timer);
  float t = cutGetTimerValue(timer);
  fprintf(stderr, "gpu__integral_kernel3 took %f ms\n", t);
  double *stream_weight = new double[ap->number_streams];
  double exp_weight = exp(background_weight);
  double sum_exp_weights = exp_weight; 
  double bg_weight = exp_weight/background_integral;
  for (int i = 0; i < ap->number_streams; i++) 
    {
      exp_weight = exp(stream_weights(i));
      sum_exp_weights += exp_weight;
      stream_weight[i] = exp_weight/stream_integrals[i];
    }

  GPU_PRECISION f_background_weight[1];
  GPU_PRECISION *f_stream_weight = new GPU_PRECISION[ap->number_streams];
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
  //printf("total_blocks: %d\n", total_blocks);
  for(int mu_step = 0; mu_step < total_blocks; mu_step += NUM_BLOCKS) 
    {
      //printf("mu_step:%d\n", mu_step);
      dim3 dimGrid(min(NUM_BLOCKS, total_blocks - mu_step), R_INCREMENT);
      gpu__zero_likelihood<<<dimGrid, dimBlock>>>
	(mu_step * NUM_THREADS, device_probability);
#ifdef SINGLE_PRECISION
      gpu__zero_likelihood<<<dimGrid, dimBlock>>>
	(mu_step * NUM_THREADS, device_probability_correction);
#endif
      cudaError_t err;
      err = cudaThreadSynchronize();
      if(err != cudaSuccess)
	{
	  fprintf(stderr, "Error executing gpu__zero_likelihood error message: %s\n", 
		  cudaGetErrorString(err));
	  exit(1);
	}
    }
  //printf("got here 3\n");



  //	printf("num streams:%u\n", number_streams);
  for(int mu_step = 0; mu_step < total_blocks; mu_step += NUM_BLOCKS) 
    {
      dim3 dimGrid(min(NUM_BLOCKS, total_blocks - mu_step), R_INCREMENT);
      int convolve = ap->convolve;
      int number_stars = get_total_threads(sp->number_stars);
      int offset = mu_step * NUM_THREADS;
      //printf("mu_step:%d, offset:%d\n", mu_step, offset);
      EXECUTE_LIKELIHOOD_KERNEL;
      cudaError_t err;
      err = cudaThreadSynchronize();
      if(err != cudaSuccess)
	{
	  fprintf(stderr, "Error executing gpu__likelihood_kernel error message: %s\n", 
		  cudaGetErrorString(err));
	  exit(1);
	}
    }
  cutStopTimer(timer);
  t = cutGetTimerValue(timer);
  fprintf(stderr, "gpu__likelihood_kernel took %f ms\n", t);
  cpu__sum_likelihood(&likelihood);
  likelihood /= sp->number_stars;
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
