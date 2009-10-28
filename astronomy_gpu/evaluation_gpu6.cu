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
			  //#define MU_STEP_SIZE 1024
#define MU_STEP_SIZE 16000
#define NU_STEP_SIZE 64

#ifdef SINGLE_PRECISION
#define GPU_PRECISION float
#endif
#ifdef DOUBLE_PRECISION
#define GPU_PRECISION double
#endif

int	number_threads = 256;

int	sgr_coordinates;
int	mu_increment = 1;
int	wedge;
int	convolve;
int	number_streams;
int	number_integrals;

int	*r_steps;
int	*mu_steps;
int	*nu_steps;
double     *r_min;
double *r_step_size;
double *mu_min;
double *mu_step_size;
double *nu_min;
double *nu_step_size;

int	*sizeof_V;
GPU_PRECISION	**device__V;		//V				-- float[nu][r]

int	*sizeof_lb;
GPU_PRECISION	**device__lb;		//sinb, sinl, cosb, cosl	-- float[nu][mu][4]
GPU_PRECISION **device__sinb;
GPU_PRECISION **device__sinl;
GPU_PRECISION **device__cosb;
GPU_PRECISION **device__cosl;

int	*integral_size;
GPU_PRECISION	**host__background_integrals;
GPU_PRECISION	**host__stream_integrals;
GPU_PRECISION	**device__background_integrals;
GPU_PRECISION	**device__stream_integrals;

#ifdef SINGLE_PRECISION
GPU_PRECISION **device__background_correction;
GPU_PRECISION **device__stream_correction;
#endif

__device__ __constant__ GPU_PRECISION constant__inverse_fstream_sigma_sq2[8];
__device__ __constant__ GPU_PRECISION constant__fstream_a[12];
__device__ __constant__ GPU_PRECISION constant__fstream_c[12];

__device__ __constant__ GPU_PRECISION constant__dx[MAX_CONVOLVE];
__device__ __constant__ GPU_PRECISION constant__qgaus_W[MAX_CONVOLVE];

__device__ __constant__ GPU_PRECISION constant__background_weight[1];
__device__ __constant__ GPU_PRECISION constant__stream_weight[4];

int	number_stars;
GPU_PRECISION	*device__stars;

int	probability_size;
GPU_PRECISION	*device__probability;
GPU_PRECISION	*device__probability_correction;
GPU_PRECISION	*host__probability;

#define kernel3__mu_step	(mu_offset + blockIdx.x)
#define kernel3__mu_steps	(mu_steps)
//TODO check if blockIdx.y is correct
#define kernel3__r_step		(in_step + blockIdx.y)
#define kernel3__r_steps	in_steps
#define kernel3__nu_step	threadIdx.x
#define kernel3__nu_steps	blockDim.x

//extern "C" void gpu__initialize(ASTRONOMY_PARAMETERS *ap, STAR_POINTS *sp);
//extern "C" double gpu__likelihood(double *parameters);
//extern "C" void gpu__free_constants();

extern __shared__ GPU_PRECISION shared_mem[];

#ifdef SINGLE_PRECISION
#include "evaluation_gpu6_float.cu"
#else
#include "evaluation_gpu6_double.cu"
#endif

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

void gpu__initialize(	int ap_sgr_coordinates, int ap_wedge, int ap_convolve, int ap_number_streams, int ap_number_integrals, 
			int *in__r_steps, double *in__r_min, double *in__r_step_size,
			int *in__mu_steps, double *in__mu_min, double *in__mu_step_size,
			int *in__nu_steps, double *in__nu_min, double *in__nu_step_size,
			int in__number_stars, double **stars) { 
  int i, j;

	sgr_coordinates = ap_sgr_coordinates;
	wedge = ap_wedge;
	convolve = ap_convolve;
	number_streams = ap_number_streams;
	number_integrals = ap_number_integrals;

	sizeof_V = (int*)malloc(number_integrals * sizeof(int));
	sizeof_lb = (int*)malloc(number_integrals * sizeof(int));

	integral_size = (int*)malloc(number_integrals * sizeof(int));

	device__background_integrals = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	device__stream_integrals = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	host__background_integrals = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	host__stream_integrals = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
#ifdef SINGLE_PRECISION
	device__background_correction = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	device__stream_correction = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
#endif

	device__V = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	device__lb = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	device__sinb = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	device__sinl = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	device__cosb = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));
	device__cosl = (GPU_PRECISION**)malloc(number_integrals * sizeof(GPU_PRECISION*));

//	printf("calculating integral constants\n");
	r_steps = (int*)malloc(number_integrals * sizeof(int));
	mu_steps = (int*)malloc(number_integrals * sizeof(int));
	nu_steps = (int*)malloc(number_integrals * sizeof(int));
	nu_min = (double*) malloc(number_integrals * sizeof(double));
	mu_min = (double*) malloc(number_integrals * sizeof(double));
	nu_step_size = (double*) malloc(number_integrals * sizeof(double));
	mu_step_size = (double*) malloc(number_integrals * sizeof(double));
	r_min = (double*) malloc(number_integrals * sizeof(double));
	r_step_size = (double*) malloc(number_integrals * sizeof(double));

	printf("%d integrals %d streams\n", number_integrals, number_streams);

	allocate_cu_arrays(number_integrals);
	for (i = 0; i < number_integrals; i++) {
		r_steps[i] = in__r_steps[i];
		mu_steps[i] = in__mu_steps[i];
		nu_steps[i] = in__nu_steps[i];
		r_min[i] = in__r_min[i];
		r_step_size[i] = in__r_step_size[i];
		nu_min[i] = in__nu_min[i];
		nu_step_size[i] = in__nu_step_size[i];
		mu_min[i] = in__mu_min[i];
		mu_step_size[i] = in__mu_step_size[i];

		sizeof_V[i] = in__nu_steps[i] * in__r_steps[i];
		sizeof_lb[i] = in__mu_steps[i] * in__nu_steps[i] * 4;

		double *cpu__lb;
		double *irv, *reff_xr_rp3, **qw_r3_N, **r_point;
		double *ids, *nus;

		irv		= (double*)malloc(sizeof(double) * r_steps[i]);
		reff_xr_rp3	= (double*)malloc(sizeof(double) * r_steps[i]);
		qw_r3_N		= (double**)malloc(sizeof(double*) * r_steps[i]);
		r_point		= (double**)malloc(sizeof(double*) * r_steps[i]);
		ids		= (double*)malloc(sizeof(double) * nu_steps[i]);
		nus		= (double*)malloc(sizeof(double) * nu_steps[i]);
		printf("populate_cpu__lb\n");
		populate_cpu__lb(sgr_coordinates, wedge, mu_steps[i], mu_min[i],   
				 mu_step_size[i], nu_steps[i], nu_min[i], nu_step_size[i], &cpu__lb);
		printf("cpu__r_constants\n");
		cpu__r_constants(convolve, r_steps[i], r_min[i], r_step_size[i], mu_steps[i], mu_min[i], mu_step_size[i], 
				 nu_steps[i], nu_min[i], nu_step_size[i], irv, r_point, qw_r3_N, reff_xr_rp3, nus, ids);

		GPU_PRECISION *host__V			= (GPU_PRECISION*)malloc(sizeof_V[i] * sizeof(GPU_PRECISION));
		GPU_PRECISION *host__lb			= (GPU_PRECISION*)malloc(sizeof_lb[i] * sizeof(GPU_PRECISION));
		unsigned int sizeof_quarter_lb = sizeof_lb[i] / 4 * sizeof(GPU_PRECISION);
		printf("allocating %d for each lb\n", sizeof_quarter_lb);
		GPU_PRECISION *host__sinb = (GPU_PRECISION*)malloc(sizeof_quarter_lb);
		GPU_PRECISION *host__sinl = (GPU_PRECISION*)malloc(sizeof_quarter_lb);
		GPU_PRECISION *host__cosb = (GPU_PRECISION*)malloc(sizeof_quarter_lb);
		GPU_PRECISION *host__cosl = (GPU_PRECISION*)malloc(sizeof_quarter_lb);
		cutilSafeCall( cudaMalloc((void**) &(device__sinb[i]), sizeof_quarter_lb ));
		cutilSafeCall( cudaMalloc((void**) &(device__sinl[i]), sizeof_quarter_lb ));
		cutilSafeCall( cudaMalloc((void**) &(device__cosb[i]), sizeof_quarter_lb ));
		cutilSafeCall( cudaMalloc((void**) &(device__cosl[i]), sizeof_quarter_lb ));

		long constants_size = 0;
		constants_size += sizeof_V[i] * sizeof(GPU_PRECISION);
		constants_size += sizeof_lb[i] * sizeof(GPU_PRECISION);

//		printf("sizeof_V[%d]: %d\n", i, sizeof_V[i] * sizeof(GPU_PRECISION));
//		printf("sizeof_r_constants[%d]: %d\n", i, sizeof_r_constants[i] * sizeof(GPU_PRECISION));
//		printf("sizeof_lb[%d]: %d\n", i, sizeof_lb[i] * sizeof(GPU_PRECISION));

		printf("Allocating %ld bytes for constants on GPU.\n", constants_size);

		int k;
		for (k = 0; k < nu_steps[i]; k++) {
		  for (j = 0; j < r_steps[i]; j++) {
		    host__V[(j * nu_steps[i]) + k] = reff_xr_rp3[j] * irv[j] * ids[k];
		  }
		}
			       int iteration = 0;
		for (j = 0; j < sizeof_lb[i]; j+=4) {
		  host__sinb[iteration] = cpu__lb[j];
		  host__sinl[iteration] = cpu__lb[j+1];
		  host__cosb[iteration] = cpu__lb[j+2];
		  host__cosl[iteration] = cpu__lb[j+3];
		  ++iteration;
		  host__lb[j] = (GPU_PRECISION)cpu__lb[j];
		}

		//setup_texture(in__mu_steps[i],
		//	      in__nu_steps[i],
		//	      i, host__lb);
		setup_r_point_texture(r_steps[i], convolve,
				      i, r_point);
		setup_qw_r3_N_texture(r_steps[i], convolve,
				      i, qw_r3_N);

//		printf("freeing cpu constants\n");
		free(cpu__lb);
		free(irv);
		free(reff_xr_rp3);
		for(j = 0;j<r_steps[i];++j)
		  {
		    free(qw_r3_N[j]);
		    free(r_point[j]);
		  }
		free(qw_r3_N);
		free(r_point);
		free(ids);
		free(nus);

		printf("device malloc\n");

		cutilSafeCall( cudaMalloc((void**) &(device__V[i]), sizeof_V[i] * sizeof(GPU_PRECISION)) );
		cutilSafeCall( cudaMalloc((void**) &(device__lb[i]), sizeof_lb[i] * sizeof(GPU_PRECISION)) );

//		printf("device memcpy\n");

		cutilSafeCall( cudaMemcpy(device__V[i], host__V, sizeof_V[i] * sizeof(GPU_PRECISION), cudaMemcpyHostToDevice) );
		cutilSafeCall( cudaMemcpy(device__lb[i], host__lb, sizeof_lb[i] * sizeof(GPU_PRECISION), cudaMemcpyHostToDevice) );
		cutilSafeCall( cudaMemcpy(device__sinb[i], host__sinb, sizeof_quarter_lb, cudaMemcpyHostToDevice) );
		cutilSafeCall( cudaMemcpy(device__sinl[i], host__sinl, sizeof_quarter_lb, cudaMemcpyHostToDevice) );
		cutilSafeCall( cudaMemcpy(device__cosb[i], host__cosb, sizeof_quarter_lb, cudaMemcpyHostToDevice) );
		cutilSafeCall( cudaMemcpy(device__cosl[i], host__cosl, sizeof_quarter_lb, cudaMemcpyHostToDevice) );

		free(host__lb);
		free(host__V);
		free(host__sinb);
		free(host__sinl);
		free(host__cosb);
		free(host__cosl);

		integral_size[i] = R_INCREMENT * in__nu_steps[i] * in__mu_steps[i];
//		printf("Allocating %d bytes for integral data on GPU\n", (number_streams + 1) * integral_size[i] * sizeof(GPU_PRECISION));

		cutilSafeCall( cudaMalloc((void**) &device__background_integrals[i], integral_size[i] * sizeof(GPU_PRECISION)) );
		cutilSafeCall( cudaMalloc((void**) &device__stream_integrals[i], number_streams * integral_size[i] * sizeof(GPU_PRECISION)) );
#ifdef SINGLE_PRECISION
		cutilSafeCall( cudaMalloc((void**) &device__background_correction[i], integral_size[i] * sizeof(GPU_PRECISION)) );
		cutilSafeCall( cudaMalloc((void**) &device__stream_correction[i], number_streams * integral_size[i] * sizeof(GPU_PRECISION)) );
#endif
		host__background_integrals[i] = (GPU_PRECISION*)malloc(integral_size[i] * sizeof(GPU_PRECISION));
		host__stream_integrals[i] = (GPU_PRECISION*)malloc(number_streams * integral_size[i] * sizeof(GPU_PRECISION));
	}
	

//	printf("initializing constants for %d stars\n", number_stars);

	number_stars = in__number_stars;
	GPU_PRECISION *host__stars = (GPU_PRECISION*)malloc(number_stars * 5 * sizeof(GPU_PRECISION));
	for (i = 0; i < number_stars; i++) {
		host__stars[i] = (GPU_PRECISION)sin(stars[i][1] * D_DEG2RAD);
		host__stars[i + number_stars*1] = (GPU_PRECISION)sin(stars[i][0] * D_DEG2RAD);
		host__stars[i + number_stars*2] = (GPU_PRECISION)cos(stars[i][1] * D_DEG2RAD);
		host__stars[i + number_stars*3] = (GPU_PRECISION)cos(stars[i][0] * D_DEG2RAD);
		host__stars[i + number_stars*4] = (GPU_PRECISION)stars[i][2];
	}
//	printf("allocating %d bytes for device__stars\n", number_stars * 5 * sizeof(GPU_PRECISION));
	cutilSafeCall( cudaMalloc((void**) &device__stars, number_stars * 5 * sizeof(GPU_PRECISION)) );
	cutilSafeCall( cudaMemcpy(device__stars, host__stars, number_stars * 5 * sizeof(GPU_PRECISION), cudaMemcpyHostToDevice) );

	free(host__stars);

//	printf("initializing qgaus constants\n");

	double *d_qgaus_W = (double*)malloc(convolve * sizeof(double));
	double *d_qgaus_X = (double*)malloc(convolve * sizeof(double));

	d_gauss_legendre(-1.0, 1.0, d_qgaus_X, d_qgaus_W, convolve);
	GPU_PRECISION *host__dx = (GPU_PRECISION*)malloc(convolve * sizeof(GPU_PRECISION));
	GPU_PRECISION *host__qgaus_W = (GPU_PRECISION*)malloc(convolve * sizeof(GPU_PRECISION));
	for (i = 0; i < convolve; i++) {
		host__dx[i] = (GPU_PRECISION)(3.0 * d_stdev * d_qgaus_X[i]);
		host__qgaus_W[i] = (GPU_PRECISION)d_qgaus_W[i];
	}
	free(d_qgaus_W);
	free(d_qgaus_X);

	cutilSafeCall( cudaMemcpyToSymbol(constant__dx, host__dx, convolve * sizeof(GPU_PRECISION), 0, cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpyToSymbol(constant__qgaus_W, host__qgaus_W, convolve * sizeof(GPU_PRECISION), 0, cudaMemcpyHostToDevice) );

	//	free(host__dx);
	//	free(host__qgaus_W);

//	printf("mallocing GPU bg and stream probability: %d bytes\n", number_threads * sizeof(GPU_PRECISION));

	probability_size = number_threads;
	cutilSafeCall( cudaMalloc((void**) &device__probability, probability_size * sizeof(GPU_PRECISION)) );
	cutilSafeCall( cudaMalloc((void**) &device__probability_correction, probability_size * sizeof(GPU_PRECISION)) );

//	printf("mallocing host bg and stream probability\n");

	host__probability = (GPU_PRECISION*)malloc(probability_size * sizeof(GPU_PRECISION));
	free(host__dx);
	free(host__qgaus_W);

	unsigned int total, free;
	cuMemGetInfo(&free, &total);
	printf("Used %d/%d memory, %d remaining\n", (total-free) / 1024, total/1024, free/1024);
}

void gpu__free_constants() {
	int i;
	for (i = 0; i < number_integrals; i++) {
		cutilSafeCall( cudaFree(device__V[i]) );
		cutilSafeCall( cudaFree(device__lb[i]) );
		cutilSafeCall( cudaFree(device__background_integrals[i]) );
		cutilSafeCall( cudaFree(device__stream_integrals[i]) );
#ifdef SINGLE_PRECISION
		cutilSafeCall( cudaFree(device__background_correction[i]) );
		cutilSafeCall( cudaFree(device__stream_correction[i]) );
#endif
		free(host__background_integrals[i]);
	}

	cutilSafeCall( cudaFree(device__stars) );

	free(host__background_integrals);
	free(device__V);
	free(device__lb);
	free(device__background_integrals);
	free(device__stream_integrals);
#ifdef SINGLE_PRECISION
	free(device__background_correction);
	free(device__stream_correction);
#endif

	free(r_steps);
	free(mu_steps);
	free(nu_steps);

	free(sizeof_V);
	free(sizeof_lb);
	free(integral_size);

	cutilSafeCall( cudaFree(device__probability) );
	free(host__probability);

}

template <unsigned int number_streams>
__global__ void gpu__zero_integrals(int offset, int mu_steps, GPU_PRECISION *background_integrals, GPU_PRECISION *stream_integrals) {
  int pos = threadIdx.x + ((offset + blockIdx.x) * blockDim.x) + (blockIdx.y * blockDim.x * mu_steps);

	background_integrals[pos] = 0;
	for (int i = 0; i < number_streams; i++) stream_integrals[(i * gridDim.y * mu_steps * blockDim.x) + pos] = 0;
}

void cpu__sum_integrals(int iteration, double *background_integral, double *stream_integrals) {
  int i, j;

	cutilSafeCall( cudaMemcpy(host__background_integrals[iteration], device__background_integrals[iteration], integral_size[iteration] * sizeof(GPU_PRECISION), cudaMemcpyDeviceToHost) );

	double sum = 0.0;
	for (i = 0; i < integral_size[iteration]; i++) {
		sum += (double)(host__background_integrals[iteration][i]);
		//printf("background_integral[%d/%d]: %.15f\n", i, integral_size[iteration], host__background_integrals[iteration][i]);	  
		//printf("(sum(background_integral[%d/%d]: %.15lf\n", i, integral_size[iteration], sum);
	}
	if (iteration == 0) *background_integral = sum;
	else *background_integral -= sum;

 	cutilSafeCall( cudaMemcpy(host__stream_integrals[iteration], device__stream_integrals[iteration], number_streams * integral_size[iteration] * sizeof(GPU_PRECISION), cudaMemcpyDeviceToHost) );
	for (i = 0; i < number_streams; i++) {
		sum = 0.0;
		for (j = 0; j < integral_size[iteration]; j++) {

		  sum += (double)(host__stream_integrals[iteration]
				  [j + (i * integral_size[iteration])]);
		  //printf("stream_integral: %.15f\n", host__stream_integrals[iteration][j + (i * integral_size[iteration])]);
		}
		if (iteration == 0) stream_integrals[i] = sum;
		else stream_integrals[i] -= sum;
	}
}

/********
 *	Likelihood calculation
 ********/

template <unsigned int number_streams>
__global__ void gpu__zero_likelihood(int block_size, GPU_PRECISION *device__probability) {
	device__probability[threadIdx.x] = 0;
}

void cpu__sum_likelihood(int block_size, double *probability) {
	int i;

	cutilSafeCall( cudaMemcpy(host__probability, device__probability, probability_size * sizeof(GPU_PRECISION), cudaMemcpyDeviceToHost) );

//	*probability = 0.0;

	for (i = 0; i < block_size; i++) {
		*probability += host__probability[i];
//		printf("bg_prob: %.15f\n", host__probability[i]);
	}
}

#ifndef SINGLE_PRECISION
template <unsigned int number_streams>
__global__ void gpu__likelihood_kernel(	int offset, int convolve,
					GPU_PRECISION q_squared_inverse, GPU_PRECISION r0,
					GPU_PRECISION coeff, 
					GPU_PRECISION *device__stars,
					int number_stars,
					GPU_PRECISION *probability) {
#else
  template <unsigned int number_streams>
    __global__ void gpu__likelihood_kernel(int offset, int convolve,
					   GPU_PRECISION q_squared_inverse, GPU_PRECISION r0,
					   GPU_PRECISION coeff, 
					   GPU_PRECISION *device__stars,
					   int number_stars,
					   GPU_PRECISION *probability, GPU_PRECISION *probability_correction) {
#endif
	int i;
	int pos = (offset + threadIdx.x);
	GPU_PRECISION sinb = device__stars[pos];
	GPU_PRECISION sinl = device__stars[pos + number_stars];
	GPU_PRECISION cosb = device__stars[pos + number_stars*2];
	GPU_PRECISION cosl = device__stars[pos + number_stars*3];
	GPU_PRECISION coords = device__stars[pos + number_stars*4];

	GPU_PRECISION rg, xyz0, xyz1, xyz2;
	GPU_PRECISION dotted, sxyz0, sxyz1, sxyz2;

	GPU_PRECISION gPrime = 5.0 * (log10(coords * 1000.0) - 1.0) + d_absm;
	GPU_PRECISION exponent = exp(sigmoid_curve_1 * (gPrime - sigmoid_curve_2));
	GPU_PRECISION reff_value = sigmoid_curve_0 / (exponent + 1);
	GPU_PRECISION rPrime3 = coords * coords * coords;

#ifdef __DEVICE_EMULATION__
	printf("%d gPrime:%.30lf exponent:%.30lf reff_value:%.30lf rPrime3:%.30lf\n",
	       pos, gPrime, exponent, reff_value, rPrime3);
#endif

	GPU_PRECISION reff_xr_rp3 = reff_value * d_xr / rPrime3;

	GPU_PRECISION r_point, qw_r3_N;
	GPU_PRECISION zp, rs, g;

#ifndef SINGLE_PRECISION
	GPU_PRECISION bg_int = 0.0;
	GPU_PRECISION st_int[number_streams];
	for (i = 0; i < number_streams; i++) st_int[i] = 0.0;
#else
	GPU_PRECISION bg_int, bg_int_correction;
	bg_int = 0.0;
	bg_int_correction = 0.0; 

	GPU_PRECISION st_int[number_streams];
	GPU_PRECISION st_int_correction[number_streams];
	for (i = 0; i < number_streams; i++) {
		st_int[i] = 0;
		st_int_correction[i] = 0.0;
	}

	GPU_PRECISION corrected_next_term, new_sum;
#endif

	for (i = 0; i < convolve; i++) {
		g = gPrime + constant__dx[i];
#ifdef SINGLE_PRECISION
		r_point = pow(10.0f, (g - f_absm)/5.0f + 1.0f) / 1000.0f;
#else
		r_point = pow(10.0, (g - d_absm) * .2 + 1.0) * .001;
#endif
		rPrime3 = r_point * r_point * r_point;

		qw_r3_N = constant__qgaus_W[i] * rPrime3 * coeff * exp( -((g - gPrime) * (g - gPrime) / (2 * d_stdev * d_stdev)) );

		xyz2 = r_point * sinb;
		zp = r_point * cosb;
		xyz0 = zp * cosl - d_lbr_r;
		xyz1 = zp * sinl;

		rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_squared_inverse);
		rs = rg + r0;


#ifndef SINGLE_PRECISION
		bg_int += (qw_r3_N / (rg * rs * rs * rs));
#else
		corrected_next_term = (qw_r3_N / (rg * rs * rs * rs)) - bg_int_correction;
		new_sum = bg_int + corrected_next_term;
		bg_int_correction = (new_sum - bg_int) - corrected_next_term;
		bg_int = new_sum;
#endif

		for (int j = 0; j < number_streams; j++) {
			pos = (j * 3);
 			sxyz0 = xyz0 - constant__fstream_c[pos + 0];
 			sxyz1 = xyz1 - constant__fstream_c[pos + 1];
 			sxyz2 = xyz2 - constant__fstream_c[pos + 2];

 			dotted = constant__fstream_a[pos + 0] * sxyz0 
			  + constant__fstream_a[pos + 1] * sxyz1
			  + constant__fstream_a[pos + 2] * sxyz2;
			
 			sxyz0 -= dotted * constant__fstream_a[pos + 0];
 			sxyz1 -= dotted * constant__fstream_a[pos + 1];
 			sxyz2 -= dotted * constant__fstream_a[pos + 2];
			
#ifndef SINGLE_PRECISION
			st_int[j] += qw_r3_N * exp(-((sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2)) * constant__inverse_fstream_sigma_sq2[j]);
#else
			corrected_next_term = (qw_r3_N * exp(-((sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2)) * constant__inverse_fstream_sigma_sq2[j])) - st_int_correction[j];
			new_sum = st_int[j] + corrected_next_term;
			st_int_correction[j] = (new_sum - st_int[j]) - corrected_next_term;
			st_int[j] = new_sum;
#endif
		}
	}
	GPU_PRECISION probability_sum = 0.0;
	probability_sum += bg_int * constant__background_weight[0];
	//pragma unroll 1 makes the loop not unroll,
	//when it unrolls it causes a launch failure when trying
	//to access constant__stream_weight[i], when i is 1
	#pragma unroll 1
	for (i = 0; i < number_streams; i++) {
		probability_sum += st_int[i] * constant__stream_weight[i];
	}
	probability_sum *= reff_xr_rp3;
//	printf("bg_prob %.15f st_prob[0]: %.15f st_prob[1]: %.15f, prob_sum: %.15f\n", (bg_int * reff_xr_rp3), (st_int[0] * reff_xr_rp3), (st_int[1] * reff_xr_rp3), probability_sum);

	if (probability_sum == 0.0) probability_sum = -238.0;
	else probability_sum = log10(probability_sum);

#ifndef SINGLE_PRECISION
	probability[threadIdx.x] += probability_sum;
#else
	pos = threadIdx.x;
	corrected_next_term = probability_sum - probability_correction[pos];
	new_sum = probability[pos] + corrected_next_term;
	probability_correction[pos] = (new_sum - probability[pos]) - corrected_next_term;
	probability[pos] = new_sum;
#endif
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
//#define delta parameters[4]

double gpu__likelihood(double *parameters) {
	int i, j;

	double stream_c[3], lbr[3];
	GPU_PRECISION *fstream_a = (GPU_PRECISION*) malloc(sizeof(GPU_PRECISION) * number_streams * 3);
	GPU_PRECISION *fstream_c = (GPU_PRECISION*) malloc(sizeof(GPU_PRECISION) * number_streams * 3);
	GPU_PRECISION *fstream_sigma_sq2 = (GPU_PRECISION*) malloc(sizeof(GPU_PRECISION) * number_streams);

	for (i = 0; i < number_streams; i++) {
		fstream_sigma_sq2[i] = 1 / (GPU_PRECISION)(2.0 * stream_parameters(i,4) * stream_parameters(i,4));

		fstream_a[(i * 3)] = (GPU_PRECISION)( sin(stream_parameters(i,2)) * cos(stream_parameters(i,3)) );
		fstream_a[(i * 3) + 1] = (GPU_PRECISION)( sin(stream_parameters(i,2)) * sin(stream_parameters(i,3)) );
		fstream_a[(i * 3) + 2] = (GPU_PRECISION)( cos(stream_parameters(i,2)) );

		if (sgr_coordinates == 0) {
			gc_eq_gal(wedge, stream_parameters(i,0) * D_DEG2RAD, 0 * D_DEG2RAD, &(lbr[0]), &(lbr[1]));
		} else {
			gc_sgr_gal(wedge, stream_parameters(i,0) * D_DEG2RAD, 0 * D_DEG2RAD, &(lbr[0]), &(lbr[1]));
		}
		lbr[2] = stream_parameters(i,1);
		d_lbr2xyz(lbr, stream_c);

		fstream_c[(i * 3)] = (GPU_PRECISION)stream_c[0]; 
		fstream_c[(i * 3) + 1] = (GPU_PRECISION)stream_c[1];
		fstream_c[(i * 3) + 2] = (GPU_PRECISION)stream_c[2];
	}

	setup_constant_textures(fstream_a, fstream_c, fstream_sigma_sq2,  number_streams);
	cutilSafeCall( cudaMemcpyToSymbol(constant__inverse_fstream_sigma_sq2, fstream_sigma_sq2, number_streams * sizeof(GPU_PRECISION), 0, cudaMemcpyHostToDevice) ); 
	cutilSafeCall( cudaMemcpyToSymbol(constant__fstream_a, fstream_a, number_streams * 3 * sizeof(GPU_PRECISION), 0, cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpyToSymbol(constant__fstream_c, fstream_c, number_streams * 3 * sizeof(GPU_PRECISION), 0, cudaMemcpyHostToDevice) );

	double background_integral = 0.0;
	double *stream_integrals = (double*)malloc(number_streams * sizeof(double));
	for (i = 0; i < number_streams; i++) stream_integrals[i] = 0.0;

	double coeff = 1.0 / (d_stdev * sqrt(2.0 * D_PI));

	unsigned int work_amount = 0;
	unsigned int timer = 0;
	cutCreateTimer(&timer);
	cutResetTimer(timer);
	cutStartTimer(timer);
	unsigned int total_work = 0;
	unsigned int current_work = 0;
	for (i = 0; i < number_integrals; i++) {
		total_work += r_steps[i];
	}

	double q_squared_inverse = 1 / (q * q);
	for (i = 0; i < number_integrals; i++) {
		bind_texture(i);
//		printf("mu_steps[i] is %u nu_steps[i] is %u rsteps[i] is %u convolve is %u number_streams is %u\n", 
//		       mu_steps[i], nu_steps[i], r_steps[i], convolve, number_streams);
		work_amount += mu_steps[i] * number_streams;

		
		//split mu_steps in order to get a smaller grid size to 
		//reduce the UI lag that occurs
		for(int mu_step = 0; mu_step < mu_steps[i]; mu_step += MU_STEP_SIZE) {
		  dim3 dimGrid(min(MU_STEP_SIZE, mu_steps[i] - mu_step), R_INCREMENT);
		  switch(number_streams) {
		  case 1:gpu__zero_integrals<1><<<dimGrid, nu_steps[i]>>>(mu_step, mu_steps[i], device__background_integrals[i], device__stream_integrals[i]);break;
		  case 2:gpu__zero_integrals<2><<<dimGrid, nu_steps[i]>>>(mu_step, mu_steps[i], device__background_integrals[i], device__stream_integrals[i]);break;
		  case 3:gpu__zero_integrals<3><<<dimGrid, nu_steps[i]>>>(mu_step, mu_steps[i], device__background_integrals[i], device__stream_integrals[i]);break;
		  case 4:gpu__zero_integrals<4><<<dimGrid, nu_steps[i]>>>(mu_step, mu_steps[i], device__background_integrals[i], device__stream_integrals[i]);break;
		  };
#ifdef SINGLE_PRECISION
		  switch(number_streams) {
		  case 1:gpu__zero_integrals<1><<<dimGrid, nu_steps[i]>>>(mu_step, mu_steps[i], device__background_correction[i], device__stream_correction[i]);break;
		  case 2:gpu__zero_integrals<2><<<dimGrid, nu_steps[i]>>>(mu_step, mu_steps[i], device__background_correction[i], device__stream_correction[i]);break;
		  case 3:gpu__zero_integrals<3><<<dimGrid, nu_steps[i]>>>(mu_step, mu_steps[i], device__background_correction[i], device__stream_correction[i]);break;
		  case 4:gpu__zero_integrals<4><<<dimGrid, nu_steps[i]>>>(mu_step, mu_steps[i], device__background_correction[i], device__stream_correction[i]);break;
		  };
#endif
		}
		cudaError_t err;
		err = cudaThreadSynchronize();
		if(err != cudaSuccess)
		  {
		    fprintf(stderr, "Error executing gpu__zero_integrals error message: %s\n", 
			    cudaGetErrorString(err));
		  }

		unsigned int sh_mem_multiple;
#ifdef SINGLE_PRECISION
		sh_mem_multiple = 2;
#else
		sh_mem_multiple = 1;
#endif
		int num_blocks = (mu_steps[i] * nu_steps[i]) / NU_STEP_SIZE;
		unsigned int shared_mem_size;
		int nu_step_size = NU_STEP_SIZE;
		if (num_blocks * NU_STEP_SIZE == mu_steps[i] * nu_steps[i])
		  {
		    //NU_STEP_SIZE divides the number of steps
		  }
		else
		  {
		    //NU_STEP_SIZE does not divide the number of steps
		    //so default to the old behaviur
		    num_blocks = mu_steps[i];
		    nu_step_size = nu_steps[i];
		  }
		shared_mem_size = sh_mem_multiple * (nu_step_size * sizeof(GPU_PRECISION) * number_streams);
		shared_mem_size += 192;
		printf("%d blocks %d threads\n", mu_steps[i], nu_steps[i]);
		printf("Allocating %d bytes for shared memory\n", shared_mem_size);
		
		dim3 dimBlock(nu_step_size);
		printf("Using %d blocks and %d threads\n", num_blocks, nu_step_size);
		//r_steps[i] = 1;
		for(j = 0;j<r_steps[i];j+= R_INCREMENT)
		  {
		    boinc_fraction_done(current_work / (double) total_work);
		    for(int mu_step = 0; mu_step < num_blocks; mu_step += MU_STEP_SIZE) {
		      dim3 dimGrid(min(MU_STEP_SIZE, num_blocks - mu_step), R_INCREMENT);
#ifndef SINGLE_PRECISION
			switch(number_streams) {
			case 1:	gpu__integral_kernel3<1, MAX_CONVOLVE><<<dimGrid, dimBlock, shared_mem_size>>>(mu_step, mu_steps[i],
													       j, r_steps[i],
													       nu_steps[i],
													       q_squared_inverse, r0,
													       device__sinb[i],
													       device__sinl[i],
													       device__cosb[i],
													       device__cosl[i],
													       device__V[i],
													       device__background_integrals[i], device__stream_integrals[i]);
			  break;
			case 2:	gpu__integral_kernel3<2, MAX_CONVOLVE><<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],
													       j, r_steps[i],
													       nu_steps[i],
													       q_squared_inverse, r0,
													       device__sinb[i],
													       device__sinl[i],
													       device__cosb[i],
													       device__cosl[i],
													       device__V[i],
													       device__background_integrals[i], device__stream_integrals[i]);
			  break;
			case 3:	gpu__integral_kernel3<3, MAX_CONVOLVE><<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],
													       j, r_steps[i],
													       nu_steps[i],
													       q_squared_inverse, r0,
													       device__sinb[i],
													       device__sinl[i],
													       device__cosb[i],
													       device__cosl[i],
													       device__V[i],
													       device__background_integrals[i], device__stream_integrals[i]);
			  break;
			case 4:	gpu__integral_kernel3<4, MAX_CONVOLVE><<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],
													       j, r_steps[i],
													       nu_steps[i],
													       q_squared_inverse, r0,
													       device__sinb[i],
													       device__sinl[i],
													       device__cosb[i],
													       device__cosl[i],
													       device__V[i],
													       device__background_integrals[i], device__stream_integrals[i]);
			  break;
			}
			
#else
			switch(number_streams) {
			case 1:	
			  gpu__integral_kernel3<1, MAX_CONVOLVE>
			    <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],
								     j, r_steps[i], 
								     nu_steps[i],
								     q_squared_inverse, r0,
								     device__sinb[i],
								     device__sinl[i],
								     device__cosb[i],
								     device__cosl[i],
								     device__V[i],
								     device__background_integrals[i], 
								     device__background_correction[i], 
								     device__stream_integrals[i], 
								     device__stream_correction[i]);
			  break;
			case 2:	
			  gpu__integral_kernel3<2, MAX_CONVOLVE>
			    <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],
								     j, r_steps[i], 
								     nu_steps[i],
								     q_squared_inverse, r0,
								     device__sinb[i],
								     device__sinl[i],
								     device__cosb[i],
								     device__cosl[i],
								     device__V[i],
								     device__background_integrals[i], 
								     device__background_correction[i], 
								     device__stream_integrals[i], 
								     device__stream_correction[i]);
			  break;
			case 3:	
			  gpu__integral_kernel3<3, MAX_CONVOLVE>
			    <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],
								     j, r_steps[i], 
								     nu_steps[i],
								     q_squared_inverse, r0,
								     device__sinb[i],
								     device__sinl[i],
								     device__cosb[i],
								     device__cosl[i],
								     device__V[i],
								     device__background_integrals[i], 
								     device__background_correction[i], 
								     device__stream_integrals[i], 
								     device__stream_correction[i]);
			  break;
			case 4:	
			  gpu__integral_kernel3<4, MAX_CONVOLVE>
			    <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],
								     j, r_steps[i], 
								     nu_steps[i],
								     q_squared_inverse, r0,
								     device__sinb[i],
								     device__sinl[i],
								     device__cosb[i],
								     device__cosl[i],
								     device__V[i],
								     device__background_integrals[i], 
								     device__background_correction[i], 
								     device__stream_integrals[i], 
								     device__stream_correction[i]);			  
			  break;
			}
#endif
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
		cudaThreadSynchronize();
		printf("background_integral: %.30lf, stream_integral[0]: %.30lf, stream_integral[1]: %.30lf\n", background_integral, stream_integrals[0], stream_integrals[1]);
	}
	//background_integral = 0.00047190245735466022957219078826085478794993832707;
	//stream_integrals[0] = 10.71311528236298293847994500538334250450134277343750;
	//stream_integrals[1] = 513.27062657288058744597947224974632263183593750000000;
	cutStopTimer(timer);
	float t = cutGetTimerValue(timer);
	printf("gpu__integral_kernel3 took %f ms\n", t);
	int block_size;

	double *stream_weight = (double*)malloc(number_streams * sizeof(double));
	double exp_weight = exp(background_weight);
	double sum_exp_weights = exp_weight; 
	double bg_weight = exp_weight/background_integral;
	for (i = 0; i < number_streams; i++) {
		exp_weight = exp(stream_weights(i));
		sum_exp_weights += exp_weight;
		stream_weight[i] = exp_weight/stream_integrals[i];
	}

	GPU_PRECISION f_background_weight[1];
	GPU_PRECISION *f_stream_weight = (GPU_PRECISION*)malloc(number_streams * sizeof(GPU_PRECISION));
	f_background_weight[0] = (GPU_PRECISION)( bg_weight / sum_exp_weights );
	for (i = 0; i < number_streams; i++) {
		f_stream_weight[i] = (GPU_PRECISION)( stream_weight[i] / sum_exp_weights );
	}

	cutilSafeCall( cudaMemcpyToSymbol(constant__background_weight, f_background_weight, 1 * sizeof(GPU_PRECISION), 0, cudaMemcpyHostToDevice) );
	cutilSafeCall( cudaMemcpyToSymbol(constant__stream_weight, f_stream_weight, number_streams * sizeof(GPU_PRECISION), 0, cudaMemcpyHostToDevice) );

	double likelihood = 0.0;
	cutResetTimer(timer);
	cutStartTimer(timer);
	switch(number_streams) {
	case 1: gpu__zero_likelihood<1><<<1, number_threads>>>(number_threads, device__probability);break;
	case 2:	gpu__zero_likelihood<2><<<1, number_threads>>>(number_threads, device__probability);break;
	case 3:	gpu__zero_likelihood<3><<<1, number_threads>>>(number_threads, device__probability);break;
	case 4:	gpu__zero_likelihood<4><<<1, number_threads>>>(number_threads, device__probability);break;
	};
#ifdef SINGLE_PRECISION
	switch(number_streams) {
	case 1: gpu__zero_likelihood<1><<<1, number_threads>>>(number_threads, device__probability_correction);break;
	case 2: gpu__zero_likelihood<2><<<1, number_threads>>>(number_threads, device__probability_correction);break;
	case 3: gpu__zero_likelihood<3><<<1, number_threads>>>(number_threads, device__probability_correction);break;
	case 4: gpu__zero_likelihood<4><<<1, number_threads>>>(number_threads, device__probability_correction);break;
	};
#endif

	cudaError_t err;
	err = cudaThreadSynchronize();
	if(err != cudaSuccess)
	  {
	    fprintf(stderr, "Error executing gpu__zero_likelihood error message: %s\n", 
		    cudaGetErrorString(err));
	    exit(1);
	  }
//	printf("num streams:%u\n", number_streams);
	for (i = 0; i < number_stars; i += number_threads) {
		block_size = min(number_threads, number_stars - i);
#ifndef SINGLE_PRECISION
		switch (number_streams) {
			case 1:	gpu__likelihood_kernel<1><<<1, block_size>>>(	i, convolve,
										q_squared_inverse, r0,
										(GPU_PRECISION)coeff,
										device__stars,
										number_stars,
										device__probability);
			break;
			case 2:	gpu__likelihood_kernel<2><<<1, block_size>>>(	i, convolve,
										q_squared_inverse, r0,
										(GPU_PRECISION)coeff,
										device__stars,
										number_stars,
										device__probability);
			break;
			case 3:	gpu__likelihood_kernel<3><<<1, block_size>>>(	i, convolve,
										q_squared_inverse, r0,
										(GPU_PRECISION)coeff,
										device__stars,
										number_stars,
										device__probability);
			break;
			case 4:	gpu__likelihood_kernel<4><<<1, block_size>>>(	i, convolve,
										q_squared_inverse, r0,
										(GPU_PRECISION)coeff,
										device__stars,
										number_stars,
										device__probability);
			break;
		}

#else
		switch (number_streams) {
		case 1:	gpu__likelihood_kernel<1><<<1, block_size>>>(	i, convolve,
									q_squared_inverse, r0,
									(GPU_PRECISION)coeff,
									device__stars,
									number_stars,
									device__probability, device__probability_correction);
		  break;
		case 2:	gpu__likelihood_kernel<2><<<1, block_size>>>(	i, convolve,
									q_squared_inverse, r0,
									(GPU_PRECISION)coeff,
									device__stars,
									number_stars,
									device__probability, device__probability_correction);
		  break;
		case 3:	gpu__likelihood_kernel<3><<<1, block_size>>>(	i, convolve,
									q_squared_inverse, r0,
									(GPU_PRECISION)coeff,
									device__stars,
									number_stars,
									device__probability, device__probability_correction);
		  break;
		case 4:	gpu__likelihood_kernel<4><<<1, block_size>>>(	i, convolve,
									q_squared_inverse, r0,
									(GPU_PRECISION)coeff,
									device__stars,
									number_stars,
									device__probability, device__probability_correction);
		  break;
		}
#endif
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
	printf("gpu__likelihood_kernel took %f ms\n", t);
	cpu__sum_likelihood(number_threads, &likelihood);
	likelihood /= number_stars;
	printf("likelihood: %.30lf\n", likelihood);

	free(fstream_a);
	free(fstream_c);
	free(fstream_sigma_sq2);
	free(f_stream_weight);
	free(stream_integrals);
	free(stream_weight);
	cutDeleteTimer(timer);

	return likelihood;

}
