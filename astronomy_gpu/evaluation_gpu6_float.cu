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

texture<float, 3, cudaReadModeElementType> tex_device_lb;
texture<float, 2, cudaReadModeElementType> tex_fstream_a;
texture<float, 2, cudaReadModeElementType> tex_fstream_c;
texture<float, 2, cudaReadModeElementType> tex_fstream_sigma_sq2;

cudaArray **cu_arrays;

/**
   Similar to setup_texture except it deals with 2d textures
   that were previously in constant memory
 */
void setup_constant_textures(float *fstream_a, float *fstream_c, 
			     float *fstream_sigma_sq2, int number_streams)
{
  // allocate array
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  cudaArray* cu_array_a;
  cutilSafeCall(cudaMallocArray(&cu_array_a, &channelDesc, 3, number_streams)); 
  cutilSafeCall(cudaMemcpyToArray(cu_array_a, 0, 0, fstream_a, 3*number_streams, cudaMemcpyHostToDevice));

  cudaArray* cu_array_c;
  cutilSafeCall(cudaMallocArray(&cu_array_c, &channelDesc, 3, number_streams)); 
  cutilSafeCall(cudaMemcpyToArray(cu_array_c, 0, 0, fstream_c, 3*number_streams, cudaMemcpyHostToDevice));

  cudaArray* cu_array_sq2;
  cutilSafeCall(cudaMallocArray(&cu_array_sq2, &channelDesc, 2, number_streams)); 
  cutilSafeCall(cudaMemcpyToArray(cu_array_sq2, 0, 0, fstream_a, 2*number_streams, cudaMemcpyHostToDevice));
  
  // set texture parameters
  tex_fstream_a.addressMode[0] = cudaAddressModeClamp;
  tex_fstream_a.addressMode[1] = cudaAddressModeClamp;
  tex_fstream_a.filterMode = cudaFilterModePoint;
  tex_fstream_a.normalized = false;

  tex_fstream_c.addressMode[0] = cudaAddressModeClamp;
  tex_fstream_c.addressMode[1] = cudaAddressModeClamp;
  tex_fstream_c.filterMode = cudaFilterModePoint;
  tex_fstream_c.normalized = false;

  tex_fstream_sigma_sq2.addressMode[0] = cudaAddressModeClamp;
  tex_fstream_sigma_sq2.addressMode[1] = cudaAddressModeClamp;
  tex_fstream_sigma_sq2.filterMode = cudaFilterModePoint;
  tex_fstream_sigma_sq2.normalized = false;
  
  // Bind the array to the texture
  cutilSafeCall(cudaBindTextureToArray(tex_fstream_a, cu_array_a, channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_fstream_c, cu_array_c, channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_fstream_sigma_sq2, cu_array_sq2, channelDesc));
}

/**
   Sets up a 3D texture for device__lb, right now it only
   support 1 integral.  In the future it should support
   more by copying from device memory the next integral.
 */
void setup_texture(int mu_steps, int nu_steps, int current_integral, float *host__lb) {
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  //create the array
  cudaArray *cu_array;
  cudaExtent size;
  size.width = 4;
  size.height = mu_steps;
  size.depth = nu_steps;
  cutilSafeCall(cudaMalloc3DArray(&cu_array,&channelDesc,size));
  
  //copy date to the array
  cudaMemcpy3DParms copyParams = {0};
  copyParams.srcPtr   = make_cudaPitchedPtr((void*)host__lb,
					    size.width*sizeof(float), 
					    size.width, 
					    size.height);
  copyParams.dstArray = cu_array;
  copyParams.extent   = size;
  copyParams.kind     = cudaMemcpyHostToDevice;
  cutilSafeCall( cudaMemcpy3D(&copyParams) );
  cu_arrays[current_integral] = cu_array;
}

/**
   Allocates cu arrays for the tex_device_lb texture and sets up
   parts of the texture
*/
void allocate_cu_arrays(int number_integrals) {
  cu_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);
  //set texture parameters
  tex_device_lb.normalized = false;
  tex_device_lb.filterMode = cudaFilterModePoint;
  tex_device_lb.addressMode[0] = cudaAddressModeClamp;
  tex_device_lb.addressMode[1] = cudaAddressModeClamp;
}

void bind_texture(int current_integral) {
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  printf("binding the tex_device_lb texture for integral %u\n",
	 current_integral);
  cutilSafeCall(cudaBindTextureToArray(tex_device_lb, cu_arrays[current_integral], channelDesc));
}


template <unsigned int number_streams, unsigned int convolve> 
__global__ void gpu__integral_kernel3(	int in_step, int in_steps,
		float q, float r0,
		float *device__lb, float *device__V,
		float *background_integrals, float *background_correction, 
		float *stream_integrals, float *stream_correction) {
	int i, j, pos;

	float *shared__r_point = shared_mem;
	float *shared__qw_r3_N = &shared__r_point[convolve];

	//make it so thread 0 populates the shared memory
	//the old way relied on blockDim.x >= convolve
	//which is not always the case
	if (threadIdx.x == 0) {
	  for(unsigned int idx = 0;idx < convolve; ++idx) {
	        pos = (blockIdx.y * convolve * 2) + (idx * 2);

		shared__r_point[idx] = constant__r_constants[pos];
		shared__qw_r3_N[idx] = constant__r_constants[pos + 1];
	  }
	}

	float bg_int, bg_int_correction;
	bg_int = 0.0;
	bg_int_correction = 0.0; 

	float *st_int = &shared__qw_r3_N[convolve];
	float *st_int_correction = &st_int[blockDim.x * number_streams];
	for (i = 0; i < number_streams; i++) {
	  st_int[i * blockDim.x + threadIdx.x] = 0.0f;
	  st_int_correction[i * blockDim.x + threadIdx.x] = 0.0f;
	}

	__syncthreads();

	float corrected_next_term, new_sum;
	float dotted, sxyz0, sxyz1, sxyz2;

//  	float sinb = tex3D(tex_device_lb, 0,
//  			   kernel3__mu_step, kernel3__nu_step);
//  	float sinl = tex3D(tex_device_lb, 1,
//  			   kernel3__mu_step, kernel3__nu_step);
//  	float cosb = tex3D(tex_device_lb, 2,
//  			   kernel3__mu_step, kernel3__nu_step);
//  	float cosl = tex3D(tex_device_lb, 3,
//  			   kernel3__mu_step, kernel3__nu_step);
	
	float zp, rs;
	float xyz0, xyz1, xyz2;
	
	float rg;

	for (i = 0; i < convolve; i++) {
	          xyz2 = shared__r_point[i] * tex3D(tex_device_lb, 0,
						    kernel3__mu_step, kernel3__nu_step);
		  zp = shared__r_point[i] * tex3D(tex_device_lb, 2,
						  kernel3__mu_step, kernel3__nu_step);
		  xyz0 = zp * tex3D(tex_device_lb, 3,
				    kernel3__mu_step, kernel3__nu_step) - f_lbr_r;
		  xyz1 = zp * tex3D(tex_device_lb, 1,
				    kernel3__mu_step, kernel3__nu_step);
		  
		  //__fdividef providers faster fp division, with restrictions on
		  // the fact that (q*q) < 2^126 (appendix b.2.1 in nvidia programming guide)
		  rg = sqrtf(xyz0*xyz0 + xyz1*xyz1 + __fdividef((xyz2*xyz2),(q*q)));
		  rs = rg + r0;
		  

		  
		  corrected_next_term = shared__qw_r3_N[i] / (rg * rs * rs * rs) - bg_int_correction;
		  new_sum = bg_int + corrected_next_term;
		  bg_int_correction = (new_sum - bg_int) - corrected_next_term;
		  bg_int = new_sum;

		  for (j = 0; j < number_streams; j++) {
			pos = (j * 3);
			sxyz0 = xyz0 - constant__fstream_c[pos];
			sxyz1 = xyz1 - constant__fstream_c[pos + 1];
			sxyz2 = xyz2 - constant__fstream_c[pos + 2];

			dotted = constant__fstream_a[pos] * sxyz0 + constant__fstream_a[pos + 1] * sxyz1 + constant__fstream_a[pos + 2] * sxyz2;
			
			sxyz0 -= dotted * constant__fstream_a[pos];
			sxyz1 -= dotted * constant__fstream_a[pos + 1];
			sxyz2 -= dotted * constant__fstream_a[pos + 2];

			pos = j * blockDim.x + threadIdx.x;
			
			corrected_next_term = (shared__qw_r3_N[i] * exp(-((sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2)) / constant__fstream_sigma_sq2[j])) - st_int_correction[pos];
			new_sum = st_int[pos] + corrected_next_term;
			st_int_correction[pos] = (new_sum - st_int[pos]) - corrected_next_term;
			st_int[pos] = new_sum;
		  }
	}
	
	//define V down here so that one to reduce the number of registers, because a register
	//will be reused
	float V = device__V[kernel3__r_step + (kernel3__r_steps * kernel3__nu_step)];
	pos = threadIdx.x + (blockIdx.x * blockDim.x) + (blockIdx.y * gridDim.x * blockDim.x);

	corrected_next_term = (bg_int * V) - background_correction[pos];
	new_sum = background_integrals[pos] + corrected_next_term;	
	background_correction[pos] = (new_sum - background_integrals[pos]) - corrected_next_term;
	background_integrals[pos] = new_sum;
	for (i = 0; i < number_streams; i++) {
		corrected_next_term = (st_int[i * blockDim.x + threadIdx.x] * V) - stream_correction[pos];
		new_sum = stream_integrals[pos] + corrected_next_term;
		stream_correction[pos] = (new_sum - stream_integrals[pos]) - corrected_next_term;
		stream_integrals[pos] = new_sum;

		pos += (blockDim.x * gridDim.x * gridDim.y);
	}
}
