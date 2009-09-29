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
texture<float, 2, cudaReadModeElementType> tex_r_point;
texture<float, 2, cudaReadModeElementType> tex_qw_r3_N;
texture<float, 2, cudaReadModeElementType> tex_fstream_a;
texture<float, 2, cudaReadModeElementType> tex_fstream_c;
texture<float, 2, cudaReadModeElementType> tex_fstream_sigma_sq2;

cudaArray **cu_arrays;
cudaArray **cu_r_point_arrays;
cudaArray **cu_qw_r3_N_arrays;

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
  cutilSafeCall(cudaMemcpyToArray(cu_array_sq2, 0, 0, fstream_sigma_sq2, 2*number_streams, cudaMemcpyHostToDevice));
  
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

void setup_r_point_texture(int r_steps, int convolve, int current_integral, double **r_point)
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  cudaArray* cu_array;
  float *r_point_flat = (float*) malloc(sizeof(float) * r_steps * convolve);
  int i,j;
  for(i = 0;i<r_steps;++i)
    {
      for(j = 0;j<convolve;++j)
	{
	  r_point_flat[i * convolve + j] = 
	    r_point[i][j];
	}
    }
  cutilSafeCall(cudaMallocArray(&cu_array, &channelDesc, convolve, r_steps)); 
  cutilSafeCall(cudaMemcpyToArray(cu_array, 0, 0, r_point_flat, r_steps * convolve * sizeof(float), cudaMemcpyHostToDevice));
  cu_r_point_arrays[current_integral] = cu_array;
  free(r_point_flat);
}

void setup_qw_r3_N_texture(int r_steps, int convolve, int current_integral, double **qw_r3_N)
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  cudaArray* cu_array;
  float *qw_r3_N_flat = (float*) malloc(sizeof(float) * r_steps * convolve);
  int i,j;
  for(i = 0;i<r_steps;++i)
    {
      for(j = 0;j<convolve;++j)
	{
	  qw_r3_N_flat[i * convolve + j] = 
	    qw_r3_N[i][j];
	}
    }
  cutilSafeCall(cudaMallocArray(&cu_array, &channelDesc, convolve, r_steps)); 
  cutilSafeCall(cudaMemcpyToArray(cu_array, 0, 0, qw_r3_N_flat, r_steps * convolve * sizeof(float), cudaMemcpyHostToDevice));
  cu_qw_r3_N_arrays[current_integral] = cu_array;
  free(qw_r3_N_flat);
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
  size.height = nu_steps;
  size.depth = mu_steps;
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
  cu_r_point_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);
  cu_qw_r3_N_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);
  //set texture parameters
  tex_device_lb.normalized = false;
  tex_device_lb.filterMode = cudaFilterModePoint;
  tex_device_lb.addressMode[0] = cudaAddressModeClamp;
  tex_device_lb.addressMode[1] = cudaAddressModeClamp;

  tex_r_point.addressMode[0] = cudaAddressModeClamp;
  tex_r_point.addressMode[1] = cudaAddressModeClamp;
  tex_r_point.filterMode = cudaFilterModePoint;
  tex_r_point.normalized = false;

  tex_qw_r3_N.addressMode[0] = cudaAddressModeClamp;
  tex_qw_r3_N.addressMode[1] = cudaAddressModeClamp;
  tex_qw_r3_N.filterMode = cudaFilterModePoint;
  tex_qw_r3_N.normalized = false;
}

void bind_texture(int current_integral) {
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  printf("binding the tex_device_lb texture for integral %u\n",
	 current_integral);
  cutilSafeCall(cudaBindTextureToArray(tex_device_lb, cu_arrays[current_integral], channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_r_point, cu_r_point_arrays[current_integral], channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_qw_r3_N, cu_qw_r3_N_arrays[current_integral], channelDesc));
}


template <unsigned int number_streams, unsigned int convolve> 
__global__ void gpu__integral_kernel3(int mu_offset, int mu_steps,	
				      int in_step, int in_steps,
				      int nu_steps,
				      float q_squared_inverse, float r0,
				      float *device__sinb,
				      float *device__sinl,
				      float *device__cosb,
				      float *device__cosl,
				      float *device__V,
				      float *background_integrals, 
				      float *background_correction, 
				      float *stream_integrals, 
				      float *stream_correction) {
	float bg_int, bg_int_correction;
	bg_int = 0.0f;
	bg_int_correction = 0.0f; 

	float *st_int = shared_mem;
	float *st_int_correction = &st_int[blockDim.x * number_streams];
	for (int i = 0; i < number_streams; i++) {
	  st_int[i * blockDim.x + threadIdx.x] = 0.0f;
	  st_int_correction[i * blockDim.x + threadIdx.x] = 0.0f;
	}

	float sinb = device__sinb[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
	float sinl = device__sinl[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
	float cosb = device__cosb[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
	float cosl = device__cosl[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];

	double cosb_x_cosl = cosb * cosl;
	double cosb_x_sinl = cosb * sinl;

	for (int i = 0; i < convolve; i++) {
	  float xyz0, xyz1, xyz2;
	  float rs, rg;
	  xyz2 = tex2D(tex_r_point, i, in_step) * sinb;
	  xyz0 = tex2D(tex_r_point,i,in_step) * 
	    cosb_x_cosl - f_lbr_r;
	  xyz1 = tex2D(tex_r_point,i,in_step) * 
	    cosb_x_sinl;
		  
	  //__fdividef providers faster fp division, with restrictions on
	  // the fact that (q*q) < 2^126 (appendix b.2.1 in nvidia programming guide)
	  rg = sqrtf(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_squared_inverse);
	  rs = rg + r0;
		  
	  float corrected_next_term = tex2D(tex_qw_r3_N, i, in_step) /
	    (rg * rs * rs * rs) - bg_int_correction;
	  float new_sum = bg_int + corrected_next_term;
	  bg_int_correction = (new_sum - bg_int) - corrected_next_term;
	  bg_int = new_sum;
	  
	  for (int j = 0; j < number_streams; j++) {
	    float dotted, sxyz0, sxyz1, sxyz2;
	    sxyz0 = xyz0 - constant__fstream_c[j*3 + 0];
	    sxyz1 = xyz1 - constant__fstream_c[j*3 + 1];
	    sxyz2 = xyz2 - constant__fstream_c[j*3 + 2];
	    
	    dotted = constant__fstream_a[j*3 + 0] * sxyz0 
	      + constant__fstream_a[j*3 + 1] * sxyz1
	      + constant__fstream_a[j*3 + 2] * sxyz2;
	    
	    sxyz0 -= dotted * constant__fstream_a[j*3 + 0];
	    sxyz1 -= dotted * constant__fstream_a[j*3 + 1];
	    sxyz2 -= dotted * constant__fstream_a[j*3 + 2];
      
	    float xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	    float corrected_next_term = 
	      (tex2D(tex_qw_r3_N,i,in_step) 
	       * exp(-(xyz_norm) * 
		     constant__inverse_fstream_sigma_sq2[j])) - st_int_correction[j * blockDim.x + threadIdx.x];
	    float new_sum = st_int[j * blockDim.x + threadIdx.x] + corrected_next_term;
	    st_int_correction[j * blockDim.x + threadIdx.x] = 
	      (new_sum - st_int[j * blockDim.x + threadIdx.x]) - corrected_next_term;
	    st_int[j * blockDim.x + threadIdx.x] = new_sum;
	  }
	}
	
	//define V down here so that one to reduce the number of registers, because a register
	//will be reused
	int nu_step = (threadIdx.x + (blockDim.x * (blockIdx.x + mu_offset))) % nu_steps;
	float V = device__V[nu_step + (in_step * nu_steps)];
	int pos = threadIdx.x + (blockDim.x * (blockIdx.x + mu_offset));
	
	float corrected_next_term = (bg_int * V) - background_correction[pos];
	float new_sum = background_integrals[pos] + corrected_next_term;	
	background_correction[pos] = (new_sum - background_integrals[pos]) - corrected_next_term;
	background_integrals[pos] = new_sum;

	for (int i = 0; i < number_streams; i++) {
	  corrected_next_term = (st_int[i * blockDim.x + threadIdx.x] * V) - stream_correction[pos];
	  new_sum = stream_integrals[pos] + corrected_next_term;
	  stream_correction[pos] = (new_sum - stream_integrals[pos]) - corrected_next_term;
	  stream_integrals[pos] = new_sum;
	  pos += (nu_steps * mu_steps);
	}
}
