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

texture<float, 2, cudaReadModeElementType> tex_r_point;
texture<float, 2, cudaReadModeElementType> tex_qw_r3_N;
texture<float, 2, cudaReadModeElementType> tex_fstream_a;
texture<float, 2, cudaReadModeElementType> tex_fstream_c;
texture<float, 2, cudaReadModeElementType> tex_fstream_sigma_sq2;

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
  cutilSafeCall(cudaMemcpyToArray(cu_array_a, 0, 0, fstream_a, 3*number_streams* sizeof(float), cudaMemcpyHostToDevice));

  cudaArray* cu_array_c;
  cutilSafeCall(cudaMallocArray(&cu_array_c, &channelDesc, 3, number_streams)); 
  cutilSafeCall(cudaMemcpyToArray(cu_array_c, 0, 0, fstream_c, 3*number_streams * sizeof(float), cudaMemcpyHostToDevice));

  cudaArray* cu_array_sq2;
  cutilSafeCall(cudaMallocArray(&cu_array_sq2, &channelDesc, 2, number_streams)); 
  cutilSafeCall(cudaMemcpyToArray(cu_array_sq2, 0, 0, fstream_sigma_sq2, 2*number_streams* sizeof(float), cudaMemcpyHostToDevice));
  
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
   Allocates cu arrays for the tex_device_lb texture and sets up
   parts of the texture
*/
void allocate_cu_arrays(int number_integrals) {
  cu_r_point_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);
  cu_qw_r3_N_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);
  //set texture parameters
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
  cutilSafeCall(cudaBindTextureToArray(tex_r_point, cu_r_point_arrays[current_integral], channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_qw_r3_N, cu_qw_r3_N_arrays[current_integral], channelDesc));
}


template <unsigned int number_streams, unsigned int convolve> 
__global__ void gpu__integral_kernel3(int mu_offset, int mu_steps,	
				      int in_step, int in_steps,
				      int nu_steps, int total_threads,
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
  float bg_int = 0.0;
  float st_int0 = 0.0;
  float st_int1 = 0.0;
  float st_int2 = 0.0;
  float st_int3 = 0.0;

  float sinb = device__sinb[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  float sinl = device__sinl[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  float cosb = device__cosb[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  float cosl = device__cosl[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];

  float cosb_x_cosl = cosb * cosl;
  float cosb_x_sinl = cosb * sinl;

  for (int i = 0; i < convolve; i++) {
    float xyz2 = tex2D(tex_r_point,i,in_step) * sinb;
    float xyz0 = tex2D(tex_r_point,i,in_step) * cosb_x_cosl - d_lbr_r;
    float xyz1 = tex2D(tex_r_point,i,in_step) * cosb_x_sinl;

    {
      float rg = __fsqrt_rn(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) 
			 * q_squared_inverse);
      float rs = rg + r0;
      
      bg_int += (tex2D(tex_qw_r3_N,i,in_step) / (rg * rs * rs * rs));
    }
    if (number_streams >= 1)
      {
	//stream 0
	float sxyz0 = xyz0 - constant_fstream_c[0];
	float sxyz1 = xyz1 - constant_fstream_c[1];
	float sxyz2 = xyz2 - constant_fstream_c[2];

	float dotted = constant_fstream_a[0] * sxyz0 
	  + constant_fstream_a[1] * sxyz1
	  + constant_fstream_a[2] * sxyz2;
	
	sxyz0 -= dotted * constant_fstream_a[0];
	sxyz1 -= dotted * constant_fstream_a[1];
	sxyz2 -= dotted * constant_fstream_a[2];
	
	float xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	float result = (tex2D(tex_qw_r3_N,i,in_step) 
			 * exp(-(xyz_norm) * 
			       constant_inverse_fstream_sigma_sq2[0]));   
	st_int0 += result;
      }
    if (number_streams >= 2)
      {
	//stream 1
	float sxyz0 = xyz0 - constant_fstream_c[3];
	float sxyz1 = xyz1 - constant_fstream_c[4];
	float sxyz2 = xyz2 - constant_fstream_c[5];
	
	float dotted = constant_fstream_a[3] * sxyz0 
	  + constant_fstream_a[4] * sxyz1
	  + constant_fstream_a[5] * sxyz2;
	
	sxyz0 -= dotted * constant_fstream_a[3];
	sxyz1 -= dotted * constant_fstream_a[4];
	sxyz2 -= dotted * constant_fstream_a[5];
	  
	float xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	float result = (tex2D(tex_qw_r3_N,i,in_step) 
			 * exp(-(xyz_norm) * 
			       constant_inverse_fstream_sigma_sq2[1]));   
	st_int1 += result;
      }
    if (number_streams >= 3)
      {
	//stream 2
	float sxyz0 = xyz0 - constant_fstream_c[6];
	float sxyz1 = xyz1 - constant_fstream_c[7];
	float sxyz2 = xyz2 - constant_fstream_c[8];
	
	float dotted = constant_fstream_a[6] * sxyz0 
	  + constant_fstream_a[7] * sxyz1
	  + constant_fstream_a[8] * sxyz2;
	
	sxyz0 -= dotted * constant_fstream_a[6];
	sxyz1 -= dotted * constant_fstream_a[7];
	sxyz2 -= dotted * constant_fstream_a[8];
	
	float xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	float result = (tex2D(tex_qw_r3_N,i,in_step) 
			 * exp(-(xyz_norm) * 
			       constant_inverse_fstream_sigma_sq2[2]));   
	st_int2 += result;
      }
    if (number_streams >= 4)
      {
	//stream 3
	float sxyz0 = xyz0 - constant_fstream_c[9];
	float sxyz1 = xyz1 - constant_fstream_c[10];
	float sxyz2 = xyz2 - constant_fstream_c[11];
	
	float dotted = constant_fstream_a[9] * sxyz0 
	  + constant_fstream_a[10] * sxyz1
	  + constant_fstream_a[11] * sxyz2;
	
	sxyz0 -= dotted * constant_fstream_a[9];
	sxyz1 -= dotted * constant_fstream_a[10];
	sxyz2 -= dotted * constant_fstream_a[11];
	
	float xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	float result = (tex2D(tex_qw_r3_N,i,in_step) 
			 * exp(-(xyz_norm) * 
			       constant_inverse_fstream_sigma_sq2[3]));   
	st_int3 += result;
      }
  }

  //define V down here so that one to reduce the number of registers,
  //because a register will be reused
  int nu_step = (threadIdx.x + (blockDim.x * (blockIdx.x + mu_offset)))
    % nu_steps;
  float V = device__V[nu_step + (in_step * nu_steps)];
  int pos = threadIdx.x + (blockDim.x * (blockIdx.x + mu_offset));
  background_integrals[pos] += (bg_int * V);
  if (number_streams >= 1)
    {
      stream_integrals[pos] += st_int0 * V;
      pos += total_threads;
    }
  if (number_streams >= 2)
    {
      stream_integrals[pos] += st_int1 * V;
      pos += total_threads;
    }
  if (number_streams >= 3)
    {
      stream_integrals[pos] += st_int2 * V;
      pos += total_threads;
    }
  if (number_streams >= 4)
    {
      stream_integrals[pos] += st_int3 * V;
    }
}

template <unsigned int number_streams, unsigned int convolve> 
__global__ void gpu__integral_kernel3_aux(int mu_offset, int mu_steps,	
					  int in_step, int in_steps,
					  int nu_steps, int total_threads,
					  float q_squared_inverse, float r0,
					  float bg_a, float bg_b, float bg_c,
					  float *device__sinb,
					  float *device__sinl,
					  float *device__cosb,
					  float *device__cosl,
					  float *device__V,
					  float *background_integrals, 
					  float *background_correction, 
					  float *stream_integrals, 
					  float *stream_correction) {

}
