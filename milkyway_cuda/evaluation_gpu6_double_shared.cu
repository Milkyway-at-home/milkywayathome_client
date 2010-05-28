
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

texture<int2, 2, cudaReadModeElementType> tex_r_point;
texture<int2, 2, cudaReadModeElementType> tex_qw_r3_N;
texture<int2, 2, cudaReadModeElementType> tex_r_in_mag;
texture<int2, 2, cudaReadModeElementType> tex_r_in_mag2;
texture<int2, 2, cudaReadModeElementType> tex_fstream_a;
texture<int2, 2, cudaReadModeElementType> tex_fstream_c;
texture<int2, 2, cudaReadModeElementType> tex_fstream_sigma_sq2;

cudaArray **cu_r_point_arrays;
cudaArray **cu_qw_r3_N_arrays;
cudaArray **cu_r_in_mag_arrays;
cudaArray **cu_r_in_mag2_arrays;

static __inline__ __device__
double tex3D_double(texture<int2, 3> tex, int x, int y, int z)
{
  int2 val = tex3D(tex, x, y, z);
  return __hiloint2double(val.x, val.y);
}

static __inline__ __device__
double tex2D_double(texture<int2, 2> tex, int x, int y)
{
  int2 val = tex2D(tex, x, y);
  return __hiloint2double(val.x, val.y);
}

//used to get the hi/lo parts of a double
//from NVIDIA CUDA SDK
union cvt {
  double d;
  signed int i[2];
};

//converts an array of doubles
//into an array of int2s for the use
//in textures
int2*
convert(double* data, int size)
{
  int2 *ret = (int2 *) malloc(sizeof(int2) * size);
  for(unsigned int idx = 0;
      idx < size;++idx)
    {
      cvt cvt;
      cvt.d = data[idx];
      ret[idx].x = cvt.i[1];
      ret[idx].y = cvt.i[0];
    }
  return ret;
}

/**
   Similar to setup_texture except it deals with 2d textures
   that were previously in constant memory
 */
void setup_constant_textures(double *fstream_a, double *fstream_c,
			     double *fstream_sigma_sq2, int number_streams)
{
  // allocate array
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
  cudaArray* cu_array_a;
  cutilSafeCall(cudaMallocArray(&cu_array_a, &channelDesc, 3, number_streams));
  int2 *fstream_a_data = convert(fstream_a, number_streams * 3);
  cutilSafeCall(cudaMemcpyToArray(cu_array_a, 0, 0, fstream_a_data, 3*number_streams * sizeof(int2), cudaMemcpyHostToDevice));

  cudaArray* cu_array_c;
  cutilSafeCall(cudaMallocArray(&cu_array_c, &channelDesc, 3, number_streams));
  int2 *fstream_c_data = convert(fstream_c, number_streams * 3);
  cutilSafeCall(cudaMemcpyToArray(cu_array_c, 0, 0, fstream_c_data, 3*number_streams  * sizeof(int2), cudaMemcpyHostToDevice));

  cudaArray* cu_array_sq2;
  cutilSafeCall(cudaMallocArray(&cu_array_sq2, &channelDesc, 1, number_streams));
  int2 *fstream_sigma_sq2_data = convert(fstream_sigma_sq2, number_streams);
  cutilSafeCall(cudaMemcpyToArray(cu_array_sq2, 0, 0, fstream_sigma_sq2_data, number_streams  * sizeof(int2), cudaMemcpyHostToDevice));

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

  free(fstream_a_data);
  free(fstream_c_data);
  free(fstream_sigma_sq2_data);
}

void setup_r_point_texture(int r_steps, int convolve, int current_integral, double **r_point)
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
  cudaArray* cu_array;
  double *r_point_flat = (double*) malloc(sizeof(double) * r_steps * convolve);
  int i,j;
  for(i = 0;i<r_steps;++i)
    {
      for(j = 0;j<convolve;++j)
	{
	  r_point_flat[i * convolve + j] =
	    r_point[i][j];
	}
    }
  int2 *host_int2 = convert(r_point_flat, r_steps * convolve);
  cutilSafeCall(cudaMallocArray(&cu_array, &channelDesc, convolve, r_steps));
  cutilSafeCall(cudaMemcpyToArray(cu_array, 0, 0, host_int2, r_steps * convolve * sizeof(int2), cudaMemcpyHostToDevice));
  cu_r_point_arrays[current_integral] = cu_array;
  free(r_point_flat);
  free(host_int2);
}

void setup_qw_r3_N_texture(int r_steps, int convolve, int current_integral, double **qw_r3_N)
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
  cudaArray* cu_array;
  double *qw_r3_N_flat = (double*) malloc(sizeof(double) * r_steps * convolve);
  int i,j;
  for(i = 0;i<r_steps;++i)
    {
      for(j = 0;j<convolve;++j)
	{
	  qw_r3_N_flat[i * convolve + j] =
	    qw_r3_N[i][j];
	}
    }
  int2 *host_int2 = convert(qw_r3_N_flat, r_steps * convolve);
  cutilSafeCall(cudaMallocArray(&cu_array, &channelDesc, convolve, r_steps));
  cutilSafeCall(cudaMemcpyToArray(cu_array, 0, 0, host_int2, r_steps * convolve * sizeof(int2), cudaMemcpyHostToDevice));
  cu_qw_r3_N_arrays[current_integral] = cu_array;
  free(qw_r3_N_flat);
  free(host_int2);
}

void setup_r_in_mag_texture2(int r_steps, int convolve, int current_integral,
			     double **r_in_mag, double **r_in_mag2)
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
  cudaArray* cu_array;
  cudaArray* cu_array2;
  double *r_in_mag_flat = (double*) malloc(sizeof(double) * r_steps * convolve);
  double *r_in_mag2_flat = (double*) malloc(sizeof(double) * r_steps * convolve);
  int i,j;
  for(i = 0;i<r_steps;++i)
    {
      for(j = 0;j<convolve;++j)
	{
	  r_in_mag_flat[i * convolve + j] =
	    r_in_mag[i][j];
	  r_in_mag2_flat[i * convolve + j] =
	    r_in_mag2[i][j];
	}
    }
  int2 *r_in_mag_int2 = convert(r_in_mag_flat, r_steps * convolve);
  int2 *r_in_mag2_int2 = convert(r_in_mag2_flat, r_steps * convolve);
  cutilSafeCall(cudaMallocArray(&cu_array, &channelDesc, convolve, r_steps));
  cutilSafeCall(cudaMemcpyToArray(cu_array, 0, 0, r_in_mag_int2, r_steps * convolve * sizeof(int2), cudaMemcpyHostToDevice));
  cu_r_in_mag_arrays[current_integral] = cu_array;
  cutilSafeCall(cudaMallocArray(&cu_array2, &channelDesc, convolve, r_steps));
  cutilSafeCall(cudaMemcpyToArray(cu_array2, 0, 0, r_in_mag2_int2, r_steps * convolve * sizeof(int2), cudaMemcpyHostToDevice));
  cu_r_in_mag2_arrays[current_integral] = cu_array2;
  free(r_in_mag_flat);
  free(r_in_mag_int2);
  free(r_in_mag2_flat);
  free(r_in_mag2_int2);
}

/**
   Allocates cu arrays for the tex_device_lb texture and sets up
   parts of the texture
*/
void allocate_cu_arrays(int number_integrals) {
  cu_r_point_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);
  cu_qw_r3_N_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);
  cu_r_in_mag_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);
  cu_r_in_mag2_arrays = (cudaArray**) malloc(sizeof(cudaArray*) * number_integrals);

  tex_r_point.addressMode[0] = cudaAddressModeClamp;
  tex_r_point.addressMode[1] = cudaAddressModeClamp;
  tex_r_point.filterMode = cudaFilterModePoint;
  tex_r_point.normalized = false;

  tex_qw_r3_N.addressMode[0] = cudaAddressModeClamp;
  tex_qw_r3_N.addressMode[1] = cudaAddressModeClamp;
  tex_qw_r3_N.filterMode = cudaFilterModePoint;
  tex_qw_r3_N.normalized = false;

  tex_r_in_mag.addressMode[0] = cudaAddressModeClamp;
  tex_r_in_mag.addressMode[1] = cudaAddressModeClamp;
  tex_r_in_mag.filterMode = cudaFilterModePoint;
  tex_r_in_mag.normalized = false;

  tex_r_in_mag2.addressMode[0] = cudaAddressModeClamp;
  tex_r_in_mag2.addressMode[1] = cudaAddressModeClamp;
  tex_r_in_mag2.filterMode = cudaFilterModePoint;
  tex_r_in_mag2.normalized = false;
}

void bind_texture(int current_integral) {
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
  //  printf("binding the tex_device_lb texture for integral %u\n",
  //	 current_integral);
  //cutilSafeCall(cudaBindTextureToArray(tex_device_lb, cu_arrays[current_integral], channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_r_point, cu_r_point_arrays[current_integral], channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_qw_r3_N, cu_qw_r3_N_arrays[current_integral], channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_r_in_mag, cu_r_in_mag_arrays[current_integral], channelDesc));
  cutilSafeCall(cudaBindTextureToArray(tex_r_in_mag2, cu_r_in_mag2_arrays[current_integral], channelDesc));
}

// 11 DP ops (including 2 conversion) + 1 SP reciprocal, guess 10 DP flops would be an adequate count
__device__ double divd(double a, double b)	// accurate to 1 ulp, i.e the last bit of the double precision number
{	 // cuts some corners on the numbers range but is significantly faster, employs "faithful rounding"
  double r;
  double y;
  double c;
  y = (double)(1.0f/__double2float_rn(b));	// 22bit estimate of reciprocal, limits range to float range, but spares the exponent extraction
  c = 1.0 - b*y;
  y += y*c;	 // first Newton iteration => 44bit accurate
  r = a*y;	 // second iteration works directly on a/b, so it is effectively
  c = a - b*r;	 // a Newton-Markstein iteration without guaranteed round to nearest
  return(r + y*c);// should be generally accurate to 1 ulp, i.e "faithful rounding" (or close to it, depending on definition)
}// on a GT200 should be round to nearest most of the time, albeit not guaranteed
// one would have to add a second Newton iteration before the Markstein rounding step,
// but one looses the gained half bit of precision in the following additions, so the added effort doesn't make sense

// 11 DP ops including 2 conversions + 1 SP rsqrt, a count of 10 flops may be adequate
__device__ double fsqrtd(double y)	// accurate to 1 ulp, i.e the last bit of the double precision number
{	// cuts some corners on the numbers range but is significantly faster, employs "faithful rounding"
  double x, res;
  x = (double)rsqrtf((float)y);	// 22bit estimate for reciprocal square root, limits range to float range, but spares the exponent extraction
  x = x * (3.0 - y*(x*x));	// first Newton iteration (44bit accurate)
  res = x * y;	// do final iteration directly on sqrt(y) and not on the inverse
  return(res * (0.75 - 0.0625*(res*x)));
}	// same precision as division (1 ulp)


template <unsigned int number_streams, unsigned int convolve,
	  unsigned int aux_bg_profile>
__global__ void gpu__integral_kernel3(int mu_offset, int mu_steps,
				      int in_step, int in_steps,
				      int nu_steps, int total_threads,
				      double bg_a, double bg_b, double bg_c,
				      double q_squared_inverse, double r0,
				      double *device__sinb,
				      double *device__sinl,
				      double *device__cosb,
				      double *device__cosl,
				      double *device__V,
				      double *background_integrals,
				      double *background_integrals_c,
				      double *stream_integrals,
				      double *stream_integrals_c,
				      double *fstream_c,
				      double *fstream_a)
{
  double *s_fstream_c = shared_mem;
  double *s_fstream_a = &s_fstream_c[number_streams * 3];

  if (threadIdx.x < number_streams * 3)
    {
      s_fstream_a[threadIdx.x] = fstream_a[threadIdx.x];
      s_fstream_c[threadIdx.x] = fstream_c[threadIdx.x];
    }

  double bg_int = 0.0;
  double st_int0 = 0.0;
  double st_int1 = 0.0;
  double st_int2 = 0.0;
  double st_int3 = 0.0;

  double sinb = device__sinb[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  double sinl = device__sinl[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  double cosb = device__cosb[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  double cosl = device__cosl[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];

  double cosb_x_cosl = cosb * cosl;
  double cosb_x_sinl = cosb * sinl;

  for (int i = 0; i < convolve; i++) {
    double r_point = tex2D_double(tex_r_point,i,in_step);
    double xyz2 = r_point * sinb;
    double xyz0 = r_point * cosb_x_cosl - d_lbr_r;
    double xyz1 = r_point * cosb_x_sinl;

    double qw_r3_N = tex2D_double(tex_qw_r3_N,i,in_step);
    {
      double rg = fsqrtd(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2)
			 * q_squared_inverse);
      double rs = rg + r0;

      if (aux_bg_profile == 1)
	{
	  double r_in_mag = tex2D_double(tex_r_in_mag, i, in_step);
	  double r_in_mag2 = tex2D_double(tex_r_in_mag2, i, in_step);
	  double h_prob = divd(qw_r3_N , (rg * rs * rs * rs));
	  double aux_prob = qw_r3_N * ( bg_a * r_in_mag2 + bg_b * r_in_mag + bg_c );
	  bg_int += h_prob + aux_prob;
	}
      else
	{
	  bg_int += divd(qw_r3_N , (rg * rs * rs * rs));
	}
    }
    if (number_streams >= 1)
      {
	//stream 0
	double sxyz0 = xyz0 - s_fstream_c[0];
	double sxyz1 = xyz1 - s_fstream_c[1];
	double sxyz2 = xyz2 - s_fstream_c[2];

	double dotted = s_fstream_a[0] * sxyz0
	  + s_fstream_a[1] * sxyz1
	  + s_fstream_a[2] * sxyz2;

	sxyz0 -= dotted * s_fstream_a[0];
	sxyz1 -= dotted * s_fstream_a[1];
	sxyz2 -= dotted * s_fstream_a[2];

	double xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	double result = (qw_r3_N
			 * exp(-(xyz_norm) *
			       constant_inverse_fstream_sigma_sq2[0]));
	st_int0 += result;
      }
    if (number_streams >= 2)
      {
	//stream 1
	double sxyz0 = xyz0 - s_fstream_c[3];
	double sxyz1 = xyz1 - s_fstream_c[4];
	double sxyz2 = xyz2 - s_fstream_c[5];

	double dotted = s_fstream_a[3] * sxyz0
	  + s_fstream_a[4] * sxyz1
	  + s_fstream_a[5] * sxyz2;

	sxyz0 -= dotted * s_fstream_a[3];
	sxyz1 -= dotted * s_fstream_a[4];
	sxyz2 -= dotted * s_fstream_a[5];

	double xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	double result = (qw_r3_N
			 * exp(-(xyz_norm) *
			       constant_inverse_fstream_sigma_sq2[1]));
	st_int1 += result;
      }
    if (number_streams >= 3)
      {
	//stream 2
	double sxyz0 = xyz0 - s_fstream_c[6];
	double sxyz1 = xyz1 - s_fstream_c[7];
	double sxyz2 = xyz2 - s_fstream_c[8];

	double dotted = s_fstream_a[6] * sxyz0
	  + s_fstream_a[7] * sxyz1
	  + s_fstream_a[8] * sxyz2;

	sxyz0 -= dotted * s_fstream_a[6];
	sxyz1 -= dotted * s_fstream_a[7];
	sxyz2 -= dotted * s_fstream_a[8];

	double xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	double result = (qw_r3_N
			 * exp(-(xyz_norm) *
			       constant_inverse_fstream_sigma_sq2[2]));
	st_int2 += result;
      }
    if (number_streams >= 4)
      {
	//stream 3
	double sxyz0 = xyz0 - s_fstream_c[9];
	double sxyz1 = xyz1 - s_fstream_c[10];
	double sxyz2 = xyz2 - s_fstream_c[11];

	double dotted = s_fstream_a[9] * sxyz0
	  + s_fstream_a[10] * sxyz1
	  + s_fstream_a[11] * sxyz2;

	sxyz0 -= dotted * s_fstream_a[9];
	sxyz1 -= dotted * s_fstream_a[10];
	sxyz2 -= dotted * s_fstream_a[11];

	double xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
	double result = (qw_r3_N
			 * exp(-(xyz_norm) *
			       constant_inverse_fstream_sigma_sq2[3]));
	st_int3 += result;
      }
  }

  //define V down here so that one to reduce the number of registers,
  //because a register will be reused
  int nu_step = (threadIdx.x + (blockDim.x * (blockIdx.x + mu_offset)))
    % nu_steps;
  double V = device__V[nu_step + (in_step * nu_steps)];
  int pos = threadIdx.x + (blockDim.x * (blockIdx.x + mu_offset));
  double temp = background_integrals[pos];
  bg_int = bg_int * V;
  background_integrals[pos] += (bg_int);
  background_integrals_c[pos] += ((bg_int) -
				  (background_integrals[pos] - temp));
  if (number_streams >= 1)
    {
      st_int0 *= V;
      double temp = stream_integrals[pos];
      stream_integrals[pos] += st_int0;
      stream_integrals_c[pos] += (st_int0 - (stream_integrals[pos] - temp));
      pos += total_threads;
    }
  if (number_streams >= 2)
    {
      st_int1 *= V;
      double temp = stream_integrals[pos];
      stream_integrals[pos] += st_int1;
      stream_integrals_c[pos] += (st_int1 - (stream_integrals[pos] - temp));
      pos += total_threads;
    }
  if (number_streams >= 3)
    {
      st_int2 *= V;
      double temp = stream_integrals[pos];
      stream_integrals[pos] += st_int2;
      stream_integrals_c[pos] += (st_int2 - (stream_integrals[pos] - temp));
      pos += total_threads;
    }
  if (number_streams >= 4)
    {
      st_int3 *= V;
      double temp = stream_integrals[pos];
      stream_integrals[pos] += st_int3;
      stream_integrals_c[pos] += (st_int3 - (stream_integrals[pos] - temp));
    }
}

template <unsigned int number_streams, unsigned int convolve>
__global__ void gpu__integral_kernel3_aux(int mu_offset, int mu_steps,
					  int in_step, int in_steps,
					  int nu_steps, int total_threads,
					  double q_squared_inverse, double r0,
					  double bg_a, double bg_b, double bg_c,
					  double *device__sinb,
					  double *device__sinl,
					  double *device__cosb,
					  double *device__cosl,
					  double *device__V,
					  double *background_integrals,
					  double *stream_integrals) {
  double *st_int = shared_mem;
  double bg_int = 0.0;
  for (int i = 0; i < number_streams; i++) {
    st_int[i * blockDim.x + threadIdx.x] = 0.0;
  }

  double sinb = device__sinb[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  double sinl = device__sinl[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  double cosb = device__cosb[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];
  double cosl = device__cosl[threadIdx.x + ((mu_offset + blockIdx.x) * blockDim.x)];

  double cosb_x_cosl = cosb * cosl;
  double cosb_x_sinl = cosb * sinl;

  for (int i = 0; i < convolve; i++) {
    double xyz0, xyz1, xyz2;
    double rs, rg;
    xyz2 =  tex2D_double(tex_r_point,i,in_step) *
      sinb;

    xyz0 = tex2D_double(tex_r_point,i,in_step) *
      cosb_x_cosl - d_lbr_r;
    xyz1 = tex2D_double(tex_r_point,i,in_step) *
    cosb_x_sinl;

    rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_squared_inverse);
    rs = rg + r0;

    double r_in_mag = 4.2 + 5 * (log10(1000 * tex2D_double(tex_r_point,i,in_step)) - 1);
    double h_prob = tex2D_double(tex_qw_r3_N,i,in_step) / (rg * rs * rs * rs);
    double aux_prob = tex2D_double(tex_qw_r3_N,i,in_step) *
      (bg_a * r_in_mag * r_in_mag + bg_b * r_in_mag + bg_c);
    bg_int += h_prob + aux_prob;

    for (int j = 0; j < number_streams; j++) {
      double dotted, sxyz0, sxyz1, sxyz2;
      sxyz0 = xyz0 - tex2D_double(tex_fstream_c, 0, j);
      sxyz1 = xyz1 - tex2D_double(tex_fstream_c, 1, j);
      sxyz2 = xyz2 - tex2D_double(tex_fstream_c, 2, j);

      dotted = tex2D_double(tex_fstream_a, 0,j) * sxyz0
      	+ tex2D_double(tex_fstream_a, 1, j) * sxyz1
      	+ tex2D_double(tex_fstream_a, 2, j) * sxyz2;

      sxyz0 -= dotted * tex2D_double(tex_fstream_a, 0, j);
      sxyz1 -= dotted * tex2D_double(tex_fstream_a, 1, j);
      sxyz2 -= dotted * tex2D_double(tex_fstream_a, 2, j);

      double xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2);
      double result = (tex2D_double(tex_qw_r3_N,i,in_step)
      	       * exp(-(xyz_norm) * constant_inverse_fstream_sigma_sq2[j]));
      st_int[j * blockDim.x + threadIdx.x] += result;
    }
  }

  //define V down here so that one to reduce the number of registers, because a register
  //will be reused
  int nu_step = (threadIdx.x + (blockDim.x * (blockIdx.x + mu_offset))) % nu_steps;
  double V = device__V[nu_step + (in_step * nu_steps)];
  int pos = threadIdx.x + (blockDim.x * (blockIdx.x + mu_offset));
  background_integrals[pos] += (bg_int * V);
  for (int i = 0; i < number_streams; i++) {
    stream_integrals[pos] += st_int[i * blockDim.x + threadIdx.x] * V;
    pos += (nu_steps * mu_steps);
  }
}
