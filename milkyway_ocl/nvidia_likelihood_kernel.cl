/*
Copyright 2010 Anthony Waters, Travis Desell,
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

#pragma OPENCL EXTENSION cl_khr_fp64: enable

#define ABSM (4.2)
#define SIGMOID_CURVE_0 (0.9402)
#define SIGMOID_CURVE_1 (1.6171)
#define SIGMOID_CURVE_2 (23.5877)
#define STDEV (0.6)
#define XR (3.0 * STDEV)
#define LBR_R (8.5)

__kernel void
likelihood_kernel(const int convolve,
		  const int number_streams,
		  const double q_sq_inv,
		  const double r0,
		  const double coeff,
		  __global const double *bg_weight,
		  __global const double *st_weight,
		  __global const double *fstream_c,
		  __global const double *fstream_a,
		  __global const double *inv_fstream_sigma_sq2,
		  __global const double *dx,
		  __global const double *qgaus_W,
		  __global const double *stars,
		  __global double *g_probability,
		  __local double *st_int,
		  __global const double *A)
{
  double sinb = stars[get_global_id(0)];
  double sinl = stars[get_global_id(0) + get_global_size(0)];
  double cosb = stars[get_global_id(0) + get_global_size(0)*2];
  double cosl = stars[get_global_id(0) + get_global_size(0)*3];
  double coords = stars[get_global_id(0) + get_global_size(0)*4];
  double log10coords = stars[get_global_id(0) + get_global_size(0)*5];

  //double gPrime = 5.0 * (log10(coords * 1000.0) - 1.0) + ABSM;
  //replace log10(coords * 1000.0) with log10(coords) + 3 - 1
  double gPrime = 5.0 * (log10coords + 2.0) + ABSM;
  double exponent = exp(SIGMOID_CURVE_1 * (gPrime - SIGMOID_CURVE_2));
  double reff_value = SIGMOID_CURVE_0 / (exponent + 1);
  double rPrime3 = coords * coords * coords;

  double reff_xr_rp3 = reff_value * XR / rPrime3;

  double bg_int = 0.0;

  for(int i = 0;i<number_streams;++i)
    st_int[(i * get_local_size(0)) + get_local_id(0)] = 0.0;

  for(int i = 0;i<convolve;++i)
    {
      double g = gPrime + dx[i];
      double r_point = pow(10.0, (g - ABSM) * 0.2 + 1.0) * .001;
      double rPrime3 = r_point * r_point * r_point;
      double qw_r3_N = qgaus_W[i] * rPrime3 * coeff *
	exp(-((g - gPrime) * (g - gPrime) / (2 * STDEV * STDEV)));

      double xyz2 = r_point * sinb;
      double xyz0 = r_point * cosb * cosl - LBR_R;
      double xyz1 = r_point * cosb * sinl;

      //double rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_sq_inv);
      double rg = (xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_sq_inv);
      {
	double x  = (double) rsqrt((float) rg);
      	//fsqrtd
      	x = x * (3.0 - rg*(x*x));
      	double res = x * rg;
      	rg = res * (0.75 - 0.0625*(res*x));
	//
      }
      double rs = rg + r0;

      bg_int += qw_r3_N / (rg * rs * rs *rs);

      for(int j = 0;j < number_streams;++j)
	{
	  double sxyz0 = xyz0 - fstream_c[(j * 3) + 0];
	  double sxyz1 = xyz1 - fstream_c[(j * 3) + 1];
	  double sxyz2 = xyz2 - fstream_c[(j * 3) + 2];

	  double dotted = fstream_a[(j * 3) + 0] * sxyz0
	    + fstream_a[(j * 3) + 1] * sxyz1
	    + fstream_a[(j * 3) + 2] * sxyz2;

	  sxyz0 -= dotted * fstream_a[(j * 3) + 0];
	  sxyz1 -= dotted * fstream_a[(j * 3) + 1];
	  sxyz2 -= dotted * fstream_a[(j * 3) + 2];

	  double xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1)
	    + (sxyz2 * sxyz2);
	  double result = qw_r3_N *
	    exp(-(xyz_norm) * inv_fstream_sigma_sq2[j]);
	  st_int[(j * get_local_size(0)) + get_local_id(0)] += result;
	}
    }
  double probability = bg_int * *bg_weight;
  for(int i = 0;i<number_streams;++i)
    {
      probability +=
	st_int[(i * get_local_size(0)) + get_local_id(0)] * st_weight[i];
    }
  probability *= reff_xr_rp3;

  if (probability == 0.0)
    probability = -238.0;
  else
    //probability = log10(probability);
    probability = (probability);

  g_probability[get_global_id(0)] = probability;
}
