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

#define LBR_R (8.5)
#include "exp.cl"
#include "sqrt.cl"

__kernel void
integral_kernel(const int convolve,
		const int number_streams,
		const int r_step,
		const int nu_steps,
		const double q_sq_inv,
		const double r0,
		__global const double *r_point,
		__global const double *qw_r3_N,
		__global const double *g_sinb,
		__global const double *g_sinl,
		__global const double *g_cosb,
		__global const double *g_cosl,
		__global const double *g_v,
		__constant double2 *fstream_c,
		__constant double2 *fstream_a,
		__constant double2 *inv_fstream_sigma_sq2,
		__global double *g_bg_int,
		__global double *g_st_int,
		__local double2 *s_fstream_c,
		__local double2 *s_fstream_a,
		__local double2 *s_inv_fstream_sigma_sq2)
{
  double bg_int = 0.0;
  double2 st_int = {0.0, 0.0};

  if (get_local_id(0) == 0)
    {
      s_fstream_c[0] = (fstream_c[0], fstream_c[3] );
      s_fstream_c[1] = (fstream_c[1], fstream_c[4] );
      s_fstream_c[2] = (fstream_c[2], fstream_c[5] );

      s_fstream_a[0] = (fstream_a[0], fstream_a[3] );
      s_fstream_a[1] = (fstream_a[1], fstream_a[4] );
      s_fstream_a[2] = (fstream_a[2], fstream_a[5] );

      s_inv_fstream_sigma_sq2[0] = (inv_fstream_sigma_sq2[0],
				    inv_fstream_sigma_sq2[1]);
    }

  double sinb = g_sinb[get_global_id(0)];
  double sinl = g_sinl[get_global_id(0)];
  double cosb = g_cosb[get_global_id(0)];
  double cosl = g_cosl[get_global_id(0)];


  for(int i = 0;i<convolve;++i)
    {
      double xyz2 = r_point[(r_step * convolve) + i] * sinb;
      double xyz0 = r_point[(r_step * convolve) + i] * cosb *
	cosl - LBR_R;
      double xyz1 = r_point[(r_step * convolve) + i] * cosb * sinl;

      //double rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_sq_inv);
      //fsqrtd
      double rg = xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_sq_inv;
      //
      {
	//manually convert double to float (ATI Hardware can't do it)
	union {double d;int val[2]} test;
	test.d = rg;
	int sign, exponent;
	float value;
	sign = (int) ((test.val[1] & 0x80000000) >> 31);
	exponent = (int) ((test.val[1] & 0x7FF00000) >> 20);
	float mult = ldexp(1.0f, exponent - 1023 - 52);
	int val1 = test.val[1] & 0x000FFFFF;
	int val2 = test.val[0];
	unsigned long long mantissa =
	  (((unsigned long long) val1) << 32) + val2 + 0x0010000000000000;
	value = mult * mantissa;
	if (sign)
	  value = -value;
	//	printf("rg: %.15f sign:%d exp:%d value:%.15f\n",
	//	       rg,sign, exponent, value);
	//rsqrt(y) -> http://en.wikipedia.org/wiki/Fast_inverse_square_root
	double x  = (double) rsqrt(value);
	//double x  = (double) rsqrt((float) rg);
	/* //fsqrtd */
	x = x * (3.0 - rg*(x*x));
	double res = x * rg;
	rg = res * (0.75 - 0.0625*(res*x));
      }
      //
      double rs = rg + r0;

      bg_int += qw_r3_N[(r_step * convolve) + i] / (rg * rs * rs *rs);

      double2 sxyz0 = xyz0 - s_fstream_c[0];
      double2 sxyz1 = xyz1 - s_fstream_c[1];
      double2 sxyz2 = xyz2 - s_fstream_c[2];

      double2 dotted = s_fstream_a[0] * sxyz0
	+ s_fstream_a[1] * sxyz1
	+ s_fstream_a[2] * sxyz2;

      sxyz0 -= dotted * s_fstream_a[0];
      sxyz1 -= dotted * s_fstream_a[1];
      sxyz2 -= dotted * s_fstream_a[2];

      double2 xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1)
	+ (sxyz2 * sxyz2);
      //cephes dp library for exp
      double2 x = -(xyz_norm) * s_inv_fstream_sigma_sq2[0];

      int2 n = (int2)(x * 1.442695040888963407359924
		      + 0.5) - 1;//log2e, estimated floor
      double2 px = ((double) n.x, (double) n.y);

      x -= px * 6.93145751953125E-1; //C1
      x -= px * 1.42860682030941723212E-6; //C2
      double2 xx = x * x;
      px = x * (( 1.26177193074810590878E-4 * xx +
      		  3.02994407707441961300E-2) * xx +
      		9.99999999999999999910E-1); //polevl
      x = px / ((((3.00198505138664455042E-6 * xx +
      		   2.52448340349684104192E-3) * xx +
      		  2.27265548208155028766E-1) * xx +
      		 2.00000000000000000009E0) - px); //polevl(xx,Q,3)
      x = 1.0 + 2.0 * x;

      //printf("x:%.15f n:=%d\n", x, n);
      //ldexp, x * 2^n
      //ldexp(x, n);

      /* if (n < 0) */
      /* 	{ */
      /* 	  if (n < -60) */
      /* 	    x = -0.0; */
      /* 	  else */
      /* 	    { */
      /* 	      n = n * -1; */
      /* 	      if (n > 30) */
      /* 		{ */
      /* 		  x = x / (double) (2 << (30-1)); */
      /* 		  n -= 30; */
      /* 		} */
      /* 	      x = x / (double) (2 << (n-1)); */
      /* 	    } */
      /* 	} */
      /* else if (n > 0) */
      /* 	{ */
      /* 	  x = x * (double) (2 << (n-1)); */
      /* 	} */

      //printf("x:%.15f n:=%d\n", x, n);
      //double result = qw_r3_N[(r_step * convolve) + i] *
      //exp(-(xyz_norm) * inv_fstream_sigma_sq2[j]);
      double2 result = qw_r3_N[(r_step * convolve) + i] * x;
      st_int += result;
    }
  int nu_step = get_global_id(0) % nu_steps;
  double v = g_v[nu_step + (r_step * nu_steps)];
  g_bg_int[get_global_id(0)] += (bg_int * v);
  g_st_int[get_global_id(0)] += (st_int.x * v);
  g_st_int[get_global_id(0) + get_global_size(0)] += (st_int.y * v);
}
