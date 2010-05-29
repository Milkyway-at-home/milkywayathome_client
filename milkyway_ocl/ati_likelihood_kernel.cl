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
  //double exponent = exp(SIGMOID_CURVE_1 * (gPrime - SIGMOID_CURVE_2));
  double exponent = (SIGMOID_CURVE_1 * (gPrime - SIGMOID_CURVE_2));
  {
    //exp from cephes dp library
    double px = (double) (int)(exponent * 1.4426950408889634073599246
  			     + 0.5) - 1; //log2e, estimated floor
    int n = px;
    exponent -= px * 6.93145751953125E-1; //C1
    exponent -= px * 1.42860682030941723212E-6; //C2
    double xx = exponent * exponent;
    px = exponent * (( 1.26177193074810590878E-4 * xx +
  		       3.02994407707441961300E-2) * xx +
  		     9.99999999999999999910E-1); //polevl
    exponent = px / ((((3.00198505138664455042E-6 * xx +
  			2.52448340349684104192E-3) * xx +
  		       2.27265548208155028766E-1) * xx +
  		      2.00000000000000000009E0) - px); //polevl(xx,Q,3)
    exponent = 1.0 + 2.0 * exponent;
    //ldexp, x * 2^n
    if (n < 0)
      {
  	if (n < -60)
  	  exponent = -0.0;
  	else
  	  {
  	    n = n * -1;
  	    if (n > 30)
  	      {
  		exponent = exponent / (double) (2 << (30-1));
  		n -= 30;
  	      }
  	    exponent = exponent / (double) (2 << (n-1));
  	  }
      }
    else if (n > 0)
      {
  	exponent = exponent * (double) (2 << (n-1));
      }
  }
  //
  double reff_value = SIGMOID_CURVE_0 / (exponent + 1);
  double rPrime3 = coords * coords * coords;

  double reff_xr_rp3 = reff_value * XR / rPrime3;

  double bg_int = 0.0;

  for(int i = 0;i<number_streams;++i)
    st_int[(i * get_local_size(0)) + get_local_id(0)] = 0.0;

  for(int i = 0;i<convolve;++i)
    {
      double g = gPrime + dx[i];
      //double r_point = pow(10.0, (g - ABSM) * 0.2 + 1.0) * .001;
      double z = 0.0;
      //
      {
	double y = (g - ABSM) * 0.2 + 1.0;
      	int e = 4;
      	double x = 0.625;
      	int i = 10;

      	x -= 6.48419777325504820276E-1; //A[i];
      	x -= 1.26912513974441574796E-17; //B[i/2];
      	x /= 6.48419777325504820276E-1; //A[i];

      	z = x * x;
      	double w = x * (z * (((((4.97778295871696322025E-1 * x +
      				 3.73336776063286838734E0) * x) +
      			       7.69994162726912503298E0) * x) +
      			     4.66651806774358464979E0) /
      			((((x + 9.33340916416696166113E0) * x +
      			   2.79999886606328401649E1) * x +
      			  3.35994905342304405431E1) * x +
      			 1.39995542032307539578E1));
      	w = w - z * 0.5;

      	w = w + 0.44269504088896340736 * w; //log(2)ea
      	z = w + 0.44269504088896340736 * x;
      	z = z + x;

      	w = -i;
      	w = w / 16.0;
      	w += e;

      	double ya = y * 16;
      	ya = (double) (int) (ya);
      	ya = ya / 16.0;
      	double yb = y - ya;

      	double W = z * y + w * yb;
      	double Wa = W * 16;
      	Wa = (double) (int) (Wa);
      	Wa = Wa / 16.0;
      	double Wb = W - Wa;

      	W = Wa + w * ya;
      	Wa = W * 16;
      	Wa = (double) (int) (Wa);
      	Wa = Wa / 16.0;
      	double u = W - Wa;

      	W = Wb + u;
      	Wb = W * 16;
      	Wb = (double) (int) (Wb);
      	Wb = Wb / 16.0;
      	w = (Wa + Wb) * 16.0;

      	e = w;
      	Wb = W - Wb;
      	if (Wb > 0.0)
      	  {
      	    e += 1;
      	    Wb -= 0.0625;
      	  }

      	z = Wb * ((((((1.49664108433729301083E-5 * Wb +
      		       1.54010762792771901396E-4) * Wb +
      		      1.33335476964097721140E-3) * Wb +
      		     9.61812908476554225149E-3) * Wb +
      		    5.55041086645832347466E-2) * Wb +
      		   2.40226506959099779976E-1) * Wb +
      		  6.93147180559945308821E-1);

      	if (e < 0)
      	  i = 0;
      	else
      	  i = 1;
      	i = e/16 + i;
      	e = 16 * i - e;
      	w = A[e];
      	z = w + w * z;
      	//z = ldexp(z, i);
      	z = z * (2 << (i-1));
      }
      double r_point = z * .001;
      //
      double rPrime3 = r_point * r_point * r_point;
      double qwx = -((g - gPrime) * (g - gPrime) / (2 * STDEV * STDEV));
      {
      	//cephes dp library for exp
      	double px = (double) (int)(qwx * 1.4426950408889634073599246
      				   + 0.5) - 1; //log2e, estimated floor
      	int n = px;
      	qwx -= px * 6.93145751953125E-1; //C1
      	qwx -= px * 1.42860682030941723212E-6; //C2
      	double xx = qwx * qwx;
      	px = qwx * (( 1.26177193074810590878E-4 * xx +
      		      3.02994407707441961300E-2) * xx +
      		    9.99999999999999999910E-1); //polevl
      	qwx = px / ((((3.00198505138664455042E-6 * xx +
      		       2.52448340349684104192E-3) * xx +
      		      2.27265548208155028766E-1) * xx +
      		     2.00000000000000000009E0) - px); //polevl(xx,Q,3)
      	qwx = 1.0 + 2.0 * qwx;
      	//ldexp, x * 2^n
      	if (n < 0)
      	  {
      	    if (n < -60)
      	      qwx = -0.0;
      	    else
      	      {
      		n = n * -1;
      		if (n > 30)
      		  {
      		    qwx = qwx / (double) (2 << (30-1));
      		    n -= 30;
      		  }
      		qwx = qwx / (double) (2 << (n-1));
      	      }
      	  }
      	else if (n > 0)
      	  {
      	    qwx = qwx * (double) (2 << (n-1));
      	  }
      }
      //
      double qw_r3_N = qgaus_W[i] * rPrime3 * coeff * qwx;
      //exp(-((g - gPrime) * (g - gPrime) / (2 * STDEV * STDEV)));

      double xyz2 = r_point * sinb;
      double xyz0 = r_point * cosb * cosl - LBR_R;
      double xyz1 = r_point * cosb * sinl;

      //double rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_sq_inv);
      double rg = (xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_sq_inv);
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
      	//rsqrt(y) -> http://en.wikipedia.org/wiki/Fast_inverse_square_root
      	double x  = (double) rsqrt(value);
	//double x  = (double) rsqrt((float) rg);
      	/* //fsqrtd */
      	x = x * (3.0 - rg*(x*x));
      	double res = x * rg;
      	rg = res * (0.75 - 0.0625*(res*x));
      //
      }
      //rg = sqrt(rg);
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
	  //cephes dp library for exp
	  double x = -(xyz_norm) * inv_fstream_sigma_sq2[j];
	  {
	    double px = (double) (int)(x * 1.4426950408889634073599246
	  			       + 0.5) - 1; //log2e, estimated floor
	    int n = px;
	    x -= px * 6.93145751953125E-1; //C1
	    x -= px * 1.42860682030941723212E-6; //C2
	    double xx = x * x;
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
	    if (n < 0)
	      {
	  	if (n < -60)
	  	  x = -0.0;
	  	else
	  	  {
	  	    n = n * -1;
	  	    if (n > 30)
	  	      {
	  		x = x / (double) (2 << (30-1));
	  		n -= 30;
	  	      }
	  	    x = x / (double) (2 << (n-1));
	  	  }
	      }
	    else if (n > 0)
	      {
	  	x = x * (double) (2 << (n-1));
	      }
	  }
	  //double result = qw_r3_N *
	  //exp(-(xyz_norm) * inv_fstream_sigma_sq2[j]);
	  double result = qw_r3_N * x;
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
