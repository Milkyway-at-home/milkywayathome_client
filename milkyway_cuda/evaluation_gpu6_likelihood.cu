template <unsigned int number_streams, unsigned int aux_bg_profile>
__global__ void gpu__likelihood_kernel(int offset, int convolve,
				       GPU_PRECISION bg_a, GPU_PRECISION bg_b,
				       GPU_PRECISION bg_c,
				       GPU_PRECISION q_squared_inverse,
				       GPU_PRECISION r0,
				       GPU_PRECISION coeff,
				       GPU_PRECISION *device__stars,
				       int number_stars,
				       GPU_PRECISION *bg_only,
				       GPU_PRECISION *st_only,
				       GPU_PRECISION *probability)
{
  int i;
  int pos = (offset + threadIdx.x + (blockDim.x * blockIdx.x));
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

  GPU_PRECISION reff_xr_rp3 = reff_value * d_xr / rPrime3;

  GPU_PRECISION r_point, qw_r3_N;
  GPU_PRECISION zp, rs, g;

  GPU_PRECISION bg_int = 0.0;
  GPU_PRECISION st_int[number_streams];
  for (i = 0; i < number_streams; i++)
    st_int[i] = 0.0;

  for (i = 0; i < convolve; i++) {
    g = gPrime + constant_dx[i];
    r_point = pow(10.0, (g - d_absm) * .2 + 1.0) * .001;
    rPrime3 = r_point * r_point * r_point;

    qw_r3_N = constant_qgaus_W[i] * rPrime3 * coeff *
      exp(-((g - gPrime) * (g - gPrime) / (2 * d_stdev * d_stdev)));

    xyz2 = r_point * sinb;
    zp = r_point * cosb;
    xyz0 = zp * cosl - d_lbr_r;
    xyz1 = zp * sinl;

    rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_squared_inverse);
    rs = rg + r0;

    if (aux_bg_profile == 1)
      {
	double r_in_mag = g;
	double r_in_mag2 = g*g;
	double h_prob = (qw_r3_N / (rg * rs * rs * rs));
	double aux_prob = qw_r3_N * ( bg_a * r_in_mag2 + bg_b * r_in_mag + bg_c );
	bg_int += h_prob + aux_prob;
	bg_int = pos;
      }
    else
      {
	bg_int += (qw_r3_N / (rg * rs * rs * rs));
      }
    for (int j = 0; j < number_streams; j++) {
      pos = (j * 3);
      sxyz0 = xyz0 - constant_fstream_c[0 + pos];
      sxyz1 = xyz1 - constant_fstream_c[1 + pos];
      sxyz2 = xyz2 - constant_fstream_c[2 + pos];

      dotted = constant_fstream_a[0 + pos] * sxyz0
	+ constant_fstream_a[1 + pos] * sxyz1
	+ constant_fstream_a[2 + pos] * sxyz2;

      sxyz0 -= dotted * constant_fstream_a[0 + pos];
      sxyz1 -= dotted * constant_fstream_a[1 + pos];
      sxyz2 -= dotted * constant_fstream_a[2 + pos];

      st_int[j] += qw_r3_N *
	exp(-((sxyz0 * sxyz0) + (sxyz1 * sxyz1) + (sxyz2 * sxyz2)) *
	    constant_inverse_fstream_sigma_sq2[j]);
    }
  }
  pos = (offset + threadIdx.x + (blockDim.x * blockIdx.x));
  GPU_PRECISION probability_sum = bg_int * constant_background_weight[0];
  if (probability_sum * reff_xr_rp3 == 0.0)
    bg_only[pos] = -238.0;
  else
    bg_only[pos] = log10(probability_sum * reff_xr_rp3);
  //pragma unroll 1 makes the loop not unroll,
  //when it unrolls it causes a launch failure when trying
  //to access constant_stream_weight[i], when i is 1
  #pragma unroll 1
  for (i = 0; i < number_streams; i++) {
    st_only[pos] = st_int[i] * constant_stream_weight[i];
    probability_sum += st_only[pos];
    if (st_only[pos] * reff_xr_rp3 == 0.0)
      st_only[pos] = -238.0;
    else
      st_only[pos] = log10(st_only[pos] * reff_xr_rp3);
    pos += number_stars;
  }
  probability_sum *= reff_xr_rp3;

  if (probability_sum == 0.0)
    probability_sum = -238.0;
  else
    probability_sum = log10(probability_sum);
  pos = (offset + threadIdx.x + (blockDim.x * blockIdx.x));
  probability[pos] = probability_sum;
}
