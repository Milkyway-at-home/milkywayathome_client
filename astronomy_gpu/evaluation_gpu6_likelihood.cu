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
	sxyz0 = xyz0 - tex2D_double(tex_fstream_c, 0, j);
	sxyz1 = xyz1 - tex2D_double(tex_fstream_c, 1, j);
	sxyz2 = xyz2 - tex2D_double(tex_fstream_c, 2, j);

	dotted = tex2D_double(tex_fstream_a, 0,j) * sxyz0 
	  + tex2D_double(tex_fstream_a, 1, j) * sxyz1
	  + tex2D_double(tex_fstream_a, 2, j) * sxyz2;

	sxyz0 -= dotted * tex2D_double(tex_fstream_a, 0, j);
	sxyz1 -= dotted * tex2D_double(tex_fstream_a, 1, j);
	sxyz2 -= dotted * tex2D_double(tex_fstream_a, 2, j);
			
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

    if (probability_sum == 0.0) 
      probability_sum = -238.0;
    else 
      probability_sum = log10(probability_sum);

#ifndef SINGLE_PRECISION
    //probability[threadIdx.x] += probability_sum;
    probability[threadIdx.x] += probability_sum;
#else
    pos = threadIdx.x;
    corrected_next_term = probability_sum - probability_correction[pos];
    new_sum = probability[pos] + corrected_next_term;
    probability_correction[pos] = (new_sum - probability[pos]) - corrected_next_term;
    probability[pos] = new_sum;
#endif
  }
