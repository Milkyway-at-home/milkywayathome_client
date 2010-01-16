#define EXECUTE_ZERO_INTEGRALS(bg_int, st_int)				\
  switch(number_streams) {						\
  case 1:gpu__zero_integrals<1><<<dimGrid, dimBlock>>>			\
      (mu_step,								\
       mu_steps[i],							\
       nu_steps[i],							\
       (bg_int),							\
       (st_int));							\
    break;								\
  case 2:gpu__zero_integrals<2><<<dimGrid, dimBlock>>>			\
      (mu_step,								\
       mu_steps[i],							\
       nu_steps[i],							\
       (bg_int),							\
       (st_int));							\
    break;								\
  case 3:gpu__zero_integrals<3><<<dimGrid, dimBlock>>>			\
      (mu_step,								\
       mu_steps[i],							\
       nu_steps[i],							\
       (bg_int),							\
       (st_int));							\
    break;								\
  case 4:gpu__zero_integrals<4><<<dimGrid, dimBlock>>>			\
      (mu_step,								\
       mu_steps[i],							\
       nu_steps[i],							\
       (bg_int),							\
       (st_int));							\
    break;								\
  };									\
  {									\
    cudaError_t err;							\
    err = cudaThreadSynchronize();					\
    if(err != cudaSuccess)						\
      {									\
	fprintf(stderr, "Error executing gpu__zero_integrals error message: %s\n", \
		cudaGetErrorString(err));				\
      }									\
  }

#ifndef SINGLE_PRECISION
#define EXECUTE_INTEGRAL_KERNEL						\
  switch(number_streams) {						\
  case 1: gpu__integral_kernel3<1, MAX_CONVOLVE>			\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step, mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__stream_integrals[i]); \
    break;								\
  case 2: gpu__integral_kernel3<2, MAX_CONVOLVE>			\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__stream_integrals[i]); \
    break;								\
  case 3: gpu__integral_kernel3<3, MAX_CONVOLVE>			\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__stream_integrals[i]); \
    break;								\
  case 4: gpu__integral_kernel3<4, MAX_CONVOLVE>			\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__stream_integrals[i]); \
    break;								\
  }									\

#else	

								
#define EXECUTE_INTEGRAL_KERNEL						\
  switch(number_streams) {						\
  case 1:								\
    gpu__integral_kernel3<1, MAX_CONVOLVE>				\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__background_correction[i], \
					       device__stream_integrals[i], \
					       device__stream_correction[i]); \
    break;								\
  case 2:								\
    gpu__integral_kernel3<2, MAX_CONVOLVE>				\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__background_correction[i], \
					       device__stream_integrals[i], \
					       device__stream_correction[i]); \
    break;								\
  case 3:								\
    gpu__integral_kernel3<3, MAX_CONVOLVE>				\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__background_correction[i], \
					       device__stream_integrals[i], \
					       device__stream_correction[i]); \
    break;								\
  case 4:								\
    gpu__integral_kernel3<4, MAX_CONVOLVE>				\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__background_correction[i], \
					       device__stream_integrals[i], \
					       device__stream_correction[i]); \
    break;								\
  }									\

#endif								

#ifndef SINGLE_PRECISION						
#define EXECUTE_AUX_INTEGRAL_KERNEL					\
  switch(number_streams) {						\
  case 1: gpu__integral_kernel3_aux<1, MAX_CONVOLVE>			\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step, mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       BG_A, BG_B, BG_C,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__stream_integrals[i]); \
    break;								\
  case 2: gpu__integral_kernel3_aux<2, MAX_CONVOLVE>			\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       BG_A, BG_B, BG_C,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i], \
					       device__stream_integrals[i]); \
    break;								\
  case 3: gpu__integral_kernel3_aux<3, MAX_CONVOLVE>			\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       BG_A, BG_B, BG_C,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i],	\
					       device__stream_integrals[i]); \
    break;								\
  case 4: gpu__integral_kernel3_aux<4, MAX_CONVOLVE>			\
      <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					       j, r_steps[i],		\
					       nu_steps[i],		\
					       q_squared_inverse, r0,	\
					       BG_A, BG_B, BG_C,	\
					       device__sinb[i],		\
					       device__sinl[i],		\
					       device__cosb[i],		\
					       device__cosl[i],		\
					       device__V[i],		\
					       device__background_integrals[i],	\
					       device__stream_integrals[i]); \
    break;								\
  }									\

#else

#define EXECUTE_AUX_INTEGRAL_KERNEL					\
switch(number_streams) {						\
 case 1:								\
   gpu__integral_kernel3_aux<1, MAX_CONVOLVE>				\
     <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					      j, r_steps[i],		\
					      nu_steps[i],		\
					      q_squared_inverse, r0,	\
					      BG_A, BG_B, BG_C,		\
					      device__sinb[i],		\
					      device__sinl[i],		\
					      device__cosb[i],		\
					      device__cosl[i],		\
					      device__V[i],		\
					      device__background_integrals[i], \
					      device__background_correction[i], \
					      device__stream_integrals[i], \
					      device__stream_correction[i]); \
   break;								\
 case 2:								\
   gpu__integral_kernel3_aux<2, MAX_CONVOLVE>				\
     <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					      j, r_steps[i],		\
					      nu_steps[i],		\
					      q_squared_inverse, r0,	\
					      BG_A, BG_B, BG_C,		\
					      device__sinb[i],		\
					      device__sinl[i],		\
					      device__cosb[i],		\
					      device__cosl[i],		\
					      device__V[i],		\
					      device__background_integrals[i], \
					      device__background_correction[i], \
					      device__stream_integrals[i], \
					      device__stream_correction[i]); \
   break;								\
 case 3:								\
   gpu__integral_kernel3_aux<3, MAX_CONVOLVE>				\
     <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					      j, r_steps[i],		\
					      nu_steps[i],		\
					      q_squared_inverse, r0,	\
					      BG_A, BG_B, BG_C,		\
					      device__sinb[i],		\
					      device__sinl[i],		\
					      device__cosb[i],		\
					      device__cosl[i],		\
					      device__V[i],		\
					      device__background_integrals[i], \
					      device__background_correction[i], \
					      device__stream_integrals[i], \
					      device__stream_correction[i]); \
   break;								\
 case 4:								\
   gpu__integral_kernel3_aux<4, MAX_CONVOLVE>				\
     <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps[i],	\
					      j, r_steps[i],		\
					      nu_steps[i],		\
					      q_squared_inverse, r0,	\
					      BG_A, BG_B, BG_C,		\
					      device__sinb[i],		\
					      device__sinl[i],		\
					      device__cosb[i],		\
					      device__cosl[i],		\
					      device__V[i],		\
					      device__background_integrals[i], \
					      device__background_correction[i], \
					      device__stream_integrals[i], \
					      device__stream_correction[i]); \
   break;								\
 }									\

#endif
