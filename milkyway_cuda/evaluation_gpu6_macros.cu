#define EXECUTE_ZERO_INTEGRALS(bg_int, st_int)              \
    switch(number_streams) {                      \
    case 1:gpu__zero_integrals<1><<<dimGrid, dimBlock>>>          \
        (offset,                              \
         mu_steps,                            \
         nu_steps,                            \
         total_threads,                           \
         (bg_int),                            \
         (st_int));                           \
        break;                              \
    case 2:gpu__zero_integrals<2><<<dimGrid, dimBlock>>>          \
        (offset,                              \
         mu_steps,                            \
         nu_steps,                            \
         total_threads,                           \
         (bg_int),                            \
         (st_int));                           \
        break;                              \
    case 3:gpu__zero_integrals<3><<<dimGrid, dimBlock>>>          \
        (offset,                              \
         mu_steps,                            \
         nu_steps,                            \
         total_threads,                           \
         (bg_int),                            \
         (st_int));                           \
        break;                              \
    case 4:gpu__zero_integrals<4><<<dimGrid, dimBlock>>>          \
        (offset,                              \
         mu_steps,                            \
         nu_steps,                            \
         total_threads,                           \
         (bg_int),                            \
         (st_int));                           \
        break;                              \
    };                                    \
    {                                 \
        cudaError_t err;                            \
        err = cudaThreadSynchronize();                  \
        if(err != cudaSuccess)                      \
        {                                 \
            fprintf(stderr, "Error executing gpu__zero_integrals error message: %s\n", \
                    cudaGetErrorString(err));               \
        }                                 \
    }

#ifndef SINGLE_PRECISION
#ifdef GLOBAL_MEMORY
#define EXECUTE_INTEGRAL_KERNEL                     \
    switch(number_streams) {                      \
    case 1: gpu__integral_kernel3<1, MAX_CONVOLVE, aux_bg_profile>    \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step, mu_steps,   \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_int_c[i],  \
                device_st_int[i],    \
                device_st_int_c[i],  \
                device_fstream_c,    \
                device_fstream_a);   \
        break;                              \
    case 2: gpu__integral_kernel3<2, MAX_CONVOLVE, aux_bg_profile>    \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_int_c[i],  \
                device_st_int[i],    \
                device_st_int_c[i],  \
                device_fstream_c,    \
                device_fstream_a);   \
        break;                              \
    case 3: gpu__integral_kernel3<3, MAX_CONVOLVE, aux_bg_profile>    \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_int_c[i],  \
                device_st_int[i],    \
                device_st_int_c[i],  \
                device_fstream_c,    \
                device_fstream_a);   \
        break;                              \
    case 4: gpu__integral_kernel3<4, MAX_CONVOLVE, aux_bg_profile>    \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_int_c[i],  \
                device_st_int[i],    \
                device_st_int_c[i],  \
                device_fstream_c,    \
                device_fstream_a);   \
        break;                              \
    }
\
#else
#define EXECUTE_INTEGRAL_KERNEL                     \
    switch(number_streams) {                      \
    case 1: gpu__integral_kernel3<1, MAX_CONVOLVE, aux_bg_profile>    \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step, mu_steps,   \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_st_int[i]);   \
        break;                              \
    case 2: gpu__integral_kernel3<2, MAX_CONVOLVE, aux_bg_profile>    \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_st_int[i]);   \
        break;                              \
    case 3: gpu__integral_kernel3<3, MAX_CONVOLVE, aux_bg_profile>    \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_st_int[i]);   \
        break;                              \
    case 4: gpu__integral_kernel3<4, MAX_CONVOLVE, aux_bg_profile>    \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_st_int[i]);   \
        break;                              \
    }
\
#endif
#else    //single precision
#ifdef GLOBAL_MEMORY
#define EXECUTE_INTEGRAL_KERNEL                     \
    switch(number_streams) {                      \
    case 1:                               \
        gpu__integral_kernel3<1, MAX_CONVOLVE, aux_bg_profile>      \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_correction[i], \
                device_st_int[i],    \
                device_st_correction[i], \
                device_fstream_c,    \
                device_fstream_a);   \
        break;                              \
    case 2:                               \
        gpu__integral_kernel3<2, MAX_CONVOLVE, aux_bg_profile>      \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_correction[i], \
                device_st_int[i],    \
                device_st_correction[i], \
                device_fstream_c,    \
                device_fstream_a);   \
        break;                              \
    case 3:                               \
        gpu__integral_kernel3<3, MAX_CONVOLVE, aux_bg_profile>      \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_correction[i], \
                device_st_int[i],    \
                device_st_correction[i], \
                device_fstream_c,    \
                device_fstream_a);   \
        break;                              \
    case 4:                               \
        gpu__integral_kernel3<4, MAX_CONVOLVE, aux_bg_profile>      \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_correction[i], \
                device_st_int[i],    \
                device_st_correction[i], \
                device_fstream_c,    \
                device_fstream_a);   \
        break;                              \
    }                                 \
     
#else

#define EXECUTE_INTEGRAL_KERNEL                     \
    switch(number_streams) {                      \
    case 1:                               \
        gpu__integral_kernel3<1, MAX_CONVOLVE, aux_bg_profile>      \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_correction[i], \
                device_st_int[i],    \
                device_st_correction[i]); \
        break;                              \
    case 2:                               \
        gpu__integral_kernel3<2, MAX_CONVOLVE, aux_bg_profile>      \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_correction[i], \
                device_st_int[i],    \
                device_st_correction[i]); \
        break;                              \
    case 3:                               \
        gpu__integral_kernel3<3, MAX_CONVOLVE, aux_bg_profile>      \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_correction[i], \
                device_st_int[i],    \
                device_st_correction[i]); \
        break;                              \
    case 4:                               \
        gpu__integral_kernel3<4, MAX_CONVOLVE, aux_bg_profile>      \
        <<<dimGrid, dimBlock, shared_mem_size>>>(mu_step,mu_steps,    \
                j, r_steps,      \
                nu_steps, total_threads, \
                bg_a, bg_b, bg_c,    \
                q_squared_inverse, r0,   \
                device_sinb[i],      \
                device_sinl[i],      \
                device_cosb[i],      \
                device_cosl[i],      \
                device_V[i],     \
                device_bg_int[i],    \
                device_bg_correction[i], \
                device_st_int[i],    \
                device_st_correction[i]); \
        break;                              \
    }                                 \
     
#endif
#endif

#ifndef SINGLE_PRECISION
#define EXECUTE_LIKELIHOOD_KERNEL                   \
    \
    switch (ap->number_streams) {                     \
    case 1: gpu__likelihood_kernel<1, aux_bg_profile><<<dimGrid, dimBlock>>> \
        (offset, convolve,                        \
         bg_a, bg_b, bg_c,                        \
         q_squared_inverse, r0,                       \
         (GPU_PRECISION)coeff,                        \
         device_stars,                            \
         number_stars,                            \
         device_bg_only, device_st_only,                  \
         device_probability);                     \
        break;                              \
    case 2: gpu__likelihood_kernel<2, aux_bg_profile><<<dimGrid, dimBlock>>> \
        (offset, convolve,                        \
         bg_a, bg_b, bg_c,                        \
         q_squared_inverse, r0,                       \
         (GPU_PRECISION)coeff,                        \
         device_stars,                            \
         number_stars,                            \
         device_bg_only, device_st_only,                  \
         device_probability);                     \
        break;                              \
    case 3: gpu__likelihood_kernel<3, aux_bg_profile><<<dimGrid, dimBlock>>> \
        (offset, convolve,                        \
         bg_a, bg_b, bg_c,                        \
         q_squared_inverse, r0,                       \
         (GPU_PRECISION)coeff,                        \
         device_stars,                            \
         number_stars,                            \
         device_bg_only, device_st_only,                  \
         device_probability);                     \
        break;                              \
    case 4: gpu__likelihood_kernel<4, aux_bg_profile><<<dimGrid, dimBlock>>> \
        (offset, convolve,                        \
         bg_a, bg_b, bg_c,                        \
         q_squared_inverse, r0,                       \
         (GPU_PRECISION)coeff,                        \
         device_stars,                            \
         number_stars,                            \
         device_bg_only, device_st_only,                  \
         device_probability);                     \
        break;                              \
    }                                 \
     
#else
#define EXECUTE_LIKELIHOOD_KERNEL                   \
    switch (ap->number_streams) {                     \
    case 1: gpu__likelihood_kernel<1, aux_bg_profile><<<dimGrid, dimBlock>>> \
        (offset, convolve,                        \
         bg_a, bg_b, bg_c,                        \
         q_squared_inverse, r0,                       \
         (GPU_PRECISION)coeff,                        \
         device_stars,                            \
         number_stars,                            \
         device_bg_only, device_st_only,                  \
         device_probability);                     \
        break;                              \
    case 2: gpu__likelihood_kernel<2, aux_bg_profile><<<dimGrid, dimBlock>>> \
        (offset, convolve,                        \
         bg_a, bg_b, bg_c,                        \
         q_squared_inverse, r0,                       \
         (GPU_PRECISION)coeff,                        \
         device_stars,                            \
         number_stars,                            \
         device_bg_only, device_st_only,                  \
         device_probability);                     \
        break;                              \
    case 3: gpu__likelihood_kernel<3, aux_bg_profile><<<dimGrid, dimBlock>>> \
        (offset, convolve,                        \
         bg_a, bg_b, bg_c,                        \
         q_squared_inverse, r0,                       \
         (GPU_PRECISION)coeff,                        \
         device_stars,                            \
         number_stars,                            \
         device_bg_only, device_st_only,                  \
         device_probability);                     \
        break;                              \
    case 4: gpu__likelihood_kernel<4, aux_bg_profile><<<dimGrid, dimBlock>>> \
        (offset, convolve,                        \
         bg_a, bg_b, bg_c,                        \
         q_squared_inverse, r0,                       \
         (GPU_PRECISION)coeff,                        \
         device_stars,                            \
         number_stars,                            \
         device_bg_only, device_st_only,                  \
         device_probability);                     \
        break;                              \
    }                                 \
     
#endif


