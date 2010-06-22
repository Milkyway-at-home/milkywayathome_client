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

#include "evaluation_ocl.h"
#include "evaluation_ocl_priv.h"

int determine_work_size(int local_work_size, int desired_size)
{
    if (desired_size % local_work_size != 0)
    {
        return local_work_size * ceil(desired_size / (double) local_work_size);
    }
    else
    {
        return desired_size;
    }
}

void execute_zero_integral_kernel(ocl_mem_t* ocl_mem,
                                  int integral,
                                  int number_streams,
                                  size_t* global_work_size)
{
    cl_int err;
    const char* program_source = read_kernel(ZERO_INTEGRAL_KERNEL);
    cl_program program = clCreateProgramWithSource(ocl_mem->context, 1,
                         &program_source, 0, &err);
    check_error(err);
    build_kernel(program, ocl_mem->devices, "zero_integral.bin");
    //create kernel
    cl_kernel kernel = clCreateKernel(program,
                                      ZERO_INTEGRAL_KERNEL_NAME, &err);
    check_error(err);
    check_error(clSetKernelArg(kernel, 0, sizeof(cl_int),
                               (void*)&(number_streams)));
    check_error(clSetKernelArg(kernel, 1, sizeof(cl_mem),
                               (void*)&ocl_mem->bg_int[integral]));
    check_error(clSetKernelArg(kernel, 2, sizeof(cl_mem),
                               (void*)&ocl_mem->st_int[integral]));
    check_error(clEnqueueNDRangeKernel(ocl_mem->queue, kernel, 1, 0,
                                       global_work_size,
                                       0, 0, 0, 0));
    check_error(clReleaseKernel(kernel));
    check_error(clReleaseProgram(program));
    delete [] program_source;

}

void execute_integral_kernel(ocl_mem_t* ocl_mem,
                             double* parameters,
                             ASTRONOMY_PARAMETERS* ap,
                             double* st_int,
                             double* bg_int)
{
    cl_int err;
    const char* program_source;
    if (ap->number_streams == 2)
        program_source = read_kernel(NVIDIA_INTEGRAL_KERNEL2);
    else
        program_source = read_kernel(NVIDIA_INTEGRAL_KERNEL);
    if (ocl_mem->platform == ATI)
        program_source = read_kernel(ATI_INTEGRAL_KERNEL);
    cl_program program = clCreateProgramWithSource(ocl_mem->context, 1,
                         &program_source, 0, &err);
    check_error(err);
    build_kernel(program, ocl_mem->devices, "integral.bin");
    //create kernel
    cl_kernel kernel = clCreateKernel(program, INTEGRAL_KERNEL_NAME, &err);
    check_error(err);
    double q_sq_inv = 1 / (parameters[0] * parameters[0]);
    double r0 = (double) parameters[1];
    check_error(clSetKernelArg(kernel, 0, sizeof(cl_int),
                               (void*)&(ap->convolve)));
    check_error(clSetKernelArg(kernel, 1, sizeof(cl_int),
                               (void*)&(ap->number_streams)));
    check_error(clSetKernelArg(kernel, 4, sizeof(cl_double),
                               (void*)&q_sq_inv));
    check_error(clSetKernelArg(kernel, 5, sizeof(cl_double), (void*)&r0));
    check_error(clSetKernelArg(kernel, 13, sizeof(cl_mem),
                               (void*)&ocl_mem->fstream_c));
    check_error(clSetKernelArg(kernel, 14, sizeof(cl_mem),
                               (void*)&ocl_mem->fstream_a));
    check_error(clSetKernelArg(kernel, 15, sizeof(cl_mem),
                               (void*)&ocl_mem->inv_fstream_sigma_sq2));
    for (int i = 0; i < ap->number_streams; ++i)
        st_int[i] = 0.0;
    *bg_int = 0.0;
    for (int i = 0; i < ap->number_integrals; ++i)
    {
        int desired_work_size = ap->integral[i]->mu_steps *
                                ap->integral[i]->nu_steps;
        size_t global_work_size[] = {determine_work_size(LOCAL_WORK_SIZE,
                                     desired_work_size)
                                    };
        size_t local_work_size[] = {LOCAL_WORK_SIZE};
        printf("Executing integral %d of size %d\n", i,
               (int)global_work_size[0]);
        execute_zero_integral_kernel(ocl_mem, i,
                                     ap->number_streams, global_work_size);
        check_error(clSetKernelArg(kernel, 3, sizeof(cl_int),
                                   (void*)&(ap->integral[i]->nu_steps)));
        check_error(clSetKernelArg(kernel, 6, sizeof(cl_mem),
                                   (void*)&ocl_mem->r_point[i]));
        check_error(clSetKernelArg(kernel, 7, sizeof(cl_mem),
                                   (void*)&ocl_mem->qw_r3_N[i]));
        check_error(clSetKernelArg(kernel, 8, sizeof(cl_mem),
                                   (void*)&ocl_mem->sinb[i]));
        check_error(clSetKernelArg(kernel, 9, sizeof(cl_mem),
                                   (void*)&ocl_mem->sinl[i]));
        check_error(clSetKernelArg(kernel, 10, sizeof(cl_mem),
                                   (void*)&ocl_mem->cosb[i]));
        check_error(clSetKernelArg(kernel, 11, sizeof(cl_mem),
                                   (void*)&ocl_mem->cosl[i]));
        check_error(clSetKernelArg(kernel, 12, sizeof(cl_mem),
                                   (void*)&ocl_mem->v[i]));
        check_error(clSetKernelArg(kernel, 16, sizeof(cl_mem),
                                   (void*)&ocl_mem->bg_int[i]));
        check_error(clSetKernelArg(kernel, 17, sizeof(cl_mem),
                                   (void*)&ocl_mem->st_int[i]));
        double* st_int_results =
            new double[global_work_size[0] * ap->number_streams];
        double* bg_int_results = new double[global_work_size[0]];
        //make the queue empty before executing the kernel
        check_error(clFinish(ocl_mem->queue));
        cl_event event;
        unsigned long start_time = 0;
        //ap->integral[i]->r_steps = 1;
        for (int r_step = 0; r_step < ap->integral[i]->r_steps; ++r_step)
        {
            if (r_step % 50 == 0)
                printf("Executing r_step %d\n", r_step);
            check_error(clSetKernelArg(kernel, 2, sizeof(cl_int),
                                       (void*)&r_step));
            if (ap->number_streams != 2 || ocl_mem->platform == ATI)
            {
                check_error(clSetKernelArg(kernel, 18,
                                           sizeof(double) * 3 *
                                           ap->number_streams, 0));
                check_error(clSetKernelArg(kernel, 19,
                                           sizeof(double) * 3 *
                                           ap->number_streams, 0));
                check_error(clSetKernelArg(kernel, 20,
                                           sizeof(double) * LOCAL_WORK_SIZE *
                                           ap->number_streams, 0));
            }
            // execute kernel
            check_error(clEnqueueNDRangeKernel(ocl_mem->queue,
                                               kernel,
                                               1,
                                               NULL,
                                               global_work_size,
                                               NULL,
                                               //local_work_size,
                                               0,
                                               NULL,
                                               &event));
            // block until the kernel completes
            check_error(clFinish(ocl_mem->queue));
            if (r_step == 0)
            {
                check_error(clGetEventProfilingInfo(event,
                                                    CL_PROFILING_COMMAND_START,
                                                    sizeof(start_time),
                                                    &start_time, 0));
            }
        }
        unsigned long end_time = 0;
        check_error(clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(end_time),
                                            &end_time, 0));
        unsigned long difference = end_time - start_time;
        printf("Took %f seconds\n", difference / 1e9);
        //copy results from device back to host
        check_error(clEnqueueReadBuffer(ocl_mem->queue, ocl_mem->bg_int[i],
                                        CL_TRUE, 0,
                                        global_work_size[0]
                                        * sizeof(cl_double),
                                        bg_int_results, 0, 0, 0));
        check_error(clEnqueueReadBuffer(ocl_mem->queue, ocl_mem->st_int[i],
                                        CL_TRUE, 0,
                                        ap->number_streams *
                                        global_work_size[0]
                                        * sizeof(cl_double),
                                        st_int_results, 0, 0, 0));
        for (int j = 0; j < desired_work_size; ++j)
        {
            for (int k = 0; k < ap->number_streams; ++k)
            {
                st_int[k] += st_int_results[(k * global_work_size[0]) + j];
            }
            *bg_int += bg_int_results[j];
        }
        delete [] st_int_results;
        delete [] bg_int_results;
    }
    check_error(clReleaseKernel(kernel));
    check_error(clReleaseProgram(program));
}

void execute_zero_likelihood_kernel(ocl_mem_t* ocl_mem,
                                    size_t* global_work_size)
{
    cl_int err;
    const char* program_source = read_kernel(ZERO_LIKELIHOOD_KERNEL);
    cl_program program = clCreateProgramWithSource(ocl_mem->context, 1,
                         &program_source, 0, &err);
    check_error(err);
    build_kernel(program, ocl_mem->devices, "zero_likelihood.bin");
    //create kernel
    cl_kernel kernel = clCreateKernel(program,
                                      ZERO_LIKELIHOOD_KERNEL_NAME, &err);
    check_error(err);
    check_error(clSetKernelArg(kernel, 0, sizeof(cl_mem),
                               (void*)&ocl_mem->probability));
    check_error(clEnqueueNDRangeKernel(ocl_mem->queue, kernel, 1, 0,
                                       global_work_size,
                                       0, 0, 0, 0));
    check_error(clReleaseKernel(kernel));
    check_error(clReleaseProgram(program));
    delete [] program_source;
}

double execute_likelihood_kernel(ocl_mem_t* ocl_mem,
                                 double* parameters,
                                 ASTRONOMY_PARAMETERS* ap,
                                 int number_stars,
                                 size_t* global_work_size)
{
    cl_int err;
    const char* program_source;
    if (ocl_mem->platform == ATI)
        program_source = read_kernel(ATI_LIKELIHOOD_KERNEL);
    else
        program_source = read_kernel(NVIDIA_LIKELIHOOD_KERNEL);
    cl_program program = clCreateProgramWithSource(ocl_mem->context, 1,
                         &program_source, 0, &err);
    check_error(err);
    build_kernel(program, ocl_mem->devices, "likelihood.bin");
    //create kernel
    cl_kernel kernel = clCreateKernel(program, LIKELIHOOD_KERNEL_NAME, &err);
    check_error(err);
    double q_sq_inv = 1 / (parameters[0] * parameters[0]);
    double r0 = (double) parameters[1];
    double coeff = 1.0 / (d_stdev * sqrt(2.0 * D_PI));
    check_error(clSetKernelArg(kernel, 0, sizeof(cl_int),
                               (void*)&(ap->convolve)));
    check_error(clSetKernelArg(kernel, 1, sizeof(cl_int),
                               (void*)&(ap->number_streams)));
    check_error(clSetKernelArg(kernel, 2, sizeof(cl_double),
                               (void*)&q_sq_inv));
    check_error(clSetKernelArg(kernel, 3, sizeof(cl_double),
                               (void*)&r0));
    check_error(clSetKernelArg(kernel, 4, sizeof(cl_double),
                               (void*)&coeff));
    check_error(clSetKernelArg(kernel, 5, sizeof(cl_mem),
                               (void*)&ocl_mem->bg_weight));
    check_error(clSetKernelArg(kernel, 6, sizeof(cl_mem),
                               (void*)&ocl_mem->st_weight));
    check_error(clSetKernelArg(kernel, 7, sizeof(cl_mem),
                               (void*)&ocl_mem->fstream_c));
    check_error(clSetKernelArg(kernel, 8, sizeof(cl_mem),
                               (void*)&ocl_mem->fstream_a));
    check_error(clSetKernelArg(kernel, 9, sizeof(cl_mem),
                               (void*)&ocl_mem->inv_fstream_sigma_sq2));
    check_error(clSetKernelArg(kernel, 10, sizeof(cl_mem),
                               (void*)&ocl_mem->dx));
    check_error(clSetKernelArg(kernel, 11, sizeof(cl_mem),
                               (void*)&ocl_mem->qgaus_W));
    check_error(clSetKernelArg(kernel, 12, sizeof(cl_mem),
                               (void*)&ocl_mem->stars));
    check_error(clSetKernelArg(kernel, 13, sizeof(cl_mem),
                               (void*)&ocl_mem->probability));
    check_error(clSetKernelArg(kernel, 14,
                               sizeof(double) * LOCAL_WORK_SIZE *
                               ap->number_streams, 0));
    check_error(clSetKernelArg(kernel, 15, sizeof(cl_mem),
                               (void*)&ocl_mem->A));
    //make the queue empty before executing the kernel
    check_error(clFinish(ocl_mem->queue));
    cl_event event;
    unsigned long start_time = 0;
    size_t local_work_size[] = {LOCAL_WORK_SIZE};
    // execute kernel
    check_error(clEnqueueNDRangeKernel(ocl_mem->queue, kernel,
                                       1,
                                       0,
                                       global_work_size,
                                       NULL,
                                       //local_work_size,
                                       0,
                                       0,
                                       &event));
    // block until the kernel completes
    check_error(clFinish(ocl_mem->queue));
    check_error(clGetEventProfilingInfo(event,
                                        CL_PROFILING_COMMAND_START,
                                        sizeof(start_time),
                                        &start_time, 0));
    unsigned long end_time = 0;
    check_error(clGetEventProfilingInfo(event,
                                        CL_PROFILING_COMMAND_END,
                                        sizeof(end_time),
                                        &end_time, 0));
    unsigned long difference = end_time - start_time;
    printf("Took %f seconds\n", difference / 1e9);
    //copy results from device back to host
    double* probability = new double[global_work_size[0]];
    check_error(clEnqueueReadBuffer(ocl_mem->queue, ocl_mem->probability,
                                    CL_TRUE, 0,
                                    global_work_size[0]
                                    * sizeof(cl_double),
                                    probability, 0, 0, 0));
    double likelihood = 0.0;
    for (int j = 0; j < number_stars; ++j)
    {
        if (probability[j] > 0)
            likelihood += log10(probability[j]);
        else
            likelihood += probability[j];
        //printf("%i: probability=%.15f\n", j, probability[j]);
    }
    likelihood /= number_stars;
    delete [] probability;
    check_error(clReleaseKernel(kernel));
    check_error(clReleaseProgram(program));
    return likelihood;
}

double ocl_evaluate(double* parameters,
                    ASTRONOMY_PARAMETERS* ap,
                    EVALUATION_STATE* es,
                    STAR_POINTS* sp)
{
    double likelihood;
    ocl_mem_t* ocl_mem;

    ocl_mem = setup_ocl_mem(ap, sp);
    likelihood = ocl_likelihood(parameters, ap, sp, ocl_mem);
    destruct_ocl(ocl_mem);

    return likelihood;
}

double ocl_likelihood(double* parameters,
                      ASTRONOMY_PARAMETERS* ap,
                      STAR_POINTS* sp,
                      ocl_mem_t* ocl_mem)
{
    setup_ocl(parameters, ocl_mem,
              ap->number_streams, ap->sgr_coordinates,
              ap->wedge);
    double* st_int = new double[ap->number_streams];
    double* bg_int = new double;
    execute_integral_kernel(ocl_mem, parameters, ap,
                            st_int, bg_int);
    //st_int[0] = 10.713115282597782;
    //st_int[1] = 513.270626572620472;
    //*bg_int = 0.000471902457355;
    for (int k = 0; k < ap->number_streams; ++k)
    {
        printf("st_int[%d] = %.15f\n", k, st_int[k]);
    }
    printf("bg_int = %.15f\n", *bg_int);
    setup_weights(ocl_mem, parameters, bg_int, st_int,
                  ap->number_streams);
    delete [] st_int;
    delete bg_int;
    int ocl_num_stars = determine_work_size(LOCAL_WORK_SIZE,
                                            sp->number_stars);
    size_t global_work_size[] = {ocl_num_stars};
    execute_zero_likelihood_kernel(ocl_mem,
                                   global_work_size);
    double likelihood = execute_likelihood_kernel(ocl_mem,
                        parameters,
                        ap,
                        sp->number_stars,
                        global_work_size);
    return likelihood;
}

