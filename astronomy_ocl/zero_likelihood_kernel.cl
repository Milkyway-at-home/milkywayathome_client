#pragma OPENCL EXTENSION cl_khr_fp64: enable

__kernel void
zero_likelihood_kernel(__global double *g_probability)
{
  g_probability[get_global_id(0)] = 0.0;
}
