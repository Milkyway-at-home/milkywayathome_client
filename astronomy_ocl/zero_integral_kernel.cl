#pragma OPENCL EXTENSION cl_khr_fp64: enable

__kernel void
zero_integral_kernel(const int number_streams,
		     __global double *g_bg_int,
		     __global double *g_st_int)
{
  g_bg_int[get_global_id(0)] = 0.0;
  for(int i = 0;i<number_streams;++i)
    g_st_int[(i * get_global_size(0)) + get_global_id(0)] = 0.0;
}
