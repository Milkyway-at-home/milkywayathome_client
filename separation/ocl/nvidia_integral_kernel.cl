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

__kernel void
integral_kernel(const int convolve,
                const int number_streams,
                const int r_step,
                const int nu_steps,
                const double q_sq_inv,
                const double r0,
                __global const double* r_point,
                __global const double* qw_r3_N,
                __global const double* g_sinb,
                __global const double* g_sinl,
                __global const double* g_cosb,
                __global const double* g_cosl,
                __global const double* g_v,
                __constant double* fstream_c,
                __constant double* fstream_a,
                __constant double* inv_fstream_sigma_sq2,
                __global double* g_bg_int,
                __global double* g_st_int,
                __local double* s_fstream_c,
                __local double* s_fstream_a,
                __local double* st_int)
{
    double bg_int = 0.0;

    for (int i = 0; i < number_streams; ++i)
    {
        st_int[(i * get_local_size(0)) + get_local_id(0)] = 0.0;
        if (get_local_id(0) == 0)
        {
            s_fstream_c[(i * 3) + 0] = fstream_c[(i * 3) + 0];
            s_fstream_c[(i * 3) + 1] = fstream_c[(i * 3) + 1];
            s_fstream_c[(i * 3) + 2] = fstream_c[(i * 3) + 2];
            s_fstream_a[(i * 3) + 0] = fstream_a[(i * 3) + 0];
            s_fstream_a[(i * 3) + 1] = fstream_a[(i * 3) + 1];
            s_fstream_a[(i * 3) + 2] = fstream_a[(i * 3) + 2];
        }
    }

    double sinb = g_sinb[get_global_id(0)];
    double sinl = g_sinl[get_global_id(0)];
    double cosb = g_cosb[get_global_id(0)];
    double cosl = g_cosl[get_global_id(0)];


    for (int i = 0; i < convolve; ++i)
    {
        double xyz2 = r_point[(r_step * convolve) + i] * sinb;
        double xyz0 = r_point[(r_step * convolve) + i] * cosb *
                      cosl - LBR_R;
        double xyz1 = r_point[(r_step * convolve) + i] * cosb * sinl;

        //double rg = sqrt(xyz0*xyz0 + xyz1*xyz1 + (xyz2*xyz2) * q_sq_inv);
        //fsqrtd
        double rg = xyz0 * xyz0 + xyz1 * xyz1 + (xyz2 * xyz2) * q_sq_inv;
        //
        {
            double x  = (double) rsqrt((float) rg);
            //fsqrtd
            x = x * (3.0 - rg * (x * x));
            double res = x * rg;
            rg = res * (0.75 - 0.0625 * (res * x));
        }
        //
        double rs = rg + r0;

        //bg_int += qw_r3_N[(r_step * convolve) + i] / (rg * rs * rs *rs);
        //divd
        {
            double a = qw_r3_N[(r_step * convolve) + i];
            double b = (rg * rs * rs * rs);
            double y = (double)(1.0f / (float) b);
            double c = 1.0 - b * y;
            y += y * c;
            double r = a * y;
            c = a - b * r;
            bg_int += r + y * c;
        }

        for (int j = 0; j < number_streams; ++j)
        {
            double sxyz0 = xyz0 - s_fstream_c[(j * 3) + 0];
            double sxyz1 = xyz1 - s_fstream_c[(j * 3) + 1];
            double sxyz2 = xyz2 - s_fstream_c[(j * 3) + 2];

            double dotted = s_fstream_a[(j * 3) + 0] * sxyz0
                            + s_fstream_a[(j * 3) + 1] * sxyz1
                            + s_fstream_a[(j * 3) + 2] * sxyz2;

            sxyz0 -= dotted * s_fstream_a[(j * 3) + 0];
            sxyz1 -= dotted * s_fstream_a[(j * 3) + 1];
            sxyz2 -= dotted * s_fstream_a[(j * 3) + 2];

            double xyz_norm = (sxyz0 * sxyz0) + (sxyz1 * sxyz1)
                              + (sxyz2 * sxyz2);
            double result = qw_r3_N[(r_step * convolve) + i] *
                            exp(-(xyz_norm) * inv_fstream_sigma_sq2[j]);
            st_int[(j * get_local_size(0)) + get_local_id(0)] += result;
        }
    }
    int nu_step = get_global_id(0) % nu_steps;
    double v = g_v[nu_step + (r_step * nu_steps)];
    g_bg_int[get_global_id(0)] += (bg_int * v);
    for (int i = 0; i < number_streams; ++i)
    {
        g_st_int[get_global_id(0) + (i * get_global_size(0))]
        += (st_int[(i * get_local_size(0)) + get_local_id(0)] * v);
    }
}

