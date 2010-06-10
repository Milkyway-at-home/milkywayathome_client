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

ocl_mem_t* setup_ocl_mem(ASTRONOMY_PARAMETERS* ap,
                         STAR_POINTS* sp)
{
    ocl_mem_t* ocl_mem = new ocl_mem_t;
    cl_platform_id* platforms = get_platforms();
    char* name = get_platform_info(platforms[0],
                                   CL_PLATFORM_NAME);
    if (strcmp(name, "ATI Stream") == 0)
        ocl_mem->platform = ATI;
    else
        ocl_mem->platform = NVIDIA;
    cl_device_id* devices = get_devices(platforms[0]);
    ocl_mem->devices = devices;
    ocl_mem->context = get_context(devices[0]);
    setup_command_queue(ocl_mem);
    ocl_mem->sinb = new cl_mem[ap->number_integrals];
    ocl_mem->sinl = new cl_mem[ap->number_integrals];
    ocl_mem->cosb = new cl_mem[ap->number_integrals];
    ocl_mem->cosl = new cl_mem[ap->number_integrals];
    ocl_mem->v = new cl_mem[ap->number_integrals];
    ocl_mem->r_point = new cl_mem[ap->number_integrals];
    ocl_mem->qw_r3_N = new cl_mem[ap->number_integrals];
    ocl_mem->bg_int = new cl_mem[ap->number_integrals];
    ocl_mem->st_int = new cl_mem[ap->number_integrals];
    ocl_mem->num_integrals = ap->number_integrals;
    cl_int err;
    for (int i = 0; i < ap->number_integrals; ++i)
    {
        //Setup sinb, sinl, cosb, cosl
        int int_size =
            ap->integral[i]->nu_steps * ap->integral[i]->mu_steps;
        //pad it to 64 work items per work group
        int ocl_int_size = determine_work_size(LOCAL_WORK_SIZE,
                                               int_size);
        double* sinb = new double[int_size];
        double* sinl = new double[int_size];
        double* cosb = new double[int_size];
        double* cosl = new double[int_size];
        populate_lb(ap->sgr_coordinates, ap->wedge,
                    ap->integral[i]->mu_steps,
                    ap->integral[i]->mu_min,
                    ap->integral[i]->mu_step_size,
                    ap->integral[i]->nu_steps,
                    ap->integral[i]->nu_min,
                    ap->integral[i]->nu_step_size,
                    sinb, sinl, cosb, cosl);
        double* fsinb = new double[ocl_int_size];
        double* fsinl = new double[ocl_int_size];
        double* fcosb = new double[ocl_int_size];
        double* fcosl = new double[ocl_int_size];
        for (int j = 0; j < int_size; ++j)
        {
            fsinb[j] = (double) sinb[j];
            fsinl[j] = (double) sinl[j];
            fcosb[j] = (double) cosb[j];
            fcosl[j] = (double) cosl[j];
        }
        for (int j = int_size; j < ocl_int_size; ++j)
        {
            fsinb[j] = 0.0;
            fsinl[j] = 0.0;
            fcosb[j] = 0.0;
            fcosl[j] = 0.0;
        }
        ocl_mem->sinb[i] = clCreateBuffer(ocl_mem->context,
                                          CL_MEM_READ_ONLY |
                                          CL_MEM_COPY_HOST_PTR,
                                          ocl_int_size * sizeof(cl_double),
                                          fsinb, &err);
        check_error(err);
        ocl_mem->sinl[i] = clCreateBuffer(ocl_mem->context,
                                          CL_MEM_READ_ONLY |
                                          CL_MEM_COPY_HOST_PTR,
                                          ocl_int_size * sizeof(cl_double),
                                          fsinl, &err);
        check_error(err);
        ocl_mem->cosb[i] = clCreateBuffer(ocl_mem->context,
                                          CL_MEM_READ_ONLY |
                                          CL_MEM_COPY_HOST_PTR,
                                          ocl_int_size * sizeof(cl_double),
                                          fcosb, &err);
        check_error(err);
        ocl_mem->cosl[i] = clCreateBuffer(ocl_mem->context,
                                          CL_MEM_READ_ONLY |
                                          CL_MEM_COPY_HOST_PTR,
                                          ocl_int_size * sizeof(cl_double),
                                          fcosl, &err);
        check_error(err);
        delete [] sinb;
        delete [] sinl;
        delete [] cosb;
        delete [] cosl;
        delete [] fsinb;
        delete [] fsinl;
        delete [] fcosb;
        delete [] fcosl;
        //Setup r constants
        double* irv = new double[ap->integral[i]->r_steps];
        double* reff_xr_rp3 = new double[ap->integral[i]->r_steps];
        double** qw_r3_N = new double*[ap->integral[i]->r_steps];
        double** r_point = new double*[ap->integral[i]->r_steps];
        double* ids = new double[ap->integral[i]->nu_steps];
        double* nus = new double[ap->integral[i]->nu_steps];
        cpu__r_constants(ap->convolve,
                         ap->integral[i]->r_steps,
                         ap->integral[i]->r_min,
                         ap->integral[i]->r_step_size,
                         ap->integral[i]->mu_steps,
                         ap->integral[i]->mu_min,
                         ap->integral[i]->mu_step_size,
                         ap->integral[i]->nu_steps,
                         ap->integral[i]->nu_min,
                         ap->integral[i]->nu_step_size,
                         irv,
                         r_point,

                         NULL, //FIXME: this is to just get it to compile
                         NULL,

                         qw_r3_N,
                         reff_xr_rp3,
                         nus,
                         ids);
        double* fr_point = new double[ap->integral[i]->r_steps * ap->convolve];
        double* fqw_r3_N = new double[ap->integral[i]->r_steps * ap->convolve];
        for (int j = 0; j < ap->integral[i]->r_steps; ++j)
        {
            for (int k = 0; k < ap->convolve; ++k)
            {
                fr_point[(j * ap->convolve) + k] =
                    r_point[j][k];
                fqw_r3_N[(j * ap->convolve) + k] =
                    qw_r3_N[j][k];
            }
        }
        ocl_mem->r_point[i] = clCreateBuffer(ocl_mem->context,
                                             CL_MEM_READ_ONLY |
                                             CL_MEM_COPY_HOST_PTR,
                                             ap->integral[i]->r_steps *
                                             ap->convolve *
                                             sizeof(cl_double),
                                             fr_point, &err);
        check_error(err);
        ocl_mem->qw_r3_N[i] = clCreateBuffer(ocl_mem->context,
                                             CL_MEM_READ_ONLY |
                                             CL_MEM_COPY_HOST_PTR,
                                             ap->integral[i]->r_steps *
                                             ap->convolve *
                                             sizeof(cl_double),
                                             fqw_r3_N, &err);
        check_error(err);
        delete [] fr_point;
        delete [] fqw_r3_N;
        double* v = new double[ap->integral[i]->nu_steps *
                               ap->integral[i]->r_steps];
        for (int j = 0; j < ap->integral[i]->nu_steps; ++j)
        {
            for (int k = 0; k < ap->integral[i]->r_steps; ++k)
            {
                v[(k * ap->integral[i]->nu_steps) + j] = (double)
                        reff_xr_rp3[k] * irv[k] * ids[j];
            }
        }
        ocl_mem->v[i] = clCreateBuffer(ocl_mem->context,
                                       CL_MEM_READ_ONLY |
                                       CL_MEM_COPY_HOST_PTR,
                                       ap->integral[i]->nu_steps *
                                       ap->integral[i]->r_steps *
                                       sizeof(cl_double),
                                       v, &err);
        check_error(err);
        delete [] v;
        delete [] irv;
        delete [] reff_xr_rp3;
        delete [] qw_r3_N;
        delete [] r_point;
        delete [] ids;
        delete [] nus;
        //Setup device only memory
        ocl_mem->bg_int[i] = clCreateBuffer(ocl_mem->context,
                                            CL_MEM_READ_WRITE,
                                            ocl_int_size *
                                            sizeof(cl_double),
                                            0, &err);
        ocl_mem->st_int[i] = clCreateBuffer(ocl_mem->context,
                                            CL_MEM_READ_WRITE,
                                            ap->number_streams *
                                            ocl_int_size *
                                            sizeof(cl_double),
                                            0, &err);
        check_error(err);
    }
    //setup memory for the likelihood kernel
    int ocl_num_stars = determine_work_size(LOCAL_WORK_SIZE,
                                            sp->number_stars);
    double* stars = new double[ocl_num_stars * 6];
    for (int i = 0; i < sp->number_stars; ++i)
    {
        stars[i + ocl_num_stars*0] = sin(sp->stars[i][1] * D_DEG2RAD);
        stars[i + ocl_num_stars*1] = sin(sp->stars[i][0] * D_DEG2RAD);
        stars[i + ocl_num_stars*2] = cos(sp->stars[i][1] * D_DEG2RAD);
        stars[i + ocl_num_stars*3] = cos(sp->stars[i][0] * D_DEG2RAD);
        stars[i + ocl_num_stars*4] = sp->stars[i][2];
        stars[i + ocl_num_stars*5] = log10(sp->stars[i][2]);
    }
    for (int i = sp->number_stars; i < ocl_num_stars; ++i)
    {
        stars[i + ocl_num_stars*0] = 0.0;
        stars[i + ocl_num_stars*1] = 0.0;
        stars[i + ocl_num_stars*2] = 0.0;
        stars[i + ocl_num_stars*3] = 0.0;
        stars[i + ocl_num_stars*4] = 0.0;
        stars[i + ocl_num_stars*5] = 0.0;
    }
    ocl_mem->stars = clCreateBuffer(ocl_mem->context,
                                    CL_MEM_READ_ONLY |
                                    CL_MEM_COPY_HOST_PTR,
                                    ocl_num_stars * 6 *
                                    sizeof(cl_double),
                                    stars, &err);
    check_error(err);
    delete [] stars;
    double* qgaus_W = new double[ap->convolve];
    double* qgaus_X = new double[ap->convolve];

    d_gauss_legendre(-1.0, 1.0, qgaus_X, qgaus_W, ap->convolve);
    double* dx = new double[ap->convolve];
    double* fqgaus_W = new double[ap->convolve];
    for (int i = 0; i < ap->convolve; ++i)
    {
        dx[i] = (3.0 * d_stdev * qgaus_X[i]);
        fqgaus_W[i] = (double) qgaus_W[i];
    }
    ocl_mem->qgaus_W = clCreateBuffer(ocl_mem->context,
                                      CL_MEM_READ_ONLY |
                                      CL_MEM_COPY_HOST_PTR,
                                      ap->convolve *
                                      sizeof(cl_double),
                                      fqgaus_W, &err);
    check_error(err);
    ocl_mem->dx = clCreateBuffer(ocl_mem->context,
                                 CL_MEM_READ_ONLY |
                                 CL_MEM_COPY_HOST_PTR,
                                 ap->convolve *
                                 sizeof(cl_double),
                                 dx, &err);
    check_error(err);
    delete [] dx;
    delete [] fqgaus_W;
    delete [] qgaus_W;
    delete [] qgaus_X;
    ocl_mem->probability = clCreateBuffer(ocl_mem->context,
                                          CL_MEM_READ_WRITE,
                                          ocl_num_stars *
                                          sizeof(cl_double),
                                          0, &err);
    check_error(err);

    //cephes DP libary
    double A[] =
    {
        1.00000000000000000000E0,
        9.57603280698573700036E-1,
        9.17004043204671215328E-1,
        8.78126080186649726755E-1,
        8.40896415253714502036E-1,
        8.05245165974627141736E-1,
        7.71105412703970372057E-1,
        7.38413072969749673113E-1,
        7.07106781186547572737E-1,
        6.77127773468446325644E-1,
        6.48419777325504820276E-1,
        6.20928906036742001007E-1,
        5.94603557501360513449E-1,
        5.69394317378345782288E-1,
        5.45253866332628844837E-1,
        5.22136891213706877402E-1,
        5.00000000000000000000E-1
    };
    ocl_mem->A = clCreateBuffer(ocl_mem->context,
                                CL_MEM_READ_ONLY |
                                CL_MEM_COPY_HOST_PTR,
                                17 *
                                sizeof(cl_double),
                                A, &err);
    check_error(err);
    return ocl_mem;
}



void setup_ocl(double* parameters, ocl_mem_t* ocl_mem,
               int number_streams, int sgr_coordinates,
               int wedge)
{
    double stream_c[3], lbr[3];
    double* fstream_a = new double[number_streams * 3];
    double* fstream_c = new double[number_streams * 3];
    double* inv_fstream_sigma_sq2 = new double[number_streams];
    for (int i = 0; i < number_streams; ++i)
    {
        inv_fstream_sigma_sq2[i] = 1 / (2.0 * STREAM_PARAMETERS(i, 4)
                                        * STREAM_PARAMETERS(i, 4));
        fstream_a[(i * 3) + 0] = sin(STREAM_PARAMETERS(i, 2)) *
                                 cos(STREAM_PARAMETERS(i, 3));
        fstream_a[(i * 3) + 1] = sin(STREAM_PARAMETERS(i, 2)) *
                                 sin(STREAM_PARAMETERS(i, 3));
        fstream_a[(i * 3) + 2] = cos(STREAM_PARAMETERS(i, 2));
        if (sgr_coordinates == 0)
        {
            gc_eq_gal(wedge, STREAM_PARAMETERS(i, 0) * D_DEG2RAD,
                      0, &(lbr[0]), &(lbr[1]));
        }
        else
        {
            gc_sgr_gal(wedge, STREAM_PARAMETERS(i, 0) * D_DEG2RAD,
                       0, &(lbr[0]), &(lbr[1]));
        }
        lbr[2] = STREAM_PARAMETERS(i, 1);
        d_lbr2xyz(lbr, stream_c);

        fstream_c[(i * 3) + 0] = stream_c[0];
        fstream_c[(i * 3) + 1] = stream_c[1];
        fstream_c[(i * 3) + 2] = stream_c[2];

        // printf("fstream_sigma_sq2[%d] = %.15f\n", i, inv_fstream_sigma_sq2[i]);
        // printf("stream %d fstream_c[0] = %.15f\n", i, fstream_c[i*3]);
        // printf("stream %d fstream_c[1] = %.15f\n", i, fstream_c[i*3 + 1]);
        // printf("stream %d fstream_c[2] = %.15f\n", i, fstream_c[i*3 + 2]);

        // printf("stream %d fstream_a[0] = %.15f\n", i, fstream_a[i*3]);
        // printf("stream %d fstream_a[1] = %.15f\n", i, fstream_a[i*3 + 1]);
        // printf("stream %d fstream_a[2] = %.15f\n", i, fstream_a[i*3 + 2]);
    }
    cl_int err;
    ocl_mem->fstream_a = clCreateBuffer(ocl_mem->context,
                                        CL_MEM_READ_ONLY |
                                        CL_MEM_COPY_HOST_PTR,
                                        number_streams * 3 *
                                        sizeof(cl_double),
                                        fstream_a, &err);
    check_error(err);
    ocl_mem->fstream_c = clCreateBuffer(ocl_mem->context,
                                        CL_MEM_READ_ONLY |
                                        CL_MEM_COPY_HOST_PTR,
                                        number_streams * 3 *
                                        sizeof(cl_double),
                                        fstream_c, &err);
    check_error(err);
    ocl_mem->inv_fstream_sigma_sq2 = clCreateBuffer(ocl_mem->context,
                                     CL_MEM_READ_ONLY |
                                     CL_MEM_COPY_HOST_PTR,
                                     number_streams *
                                     sizeof(cl_double),
                                     inv_fstream_sigma_sq2,
                                     &err);
    check_error(err);
    delete [] fstream_a;
    delete [] fstream_c;
    delete [] inv_fstream_sigma_sq2;
}

void setup_weights(ocl_mem_t* ocl_mem,
                   double* parameters,
                   double* bg_int,
                   double* st_int,
                   int number_streams)
{
    double* stream_weight = new double[number_streams];
    double exp_weight = exp(BACKGROUND_WEIGHT);
    double sum_exp_weights = exp_weight;
    double bg_weight = exp_weight / *bg_int;
    for (int i = 0; i < number_streams; ++i)
    {
        exp_weight = exp(STREAM_WEIGHTS(i));
        sum_exp_weights += exp_weight;
        stream_weight[i] = exp_weight / st_int[i];
    }

    bg_weight = bg_weight / sum_exp_weights;
    printf("bg_weight = %.15f\n", bg_weight);
    for (int i = 0; i < number_streams; ++i)
    {
        stream_weight[i] = stream_weight[i] / sum_exp_weights;
        printf("st_weight[%d] = %.15f\n", i, stream_weight[i]);
    }
    cl_int err;
    ocl_mem->bg_weight = clCreateBuffer(ocl_mem->context,
                                        CL_MEM_READ_ONLY |
                                        CL_MEM_COPY_HOST_PTR,
                                        sizeof(cl_double),
                                        &bg_weight,
                                        &err);
    check_error(err);
    ocl_mem->st_weight = clCreateBuffer(ocl_mem->context,
                                        CL_MEM_READ_ONLY |
                                        CL_MEM_COPY_HOST_PTR,
                                        number_streams *
                                        sizeof(cl_double),
                                        stream_weight,
                                        &err);
    check_error(err);
    delete [] stream_weight;
}

void destruct_ocl(ocl_mem_t* ocl_mem)
{
    for (int i = 0; i < ocl_mem->num_integrals; ++i)
    {
        check_error(clReleaseMemObject(ocl_mem->sinb[i]));
        check_error(clReleaseMemObject(ocl_mem->sinl[i]));
        check_error(clReleaseMemObject(ocl_mem->cosb[i]));
        check_error(clReleaseMemObject(ocl_mem->cosl[i]));
        check_error(clReleaseMemObject(ocl_mem->v[i]));
        check_error(clReleaseMemObject(ocl_mem->r_point[i]));
        check_error(clReleaseMemObject(ocl_mem->qw_r3_N[i]));
        check_error(clReleaseMemObject(ocl_mem->bg_int[i]));
        check_error(clReleaseMemObject(ocl_mem->st_int[i]));
    }
    check_error(clReleaseMemObject(ocl_mem->fstream_a));
    check_error(clReleaseMemObject(ocl_mem->fstream_c));
    check_error(clReleaseMemObject(ocl_mem->inv_fstream_sigma_sq2));
    check_error(clReleaseMemObject(ocl_mem->stars));
    check_error(clReleaseMemObject(ocl_mem->qgaus_W));
    check_error(clReleaseMemObject(ocl_mem->dx));
    check_error(clReleaseMemObject(ocl_mem->probability));
    check_error(clReleaseMemObject(ocl_mem->bg_weight));
    check_error(clReleaseMemObject(ocl_mem->st_weight));
    delete [] ocl_mem->sinb;
    delete [] ocl_mem->sinl;
    delete [] ocl_mem->cosb;
    delete [] ocl_mem->cosl;
    delete [] ocl_mem->v;
    delete [] ocl_mem->r_point;
    delete [] ocl_mem->qw_r3_N;
    delete [] ocl_mem->bg_int;
    delete [] ocl_mem->st_int;
    check_error(clReleaseCommandQueue(ocl_mem->queue));
    check_error(clReleaseContext(ocl_mem->context));
    delete ocl_mem;
}

