/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include <math.h>
#include <stdlib.h>

#include "pi_constants.h"
#include "gauss_legendre.h"
#include "r_constants.h"

#include "parameters.h"

void cpu__reff_gPrime(int r_steps,
                      double r_min,
                      double r_step_size,
                      double mu_step_size,
                      double* cpu__gPrime,
                      double* cpu__reff_xr_rp3_irv)
{
    int i;
    double log_r, r, next_r, coords, gp, exp_result, reff_value, rPrime3, irv;

    for (i = 0; i < r_steps; i++)
    {
        log_r = r_min + (i * r_step_size);
        r = pow(10.0, (log_r - 14.2) / 5.0);
        next_r = pow(10.0, (log_r - 14.2 + r_step_size) / 5.0);
        coords = (next_r + r) / 2.0;
        gp = 5.0 * (log10(coords * 1000) - 1.0) + d_absm;

        irv = (((next_r * next_r * next_r) - (r * r * r)) / 3.0) * mu_step_size * D_DEG2RAD;

        cpu__gPrime[i] = gp;

        exp_result = exp(sigmoid_curve_1 * (gp - sigmoid_curve_2));
        reff_value = sigmoid_curve_0 / (exp_result + 1);
        rPrime3 = coords * coords * coords;
        cpu__reff_xr_rp3_irv[i] = irv * reff_value * d_xr / rPrime3;
    }
}

void cpu__r_qw(int r_steps, int n_convolve, double coeff, double* dx, double* qgaus_W, double* gPrime, double* r_constants)
{
    int i, j, position;
    double g, rp, r3, exponent, N;

    for (i = 0; i < r_steps; i++)
    {
        for (j = 0; j < n_convolve; j++)
        {
            position = ((i * n_convolve) + j) * 2;

            g = gPrime[i] + dx[j];
            rp = pow(10.0, (g - d_absm) / 5.0 + 1.0) / 1000.0;

            r_constants[position] = rp;

            r3 = rp * rp * rp;
            exponent = (dx[j] * dx[j]) / (2 * d_stdev * d_stdev);
            N = coeff * exp(-exponent);

            r_constants[position + 1] = qgaus_W[j] * r3 * N;
        }
    }
}

void cpu__reff_V(int nu_steps, double nu_min, double nu_step_size, int r_steps, double* reff_xr_rp3_irv, double* V)
{
    int i, j;

    double nu, ids;
    for (i = 0; i < nu_steps; i++)
    {
        nu = nu_min + (i * nu_step_size);
        ids = cos((90 - nu - nu_step_size) * D_DEG2RAD) - cos((90 - nu) * D_DEG2RAD);

        for (j = 0; j < r_steps; j++)
        {
            V[(i * r_steps) + j] = reff_xr_rp3_irv[j] * ids;
        }
    }
}

void cpu__r_constants(  int n_convolve,
                        int r_steps, double r_min, double r_step_size,
                        int mu_steps, double mu_min, double mu_step_size,
                        int nu_steps, double nu_min, double nu_step_size,
                        double** cpu__V, double** cpu__r_consts)
{
    int i;
    double* cpu__gPrime, *cpu__reff_xr_rp3_irv;

    cpu__reff_xr_rp3_irv = (double*)malloc(r_steps * sizeof(double));
    cpu__gPrime = (double*)malloc(r_steps * sizeof(double));

    cpu__reff_gPrime(r_steps, r_min, r_step_size, mu_step_size, cpu__gPrime, cpu__reff_xr_rp3_irv);

    *cpu__V = (double*)malloc(nu_steps * r_steps * sizeof(double));
    cpu__reff_V(nu_steps, nu_min, nu_step_size, r_steps, cpu__reff_xr_rp3_irv, *cpu__V);

    *cpu__r_consts = (double*)malloc(2 * r_steps * n_convolve * sizeof(double));

    double* host__qgaus_W = (double*)malloc(n_convolve * sizeof(double));
    double* host__qgaus_X = (double*)malloc(n_convolve * sizeof(double));
    double* host__dx = (double*)malloc(n_convolve * sizeof(double));

    d_gauss_legendre(-1.0, 1.0, host__qgaus_X, host__qgaus_W, n_convolve);
    for (i = 0; i < n_convolve; i++)
    {
        host__dx[i] = 3 * d_stdev * host__qgaus_X[i];
    }

    double coeff = 1.0 / (d_stdev * sqrt(2.0 * D_PI));
    cpu__r_qw(r_steps, n_convolve, coeff, host__dx, host__qgaus_W, cpu__gPrime, *cpu__r_consts);

    free(cpu__reff_xr_rp3_irv);
    free(cpu__gPrime);
    free(host__qgaus_W);
    free(host__qgaus_X);
    free(host__dx);
}


void cpu__r_constants(int n_convolve, INTEGRAL* integral, double** cpu__V, double** cpu__r_consts)
{
    cpu__r_constants(   n_convolve,
                        integral->r_steps, integral->r_min, integral->r_step_size,
                        integral->mu_steps, integral->mu_min, integral->mu_step_size,
                        integral->nu_steps, integral->nu_min, integral->nu_step_size,
                        cpu__V, cpu__r_consts);
}
