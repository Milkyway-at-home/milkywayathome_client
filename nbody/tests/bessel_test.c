/*
 * Copyright (c) 2019 Rensselaer Polytechnic Institute
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
 */

#include "milkyway_util.h"
#include "nbody_bessel.h"
#include <stdio.h> 
#include <stdlib.h>
#include <time.h>

int main()
{
    srand(time(0));

    int test_fails=0;
    int wronsk_fails=0;
    int i0_fails=0;
    int i1_fails=0;
    int k0_fails=0;
    int k1_fails=0;

    const real_0 min_val = 0.5;
    const real_0 max_val = 60.0;
    const real_0 step = 0.0001;
    const int nTests = 100000;
    const real_0 ZERO_THRESHOLD = 1.0e-7;

    real_0 x;
    real_0 wronskian;
    real_0 i0_derv;
    real_0 i1_derv;
    real_0 k0_derv;
    real_0 k1_derv;
    real_0 comp_wron;
    real_0 comp_i0;
    real_0 comp_i1;
    real_0 comp_k0;
    real_0 comp_k1;

    for (int j=0; j<nTests; j++) {
        int fails = 0;
        x = ((real_0)rand()/(real_0)RAND_MAX)*(max_val - min_val) + min_val;
        //mw_printf("X = %.15f\n",x);
        //mw_printf("I0(X) = %.15f\n", besselI0(x));
        //mw_printf("I1(X) = %.15f\n", besselI1(x));
        //mw_printf("K0(X) = %.15f\n", besselK0(x));
        //mw_printf("K1(X) = %.15f\n", besselK1(x));

        wronskian = besselI0(x)*besselK1(x) + besselI1(x)*besselK0(x);
        i0_derv = (besselI0(x+step)-besselI0(x-step))/2.0/step;
        k0_derv = (besselK0(x+step)-besselK0(x-step))/2.0/step;
        i1_derv = (besselI1(x+step)-besselI1(x-step))/2.0/step;
        k1_derv = (besselK1(x+step)-besselK1(x-step))/2.0/step;

        comp_wron = mw_abs_0(wronskian - 1.0/x)*x;
        comp_i0 = mw_abs_0(i0_derv - besselI1(x))/besselI1(x);
        comp_k0 = mw_abs_0(k0_derv + besselK1(x))/besselK1(x);
        comp_i1 = mw_abs_0(i1_derv - besselI0(x) + besselI1(x)/x)/i1_derv;
        comp_k1 = -mw_abs_0(k1_derv + besselK0(x) + besselK1(x)/x)/k1_derv;

        if(comp_wron > ZERO_THRESHOLD) {
            fails += 1;
            wronsk_fails += 1;
            //mw_printf("Wronskian Error = %.15f\n",comp_wron);
        }

        if(comp_i0 > ZERO_THRESHOLD) {
            fails += 1;
            i0_fails += 1;
            //mw_printf("I0 Dervivative Error = %.15f\n",comp_i0);
        }

        if(comp_k0 > ZERO_THRESHOLD) {
            fails += 1;
            k0_fails += 1;
            //mw_printf("K0 Dervivative Error = %.15f\n",comp_k0);
        }

        if(comp_i1 > ZERO_THRESHOLD) {
            fails += 1;
            i1_fails += 1;
            //mw_printf("I1 Dervivative Error = %.15f\n",comp_i1);
        }

        if(comp_k1 > ZERO_THRESHOLD) {
            fails += 1;
            i1_fails += 1;
            //mw_printf("K1 Dervivative Error = %.15f\n",comp_k1);
        }

        if (fails != 0) {
            test_fails += 1;
        }
    }

    if (test_fails != 0) {
        mw_printf("-----------TESTS FAILED-----------\n");
        mw_printf("%d out of %d points failed.\n",test_fails,nTests);
        mw_printf("%d Wronskian tests failed.\n",wronsk_fails);
        mw_printf("%d I0 Derivative tests failed.\n", i0_fails);
        mw_printf("%d I1 Derivative tests failed.\n", i1_fails);
        mw_printf("%d K0 Derivative tests failed.\n", k0_fails);
        mw_printf("%d K1 Derivative tests failed.\n", k1_fails);
    }
    else {
        mw_printf("BESSEL TESTS PASSED!\n");
    }

    return test_fails;
}
