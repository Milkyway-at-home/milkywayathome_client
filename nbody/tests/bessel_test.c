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

    real x, xp, xm;
    real wronskian;
    real i0_derv;
    real i1_derv;
    real k0_derv;
    real k1_derv;
    real comp_wron;

    real_0 comp_i0;
    real_0 comp_i1;
    real_0 comp_k0;
    real_0 comp_k1;

    real_0 bI0;
    real_0 bK0;
    real_0 bI1;
    real_0 bK1;

    real tmp1, tmp2, tmp3;

    for (int j=0; j<nTests; j++) {
        int fails = 0;
        x = mw_real_const(((real_0)rand()/(real_0)RAND_MAX)*(max_val - min_val) + min_val);
        xp = mw_add_s(&x, step);
        xm = mw_add_s(&x, -step);
        
        //mw_printf("X = %.15f\n",x);
        //mw_printf("I0(X) = %.15f\n", besselI0(x));
        //mw_printf("I1(X) = %.15f\n", besselI1(x));
        //mw_printf("K0(X) = %.15f\n", besselK0(x));
        //mw_printf("K1(X) = %.15f\n", besselK1(x));

        tmp1 = besselI0(&x);
        tmp2 = besselK1(&x);
        tmp1 = mw_mul(&tmp1, &tmp2);
        tmp2 = besselI1(&x);
        tmp3 = besselK0(&x);
        tmp2 = mw_mul(&tmp2, &tmp3);
        wronskian = mw_add(&tmp1, &tmp2);

        tmp1 = besselI0(&xp);
        tmp2 = besselI0(&xm);
        tmp1 = mw_sub(&tmp1, &tmp2);
        i0_derv = mw_mul_s(&tmp1, 0.5/step);

        tmp1 = besselK0(&xp);
        tmp2 = besselK0(&xm);
        tmp1 = mw_sub(&tmp1, &tmp2);
        k0_derv = mw_mul_s(&tmp1, 0.5/step);

        tmp1 = besselI1(&xp);
        tmp2 = besselI1(&xm);
        tmp1 = mw_sub(&tmp1, &tmp2);
        i1_derv = mw_mul_s(&tmp1, 0.5/step);

        tmp1 = besselK1(&xp);
        tmp2 = besselK1(&xm);
        tmp1 = mw_sub(&tmp1, &tmp2);
        k1_derv = mw_mul_s(&tmp1, 0.5/step);

        tmp1 = inv(&x);
        tmp1 = mw_sub(&wronskian, &tmp1);
        tmp1 = mw_abs(&tmp1);
        comp_wron = mw_mul(&tmp1, &x);

        tmp1 = besselI0(&x);
        bI0 = showRealValue(&tmp1);
        tmp1 = besselK0(&x);
        bK0 = showRealValue(&tmp1);
        tmp1 = besselI1(&x);
        bI1 = showRealValue(&tmp1);
        tmp1 = besselK1(&x);
        bK1 = showRealValue(&tmp1);

        comp_i0 = mw_abs_0(showRealValue(&i0_derv) - bI1)/bI1;
        comp_k0 = mw_abs_0(showRealValue(&k0_derv) + bK1)/bK1;
        comp_i1 = mw_abs_0(showRealValue(&i1_derv) - bI0 + bI1/showRealValue(&x))/showRealValue(&i1_derv);
        comp_k1 = -mw_abs_0(showRealValue(&k1_derv) + bK0 + bK1/showRealValue(&x))/showRealValue(&k1_derv);

        if(showRealValue(&comp_wron) > ZERO_THRESHOLD) {
            fails += 1;
            wronsk_fails += 1;
            //mw_printf("Wronskian Error = %.15f\n",showRealValue(&comp_wron));
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
