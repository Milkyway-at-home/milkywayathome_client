/*
 *  Copyright (c) 2020 Eric Mendelsohn
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "milkyway_util.h"
#include "nbody_virial.h"

static real_0 nbDwarfDensity(real_0 r, real_0 a_b, real_0 a_d, real_0 M_b, real_0 M_d)
{
    /** FIXME: More dwarf profiles should be added here. For now, I just have two-component Plummer spheres here. **/
    real_0 b_comp = M_b/mw_pow_0(a_b,3.0) * mw_pow_0(1.0 + r*r/a_b/a_b,-2.5);
    real_0 d_comp = M_d/mw_pow_0(a_d,3.0) * mw_pow_0(1.0 + r*r/a_d/a_d,-2.5);
    return 3.0/4.0/M_PI*(b_comp + d_comp);
}

static real_0 nbDwarfIntegral(real_0 r, real_0 a_b, real_0 a_d, real_0 M_b, real_0 M_d)
{
    /** FIXME: More dwarf profiles should be added here. For now, I just have two-component Plummer spheres here. **/
    real_0 b_comp = M_b/a_b * r * r * mw_pow_0(1.0 + r*r/a_b/a_b,-0.5);
    real_0 d_comp = M_d/a_d * r * r * mw_pow_0(1.0 + r*r/a_d/a_d,-0.5);
    return 1.0/4.0/M_PI*(b_comp + d_comp);
}

real_0 nbCalculateVirial(real_0 a_b, real_0 a_d, real_0 M_b, real_0 M_d) /** General double integral used to calculate Henon length units (SPHERICAL PROFILES ONLY!)**/
{
    real_0 M_tot = M_b + M_d;
    real_0 weight[5];
    weight[0] = 0.2369268850561891;
    weight[1] = 0.4786286704993665;
    weight[2] = 0.5688888888888889;
    weight[3] = 0.4786286704993665;
    weight[4] = 0.2369268850561891;

    real_0 point[5];
    point[0] = -0.9061798459386640;
    point[1] = -0.5384693101056831;
    point[2] = 0.0;
    point[3] = 0.5384693101056831;
    point[4] = 0.9061798459386640;
    const int n = 100000;
    int i1, j1;
    const real_0 lowerlimit = 0.0;
    real_0 upperlimit, a1, b1, r1, integral1;

    //Integral goes to from zero to infinity, so we set upper limit to 1000 times larger scale radius
    if (a_b > a_d)
    {
        upperlimit = 1000.0*a_b;
    }
    else
    {
        upperlimit = 1000.0*a_d;
    }

    const real_0 width = (upperlimit-lowerlimit)/(n*1.0);

    if ((a_b <= 0.0)||(a_d <= 0.0))
    {
        mw_fail("Non-positive scale radius detected!");
    }
    else if ((M_b < 0.0)||(M_d < 0.0))
    {
        mw_fail("Negative mass detected!");
    }
    else if ((M_b == 0.0)&&(M_d == 0.0))
    {
        mw_fail("No dwarf masses detected!");
    }
    else
    {
        integral1 = 0.0;
        for (i1 = 0; i1 < n; i1++)
        {
            a1 = lowerlimit + i1*width;
            b1 = lowerlimit + (i1+1)*width;
            for (j1 = 0; j1 < 5; j1++)
            {
                r1 = width*point[j1]/2.0 + (a1+b1)/2.0;
                integral1 += weight[j1]*nbDwarfDensity(r1, a_b, a_d, M_b, M_d)*nbDwarfIntegral(r1, a_b, a_d, M_b, M_d)*width/2.0;
            }
        }
    }

    real_0 R_v = M_tot*M_tot/16.0/M_PI/M_PI/integral1;
    //mw_printf("[a_b, a_d, M_b, M_d] = [%.15f, %.15f, %.15f, %.15f]\n", a_b, a_d, M_b, M_d);
    //mw_printf("Dwarf Virial Radius = %.15f\n", R_v);
    return R_v;
}

