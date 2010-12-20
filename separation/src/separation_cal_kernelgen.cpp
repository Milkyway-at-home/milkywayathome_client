/*
Copyright (C) 2010  Matthew Arsenault

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

#include <boost/format.hpp>
#include <cal/cal.hpp>
#include <cal/il/cal_il.hpp>
#include <cal/il/cal_il_math.hpp>
#include <cal/il/cal_il_functors_math.hpp>

#include "separation_cal.h"
#include "calculated_constants.h"
#include "r_points.h"

using namespace boost;
using namespace cal;
using namespace cal::il;

static double1 sqrt_custom(double1 x)
{
    emit_comment("sqrt");
    float1 fx = cast_type<float1>(x);
    double1 y = cast_type<double1>(-rsqrt(fx));

    double1 tmp = y * y;
    tmp = mad(x, tmp, double1(-0x1.8p+1));
    y = y * tmp;

    tmp = x * y;
    y = y * tmp;

    y = mad(y,  double1(-0x1p-4), double1(0x1.8p-1));
    return y * tmp;
}


/* exp which uses fewer registers and is accurate enough */
static double1 exp_custom(double1 x)
{
    emit_comment("exp");
    double2 dtmp;
    double1 tmp1 = x * double1(0x1.71547652b82fep+0);
    dtmp.y() = fract(tmp1);

    int1 expi = cast_type<int1>(dtmp.y() - tmp1);

    tmp1 = dtmp.y() + double1(-0x1p-1);
    dtmp.x() = tmp1 * tmp1;
    double1 tmp2 = mad(dtmp.x(), double1(0x1.52b5c4d1f00b9p-16), double1(0x1.4aa4eb649a98fp-7));
    tmp2 = mad(dtmp.x(), tmp2, double1(0x1.62e42fefa39efp-1));

    double2 tmp23;
    tmp23.x() = tmp1 * tmp2;
    tmp23.y() = mad(dtmp.x(), double1(0x1.657cdd06316dcp-23), double1(0x1.3185f478ff1ebp-12));
    tmp23.y() = mad(dtmp.x(), tmp23.y(), double1(0x1.bf3e7389ceff9p-5));
    tmp23.y() = mad(dtmp.x(), tmp23.y(), double1(0x1p+0));

    double1 tmp4 = mad(tmp23.x(), double1(-0x1p-1), tmp23.y());
    tmp1 = tmp23.x() / tmp4;
    tmp1 = mad(tmp1, double1(0x1.6a09e667f3bcdp+0), double1(0x1.6a09e667f3bcdp+0));

    return ldexp(tmp1, expi);
}

static void createSeparationKernelCore(global<double2>& bgOut,
                                       global<double2>& streamOut,
                                       input2d<double2>& rPts,
                                       input1d<double2>& rConsts,
                                       input2d<double2>& lbts,
                                       const AstronomyParameters* ap,
                                       const IntegralArea* ia,
                                       const StreamConstants* sc)
{
    unsigned int j;
    std::vector<double1> streamIntegrals;
    unsigned int number_streams = 3;  /* FIXME: Temporary to compare against old things */

    named_variable<float1> nu_step("cb7[0].x");

    double2 lTrig = lbts(cast_type<float1>(get_global_id(0)), nu_step);
    double2 bTrig = lbts(nu_step, nu_step); // FIXME: Where is this actually?

    float1 i = float1(0.0);
    float1 r_step = cast_type<float1>(get_global_id(1));

    /* 0 integrals and get stream constants */
    double1 bg_int = double1(0.0);
    for (j = 0; j < number_streams; ++j)
        streamIntegrals.push_back(double1(0.0));

    il_while (i < float1(ap->convolve))
    {
        /* CHECKME: could try rPts[float2something] */
        double2 rPt = rPts(r_step, i);
        i = i + float1(1.0);

        double1 zp = rPt.x() * bTrig.y();
        double1 x = mad(zp, lTrig.y(), double1(ap->m_sun_r0));
        double1 y = zp * lTrig.x();
        double1 z = rPt.x() * bTrig.x();

        double1 tmp1 = x * x;
        double1 tmp3 = z * z;
        double1 tmp2 = mad(y, y, tmp1);
        double1 tmp4 = mad(tmp3, double1(ap->q_inv_sqr), tmp2);
        double1 rg = sqrt_custom(tmp4);
        double1 rs = rg + double1(ap->r0);

        double1 tmp = rg * rs;
        tmp = tmp * rs;
        tmp = tmp * rs;

        emit_comment("bg_int increment");
        bg_int += rPt.y() / tmp;

        #if 0
        if (ap->aux_bg_profile)
        {
            r_in_mag = sg_dx[i] + gPrime(rcs[r_step]);
            tmp = mad(double1(ap->bg_b), r_in_mag, double1(ap->bg_c));
            tmp = mad(double1(ap->bg_a, r_in_mag * r_in_mag, tmp));
            bg_int += tmp * rPt.y();
        }
        #endif

        emit_comment("stream loops");
        for (j = 0; j < number_streams; ++j)
        {
            emit_comment("begin stream");
            double2 xs, ys, zs;
            xs.x() = x - double1(X(sc[j].c));
            ys.x() = y - double1(Y(sc[j].c));
            zs.x() = z - double1(Z(sc[j].c));

            /* Dot product */

            /* Juggle the second components of these to prevent usage
             * of 2 extra temp registers. */
            zs.y() = xs.x() * double1(X(sc[j].a));
            xs.y() = mad(double1(Y(sc[j].a)), ys.x(), zs.y());
            ys.y() = mad(double1(Z(sc[j].a)), zs.x(), xs.y());  /* tmp = dotted */

            xs.y() = mad(-ys.y(), double1(X(sc[j].a)), xs.x());
            ys.y() = mad(-ys.y(), double1(Y(sc[j].a)), ys.x());
            zs.y() = mad(-ys.y(), double1(Z(sc[j].a)), zs.x());

            zs.x() = xs.y() * xs.y();
            xs.x() = mad(ys.y(), ys.y(), zs.x());
            ys.x() = mad(zs.y(), zs.y(), xs.x());

            xs.x() = ys.x() * double1(sc[j].sigma_sq2_inv);
            ys.x() = exp_custom(xs.x());

            streamIntegrals[j] = mad(rPt.y(), ys.x(), streamIntegrals[j]);
        }
    }
    il_endloop

    named_variable<double1> nu_id("cb7[0].wz");
    double1 V_reff_xr_rp3 = nu_id * rConsts(r_step).x();

    emit_comment("Output index calculation");
    uint1 idx = get_global_id(0) * uint1(ia->r_steps) + get_global_id(1);

    /* FIXME: Buffer size and read register requirement */

    double1 tmp = bgOut[idx].x();
    bgOut[idx].x() = tmp + (V_reff_xr_rp3 * bg_int);

    emit_comment("Stream output");
    //for (j = 0; j < ap->number_streams; ++j)
    for (j = 0; j < number_streams; ++j)
    {
        tmp = streamOut[idx + j].x();
        streamOut[idx + j].x() = tmp + (V_reff_xr_rp3 * streamIntegrals[j]);
    }


    streamIntegrals.clear();
}

std::string createSeparationKernel(const AstronomyParameters* ap,
                                   const IntegralArea* ia,
                                   const StreamConstants* sc)
{
    std::stringstream code;

    code << "il_ps_2_0\n";
    code << "dcl_cb cb0[1]\n";
    code << "dcl_cb cb7[1]\n";

    code << "dcl_input_position_interp(linear_noperspective) vWinCoord0.xy__\n";

    Source::begin();

    global<double2> bgOut;
    global<double2> streamOut;
    input2d<double2> rPts(0);
    input2d<double2> lbts(1);  /* Actually should be double4 */
    input1d<double2> rConsts(2);

    createSeparationKernelCore(bgOut, streamOut, rPts, rConsts, lbts, ap, ia, sc);

    Source::end();

    Source::emitHeader(code);
    Source::emitCode(code);

    code << "end\n";

    return code.str();
}

const char* separationKernelSrc(const AstronomyParameters* ap,
                                const IntegralArea* ia,
                                const StreamConstants* sc)
{
    return createSeparationKernel(ap, ia, sc).c_str();
}


