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
#include <cal/cal_il.hpp>
#include <cal/cal_il_math.hpp>
#include <fstream>

#include "separation_cal.h"
#include "separation_cal_kernelgen.h"
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


/* exp which uses fewer registers and is accurate enough. CHECKME: Destroys x? */
static double1 exp_custom(double1 x)
{
    emit_comment("exp");

    double1 tmp1 = x * double1(0x1.71547652b82fep+0);
    double2 dtmp;
    dtmp.y() = fract(tmp1);

    int1 expi = cast_type<int1>(tmp1 - dtmp.y());

    tmp1 = dtmp.y() + double1(-0x1p-1);
    dtmp.x() = tmp1 * tmp1;
    x = mad(dtmp.x(), double1(0x1.52b5c4d1f00b9p-16), double1(0x1.4aa4eb649a98fp-7));
    x = mad(dtmp.x(), x, double1(0x1.62e42fefa39efp-1));

    double2 tmp23;
    tmp23.x() = tmp1 * x;
    tmp23.y() = mad(dtmp.x(), double1(0x1.657cdd06316dcp-23), double1(0x1.3185f478ff1ebp-12));
    tmp23.y() = mad(dtmp.x(), tmp23.y(), double1(0x1.bf3e7389ceff9p-5));
    tmp23.y() = mad(dtmp.x(), tmp23.y(), double1(0x1p+0));

    double1 tmp4 = mad(tmp23.x(), double1(-0x1p-1), tmp23.y());
    tmp1 = tmp23.x() / tmp4;
    tmp1 = mad(tmp1, double1(0x1.6a09e667f3bcdp+0), double1(0x1.6a09e667f3bcdp+0));

    return ldexp(tmp1, expi);
}

static void createSeparationKernelCore(input2d<double2>& bgInput,
                                       std::vector< input2d<double2> >& streamInput,
                                       input2d<double2>& rPts,
                                       input1d<double2>& rConsts,
                                       input2d<double2>& lTrigBuf,
                                       input2d<double2>& bTrigBuf,
                                       const AstronomyParameters* ap,
                                       const IntegralArea* ia,
                                       const StreamConstants* sc)
{
    unsigned int j;
    std::vector<double1> streamIntegrals;
    unsigned int number_streams = 3;  /* FIXME: Temporary to compare against old things */

    indexed_register<double1> sg_dx("cb0");
    named_variable<float1> nu_step("cb1[0].x");
    named_variable<double1> nu_id("cb1[0].zw");
    named_variable<double2> bgOut("o0");
    named_variable<float2> pos("vWinCoord0"); /* .x() = mu, .y() = r */
    std::vector< named_variable<double2> > streamOutputRegisters;

    std::stringstream regName;
    for (j = 0; j < number_streams; ++j)
    {
        regName.seekp(0);
        regName << 'o' << (j + 1);
        streamOutputRegisters.push_back(regName.str());
    }

    double2 lTrig = lTrigBuf(nu_step, pos.x());
    double2 bTrig = bTrigBuf(nu_step, pos.x());

    float1 i = float1(0.0);

    /* 0 integrals and get stream constants */
    double1 bg_int = double1(0.0);
    for (j = 0; j < number_streams; ++j)
        streamIntegrals.push_back(double1(0.0));

    il_while (i < float1(ap->convolve))
    {
        /* CHECKME: could try rPts[float2something] */
        double2 rPt = rPts(pos.y(), i);
        i = i + float1(1.0);

        double1 zp = rPt.x() * bTrig.y();
        double1 x = mad(zp, lTrig.y(), double1(ap->m_sun_r0));
        double1 y = zp * lTrig.x();
        double1 z = rPt.x() * bTrig.x();


        double1 tmp1 = x * x;
        double1 tmp2 = z * z;
        double1 tmp3 = mad(y, y, tmp1);
        double1 tmp4 = mad(double1(ap->q_inv_sqr), tmp2, tmp3);
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
            r_in_mag = sg_dx[i] + gPrime(rcs[pos.y()]);
            tmp = mad(double1(ap->bg_b), r_in_mag, double1(ap->bg_c));
            tmp = mad(double1(ap->bg_a, r_in_mag * r_in_mag, tmp));
            bg_int += tmp * rPt.y();
        }
        #endif

        emit_comment("stream loops");
        for (j = 0; j < number_streams; ++j)
        {
            emit_comment("begin stream");
            double1 xs = x - double1(X(sc[j].c));
            double1 ys = y - double1(Y(sc[j].c));
            double1 zs = z - double1(Z(sc[j].c));

            /* Dot product */
            tmp = double1(X(sc[j].a)) * xs;
            tmp = mad(double1(Y(sc[j].a)), ys.x(), tmp);
            tmp = mad(double1(Z(sc[j].a)), zs.x(), tmp);  /* tmp = dotted */

            xs = mad(-tmp, double1(X(sc[j].a)), xs);
            ys = mad(-tmp, double1(Y(sc[j].a)), ys);
            zs = mad(-tmp, double1(Z(sc[j].a)), zs);

            tmp = xs * xs;
            tmp = mad(ys, ys, tmp);
            tmp = mad(zs, zs, tmp);

            tmp = tmp * double1(-sc[j].sigma_sq2_inv);
            tmp = exp_custom(tmp);

            streamIntegrals[j] = mad(rPt.y(), tmp, streamIntegrals[j]);
        }
    }
    il_endloop

    double1 V_reff_xr_rp3 = nu_id * rConsts(pos.y()).x();

    std::vector<double2> streamRead;
    double2 bgRead = bgInput[pos];
    for (j = 0; j < number_streams; ++j)
        streamRead.push_back(streamInput[j][pos]);

    emit_comment("Output");
    bgOut.x() = bgRead.x() + (V_reff_xr_rp3 * bg_int);

    emit_comment("Stream output");
    //for (j = 0; j < ap->number_streams; ++j)
    for (j = 0; j < number_streams; ++j)
        streamOutputRegisters[j].x() = streamRead[j].x() + (V_reff_xr_rp3 * streamIntegrals[j]);

    streamIntegrals.clear();
}

std::string createSeparationKernel(const AstronomyParameters* ap,
                                   const IntegralArea* ia,
                                   const StreamConstants* sc)
{
    unsigned int numRegSgDx;
    unsigned int sgdxSize;
    unsigned int i;
    std::stringstream code;

    sgdxSize = sizeof(real) * ap->convolve;
    numRegSgDx = sgdxSize / 16;
    if (sgdxSize % 16 != 0)
    {
        warn("sg_dx size not divisible by register size\n");
        return "";
    }

    code << "il_ps_2_0\n";
    code << format("dcl_cb cb0[%u]\n") % numRegSgDx;
    code << "dcl_cb cb1[1]\n";

    for (i = 0; i < ap->number_streams; ++i)
        code << format("dcl_output_usage(generic) o%u.xy\n") % (i + 1);

    code << "dcl_input_position_interp(linear_noperspective) vWinCoord0.xy__\n";

    Source::begin();

    input2d<double2> rPts(0);
    input1d<double2> rConsts(1);
    input2d<double2> lTrig(2);
    input2d<double2> bTrig(3);

    input2d<double2> bgInput(4);
    std::vector< input2d<double2> > streamInputs;

    for (i = 0; i < ap->number_streams; ++i)
        streamInputs.push_back(input2d<double2>(5 + i));

    createSeparationKernelCore(bgInput, streamInputs, rPts, rConsts, lTrig, bTrig, ap, ia, sc);

    Source::end();

    Source::emitHeader(code);
    Source::emitCode(code);

    code << "end\n";

    return code.str();
}

char* separationKernelSrc(const AstronomyParameters* ap,
                          const IntegralArea* ia,
                          const StreamConstants* sc)
{
    std::string src = createSeparationKernel(ap, ia, sc);
     return strdup(src.c_str());
 }


