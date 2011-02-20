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


/* Doesn't give correct inf results when dividing by 0 */
static double1 div_custom(double1 x, double1 y)
{
    float1 yf = cast_type<float1>(y);
    double1 tmp = cast_type<double1>(reciprocal(yf));

    double1 tmp2 = mad(-y, tmp, double1(1.0));
    double1 tmp3 = mad(tmp, tmp2, tmp);
    double1 tmp4 = x * tmp3;
    double1 tmp5 = mad(-y, tmp4, x);

    return mad(tmp5, tmp3, tmp4);
}

static double1 sqrt_custom(double1 x)
{
    emit_comment("sqrt");
    float1 fx = cast_type<float1>(x);
    fx = rsqrt(fx);
    double1 y = cast_type<double1>(-fx);

    double1 tmp = y * y;
    tmp = mad(x, tmp, double1(-3.0));
    y = y * tmp;

    tmp = x * y;
    y = y * tmp;

    y = mad(y, double1(-0.0625), double1(0.75));
    return y * tmp;
}


/* exp which uses fewer registers and is accurate enough. CHECKME: Destroys x? */
/* FIXME: Quality variable names */
static double1 exp_custom(double1 x)
{
    emit_comment("exp");

    double1 xx = x * double1(1.0 / log(2.0));
    double1 xxFract = fract(xx);

    int1 expi = cast_type<int1>(xx - xxFract);

    double1 tmp = xxFract + double1(-0.5);
    double1 dtmp = tmp * tmp;
    double1 y = mad(dtmp, double1(0x1.52b5c4d1f00b9p-16), double1(0x1.4aa4eb649a98fp-7));
    double1 z = tmp * mad(dtmp, y, double1(0x1.62e42fefa39efp-1));

    double1 tmp3 = mad(dtmp, double1(0x1.657cdd06316dcp-23), double1(0x1.3185f478ff1ebp-12));
    tmp3 = mad(dtmp, tmp3, double1(0x1.bf3e7389ceff9p-5));
    tmp3 = mad(dtmp, tmp3, double1(1.0));

    double1 tmp4 = mad(z, double1(-0.5), tmp3);
    double1 divRes = div_custom(z, tmp4);
    divRes = mad(divRes, double1(sqrt(2.0)), double1(sqrt(2.0)));

    return ldexp(divRes, expi);
}

static double2 kahanAdd(double2 kSum, double1 x)
{
    double1 y = x - kSum.y();
    double2 tc;
    tc.x() = kSum.x() + y;
    tc.y() = (tc.x() - kSum.x()) - y;
    return tc;
}


/* Slightly faster if true, but 7 more registers */


/* For reference: ATI application register usage in 1, 2, 3 stream cases: 13, 19, 25 respectively */
/* With R_PT_SWAP = 0: 12, 14, 16 */
/* With R_PT_SWAP = 1: 12, 15, 17 */

static void createSeparationKernelCore(input2d<double2>& bgInput,
                                       std::vector< input2d<double2> >& streamInput,
                                       input2d<double2>& rPts,
                                       input1d<double2>& rConsts,
                                       input2d<double2>& lTrigBuf,
                                       input2d<double1>& bTrigBuf,
                                       const AstronomyParameters* ap,
                                       const IntegralArea* ia,
                                       const StreamConstants* sc)
{
    unsigned int j;
    unsigned int number_streams = 3;  /* FIXME: Temporary to compare against old things */

    indexed_register<double1> sg_dx("cb0");
    named_variable<float1> nu_step("cb1[0].x");
    named_variable<double1> nu_id("cb1[0].zw");
    named_variable<float2> pos("vWinCoord0"); /* .x() = mu, .y() = r */

    named_variable<double2> bgOut("o0");
    std::vector< named_variable<double2> > streamOutputRegisters;

    std::stringstream regName;
    for (j = 0; j < number_streams; ++j)
    {
        regName.seekp(0);
        regName << 'o' << (j + 1);
        streamOutputRegisters.push_back(regName.str());
    }

    double2 lTrig = lTrigBuf(nu_step, pos.x());
    double1 bSin = bTrigBuf(nu_step, pos.x());

    float2 i = float2(0.0, pos.y());

    /* 0 integrals and get stream constants */
    double1 bg_int = double1(0.0);
    std::vector<double1> streamIntegrals;
    for (j = 0; j < number_streams; ++j)
        streamIntegrals.push_back(double1(0.0));

    /* Counting down seems to save 2 registers, but is slightly slower */
    il_while (i.x() < float1((float) ap->convolve))
    {
        double2 rPt = rPts[i.xy()];
        i.x() = i.x() + float1(1.0);

        double1 x = mad(rPt.x(), lTrig.x(), double1(ap->m_sun_r0));
        double1 y = rPt.x() * lTrig.y();
        double1 z = rPt.x() * bSin;

        double1 tmp1 = x * x;
        double1 tmp2 = z * z;
        double1 tmp3 = mad(y, y, tmp1);
        double1 tmp4 = mad(double1(ap->q_inv_sqr), tmp2, tmp3);

        double1 rg = sqrt_custom(tmp4);
        double1 rs = rg + double1(ap->r0);

        emit_comment("bg_int increment");
        bg_int += div_custom(rPt.y(), rg * rs * rs * rs);

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
            double1 dotted = double1(X(sc[j].a)) * xs;
            dotted = mad(double1(Y(sc[j].a)), ys, dotted);
            dotted = mad(double1(Z(sc[j].a)), zs, dotted);

            double1 xsp = mad(dotted, double1(-X(sc[j].a)), xs);
            double1 ysp = mad(dotted, double1(-Y(sc[j].a)), ys);
            double1 zsp = mad(dotted, double1(-Z(sc[j].a)), zs);

            double1 sqrv = xsp * xsp;
            sqrv = mad(ysp, ysp, sqrv);
            sqrv = mad(zsp, zsp, sqrv);

            double1 streamExp = exp_custom(sqrv * double1(-sc[j].sigma_sq2_inv));
            emit_comment("End exp");

            streamIntegrals[j] = mad(rPt.y(), streamExp, streamIntegrals[j]);
        }
        emit_comment("End streams");
    }
    il_endloop

    double1 V_reff_xr_rp3 = nu_id * rConsts(pos.y()).x();
    std::vector<double2> streamRead;
    double2 bgRead = bgInput[pos];
    for (j = 0; j < number_streams; ++j)
        streamRead.push_back(streamInput[j][pos]);

    /* Put these multiplies together */
    bg_int *= V_reff_xr_rp3;
    for (j = 0; j < number_streams; ++j)
        streamIntegrals[j] *= V_reff_xr_rp3;

#define KAHAN 0

#if KAHAN
    bgOut = kahanAdd(bgRead, bg_int);
    for (j = 0; j < number_streams; ++j)
        streamOutputRegisters[j] = kahanAdd(streamRead[j], streamIntegrals[j]);
#else
    bg_int += bgRead.x();
    for (j = 0; j < number_streams; ++j)
        streamIntegrals[j] += streamRead[j].x();

    emit_comment("Output");
    bgOut.x() = bg_int;

    emit_comment("Stream output");
    //for (j = 0; j < ap->number_streams; ++j)
    for (j = 0; j < number_streams; ++j)
        streamOutputRegisters[j].x() = streamIntegrals[j];
#endif /* KAHAN */

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

    code << "dcl_output_usage(generic) o0\n";
    for (i = 0; i < ap->number_streams; ++i)
        code << format("dcl_output_usage(generic) o%u\n") % (i + 1);

    code << "dcl_input_position_interp(linear_noperspective) vWinCoord0.xy__\n";

    Source::begin();

    input2d<double2> rPts(0);
    input1d<double2> rConsts(1);
    input2d<double2> lTrig(2);
    input2d<double1> bTrig(3);

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


