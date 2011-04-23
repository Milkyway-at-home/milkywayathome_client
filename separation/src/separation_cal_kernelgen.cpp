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

#ifdef _MSC_VER
  #pragma warning( disable : 4522 )
#endif

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

static double1 log10_custom(double1 x)
{
    return log(x) / log(10.0);
}

/* Doesn't give correct inf results when dividing by 0 */
static double1 div_custom(double1 x, double1 y)
{
    double1 result;

    il_func(_out(result), x, y)
    {
        emit_comment("div_custom");
        //float1 yf = cast_type<float1>(y);
        //double1 tmp = cast_type<double1>(reciprocal(yf));

        double1 tmp = native_reciprocal(y);

        double1 tmp2 = mad(-y, tmp, 1.0);
        double1 tmp3 = mad(tmp, tmp2, tmp);
        double1 tmp4 = x * tmp3;
        double1 tmp5 = mad(-y, tmp4, x);
        result = mad(tmp5, tmp3, tmp4);
    }
    il_endfunc

    return result;
}

static double1 sqrt_custom(double1 x)
{
    double1 result;

    il_func(_out(result), x)
    {
        emit_comment("sqrt_custom");

        float1 fx = cast_type<float1>(x);
        double1 y = cast_type<double1>(-native_rsqrt(fx));

        double1 tmp = y * y;
        tmp = mad(x, tmp, -3.0);
        y *= tmp;

        tmp = x * y;
        y *= tmp;

        result = tmp * mad(y, -0.0625, 0.75);
    }
    il_endfunc

    return result;
}

static double1 pow_custom(double1 a, double1 b)
{
    /* TODO: Should really be using 2 instead of e, but that's more
     * work right now since exp and log exist already in CAL++. */
    return exp(log(a) * b);
}

static double1 exp10_custom(double1 x)
{
    return pow_custom(double1(10.0), x);
}

/* FIXME: Quality variable names */
static double1 exp_custom(double1 x)
{
    double1 result;

    il_func(_out(result), x)
    {
        emit_comment("exp_custom");

        double1 xx = x * (1.0 / log(2.0));
        double1 xxFract = fract(xx);
        int1 expi = cast_type<int1>(xx - xxFract);

        double1 tmp = xxFract - 0.5;
        double1 dtmp = tmp * tmp;
        double1 y = mad(dtmp, 0.000020188691287, 0.010090460718136);
        double1 z = tmp * mad(dtmp, y, 0.693147180559945);

        double1 tmp3 = mad(dtmp, 0.000000166468205, 0.000291369687659);
        tmp3 = mad(dtmp, tmp3, 0.054595208798190);
        tmp3 = mad(dtmp, tmp3, 1.0);

        double1 tmp4 = mad(z, -0.5, tmp3);
        double1 divRes = div_custom(z, tmp4);
        divRes = mad(divRes, sqrt(2.0), sqrt(2.0));

        result = ldexp(divRes, expi);
    }
    il_endfunc

    return result;
}

static double1 calcReffXrRp3(const double1 coords, const double1 gPrime)
{
    const real sigmoid_curve_params[3] = { 0.9402, 1.6171, 23.5877 };

    /* REFF */
    double1 exp_result = exp_custom(double1(sigmoid_curve_params[1]) * (gPrime - double1(sigmoid_curve_params[2])));
    double1 reff_value = double1(sigmoid_curve_params[0]) / (exp_result + 1.0);
    double1 rPrime3 = cube(coords);
    double1 reff_xr_rp3 = reff_value * xr / rPrime3;
    return reff_xr_rp3;
}

static double1 calcG(const double1 coords)
{
    return 5.0 * (log10_custom(1000.0 * coords) - 1.0) + absm;
}

static double2 calcRPoint(const AstronomyParameters* ap, const double1 dx, const double1 qgaus_W, const double1 gPrime)
{
    double2 r_pt;

    double1 g = gPrime + dx;

    /* MAG2R */
    r_pt.x() = 0.001 * exp10_custom(0.2 * (g - absm) + 1.0);

    double1 r3 = cube(r_pt.x());
    double1 exponent = sqr(g - gPrime) * inv(2.0 * sqr(stdev));
    double1 N = ap->coeff * exp(-exponent);
    r_pt.y() = qgaus_W * r3 * N;

    return r_pt;
}


static double2 kahanAdd(double2 kSum, double1 x)
{
    double1 y = x - kSum.y();
    double2 tc;
    tc.x() = kSum.x() + y;
    tc.y() = (tc.x() - kSum.x()) - y;
    return tc;
}


/* For reference: ATI application register usage in 1, 2, 3 stream cases: 13, 19, 25 respectively */
static void createSeparationKernelCore(const AstronomyParameters* ap,
                                       const StreamConstants* sc,

                                       bool likelihoodKernel,  /* integral kernel or likelihood? */

                                       /* Arguments used by integral kernel */
                                       input2d<double2>* rPts,
                                       input1d<double2>* rConsts,
                                       input2d<double2>* lTrigBuf,
                                       input2d<double1>* bTrigBuf,
                                       named_variable<float1>* nuStep,
                                       named_variable<double1>* nu_id,

                                       /* Arguments used by likelihood kernel */
                                       input1d<double2>* starsXYBuf,
                                       input1d<double1>* starsZBuf)
{
    unsigned int j;

    input2d<double2> bgInput(4);
    std::vector< input2d<double2> > streamInputs;

    for (j = 0; j < ap->number_streams; ++j)
        streamInputs.push_back(input2d<double2>(5 + j));

    indexed_register<double1> sg_dx("cb0");
    indexed_register<double1> sg_qgaus_W("cb1");
    named_variable<float2> pos("vWinCoord0"); /* .x() = r, .y() = mu in integral kernel */

    named_variable<double2> bgOut("o0");
    std::vector< named_variable<double2> > streamOutputRegisters;

    for (j = 0; j < ap->number_streams; ++j)
        streamOutputRegisters.push_back(str(format("o%u") % (j + 1)));

    double2 lTrig;
    double1 bSin;

    double2 starXY;
    double1 starZ;

    if (!likelihoodKernel)
    {
        lTrig = (*lTrigBuf)(*nuStep, pos.y());
        bSin = (*bTrigBuf)(*nuStep, pos.y());
    }
    else
    {
        starXY = (*starsXYBuf)[pos.x()];
        starZ = (*starsZBuf)[pos.x()];

        //double1 lSin = sin(
        //lTrig = starXY.x() /* l */
    }

    float2 i = float2(pos.x(), 0.0);

    double1 gPrime = calcG(starZ);

    /* 0 integrals and get stream constants */
    double1 bg_int = double1(0.0);
    std::vector<double1> streamIntegrals;
    for (j = 0; j < ap->number_streams; ++j)
        streamIntegrals.push_back(double1(0.0));

    il_whileloop
    {
        double2 rPt;

        if (!likelihoodKernel)
        {
            rPt = (*rPts)[i];
        }
        else
        {
            int1 idx = cast_type<int1>(i.y());

            /* calc_r_point */
            rPt = calcRPoint(ap, sg_dx[idx], sg_qgaus_W[idx], gPrime);
        }

        double1 x = mad(rPt.x(), lTrig.x(), ap->m_sun_r0);
        double1 y = rPt.x() * lTrig.y();
        double1 z = rPt.x() * bSin;

        double1 tmp1 = x * x;
        double1 tmp3 = mad(y, y, tmp1);
        double1 tmp2 = z * z;
        double1 tmp4 = mad(ap->q_inv_sqr, tmp2, tmp3);

        double1 rg = sqrt_custom(tmp4);
        double1 rs = rg + ap->r0;

        bg_int += div_custom(rPt.y(), rg * rs * rs * rs);

        if (ap->aux_bg_profile)
        {
            /* TODO: Checkme */
            double1 r_in_mag = sg_dx[cast_type<int1>(i.y())] + (*rConsts)[pos.x()].y();
            double1 tmp = mad(ap->bg_b, r_in_mag, ap->bg_c);
            tmp = mad(ap->bg_a, r_in_mag * r_in_mag, tmp);
            bg_int = mad(tmp, rPt.y(), bg_int);
        }

        emit_comment("stream loops");
        for (j = 0; j < ap->number_streams; ++j)
        {
            emit_comment("begin stream");
            double1 xs = x - X(sc[j].c);
            double1 ys = y - Y(sc[j].c);
            double1 zs = z - Z(sc[j].c);

            /* Dot product */
            double1 dotted = X(sc[j].a) * xs;
            dotted = mad(Y(sc[j].a), ys, dotted);
            dotted = mad(Z(sc[j].a), zs, dotted);

            xs = mad(dotted, -X(sc[j].a), xs);
            ys = mad(dotted, -Y(sc[j].a), ys);
            zs = mad(dotted, -Z(sc[j].a), zs);

            double1 sqrv = xs * xs;
            sqrv = mad(ys, ys, sqrv);
            sqrv = mad(zs, zs, sqrv);

            sqrv *= -sc[j].sigma_sq2_inv;

            streamIntegrals[j] = mad(rPt.y(), exp_custom(sqrv), streamIntegrals[j]);
        }
        emit_comment("End streams");

        i.y() = i.y() + 1.0f;
        il_breakc(i.y() >= float1((float) ap->convolve));
    }
    il_endloop

    double1 V_reff_xr_rp3;
    if (!likelihoodKernel)
        V_reff_xr_rp3 = (*nu_id) * (*rConsts)(pos.x()).x();
    else
        V_reff_xr_rp3 = calcReffXrRp3(starZ, gPrime);

    std::vector<double2> streamRead;
    double2 bgRead = bgInput[pos];
    for (j = 0; j < ap->number_streams; ++j)
        streamRead.push_back(streamInputs[j][pos]);

    /* Put these multiplies together */
    bg_int *= V_reff_xr_rp3;
    for (j = 0; j < ap->number_streams; ++j)
        streamIntegrals[j] *= V_reff_xr_rp3;

#define KAHAN 0

#if KAHAN
    bgOut = kahanAdd(bgRead, bg_int);
    for (j = 0; j < ap->number_streams; ++j)
        streamOutputRegisters[j] = kahanAdd(streamRead[j], streamIntegrals[j]);
#else
    bg_int += bgRead.x();
    for (j = 0; j < ap->number_streams; ++j)
        streamIntegrals[j] += streamRead[j].x();

    emit_comment("Output");
    bgOut.x() = bg_int;

    emit_comment("Stream output");
    for (j = 0; j < ap->number_streams; ++j)
        streamOutputRegisters[j].x() = streamIntegrals[j];
#endif /* KAHAN */

    streamIntegrals.clear();
}

static int separationKernelHeader(const AstronomyParameters* ap,
                                  std::stringstream& code)
{
    CALuint i, numRegSgDx, sgdxSize;

    sgdxSize = sizeof(real) * ap->convolve;
    numRegSgDx = sgdxSize / 16;
    if (sgdxSize % 16 != 0)
    {
        warn("sg_dx size not divisible by register size\n");
        return 1;
    }

    code << "il_ps_2_0\n";
    code << "dcl_cb cb0[1]\n";
    code << format("dcl_cb cb1[%u]\n") % numRegSgDx;


    code << "dcl_output_usage(generic) o0\n";
    for (i = 0; i < ap->number_streams; ++i)
        code << format("dcl_output_usage(generic) o%u\n") % (i + 1);

    code << "dcl_input_position_interp(linear_noperspective) vWinCoord0.xy__\n";

    return 0;
}

std::string createSeparationIntegralKernel(const AstronomyParameters* ap,
                                           const StreamConstants* sc,
                                           CALuint device)
{
    std::stringstream code;

    if (separationKernelHeader(ap, code))
        return "";

    Source::begin(device);

    input2d<double2> rPts(0);
    input1d<double2> rConsts(1);
    input2d<double2> lTrig(2);
    input2d<double1> bTrig(3);

    named_variable<float1> nuStep("cb0[0].x");
    named_variable<double1> nuId("cb0[0].zw");

    createSeparationKernelCore(ap, sc, false,
                               &rPts, &rConsts, &lTrig, &bTrig, &nuStep, &nuId,
                               NULL, NULL);

    Source::end();

    Source::emitHeader(code);
    Source::emitCode(code);

    return code.str();
}


std::string createSeparationLikelihoodKernel(const AstronomyParameters* ap,
                                             const StreamConstants* sc,
                                             CALuint device)
{
    CALuint numRegSgDx, sgdxSize;
    std::stringstream code;

    sgdxSize = sizeof(real) * ap->convolve;
    numRegSgDx = sgdxSize / 16;

    if (separationKernelHeader(ap, code))
        return "";

    code << format("dcl_cb cb2[%u]\n") % numRegSgDx;

    Source::begin(device);

    input1d<double2> starsXY(0);
    input1d<double1> starsZ(1);

    createSeparationKernelCore(ap, sc, true,
                               NULL, NULL, NULL, NULL, NULL, NULL,
                               &starsXY, &starsZ);

    Source::end();

    Source::emitHeader(code);
    Source::emitCode(code);

    return code.str();
}

char* separationIntegralKernelSrc(const AstronomyParameters* ap,
                                  const StreamConstants* sc,
                                  CALuint device)
{
    std::string src = createSeparationIntegralKernel(ap, sc, device);
    return strdup(src.c_str());
 }


char* separationLikelihoodKernelSrc(const AstronomyParameters* ap,
                                    const StreamConstants* sc,
                                    CALuint device)
{
    std::string src = createSeparationLikelihoodKernel(ap, sc, device);
    return strdup(src.c_str());
 }


