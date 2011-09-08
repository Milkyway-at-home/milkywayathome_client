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

#define USE_KAHAN 0

/* Create a flexible kernel that can handle 0 .. IL_MAX_STREAMS, but is ~10% slower (at 3 streams) */
#define FLEXIBLE_KERNEL 0
#define IL_MAX_STREAMS 4
#define IL_MAX_CONVOLVE 256



using namespace boost;
using namespace cal;
using namespace cal::il;

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

static double2 kahanAdd(double2 kSum, double1 x)
{
    double1 y = x - kSum.y();
    double2 tc;
    tc.x() = kSum.x() + y;
    tc.y() = (tc.x() - kSum.x()) - y;
    return tc;
}

/* warning: sizeof(double2) etc. will not work as expected! */
#if USE_KAHAN
typedef double2 SumType;
#define SIZEOF_SUM_TYPE uint1(16)
#else
typedef double1 SumType;
#define SIZEOF_SUM_TYPE uint1(8)
#endif

static uint1 get_global_id0_with_offset()
{
    named_variable<uint1> offset0("cb0[6].x");
    return get_global_id(0) + offset0;
}


/* For reference: ATI application register usage in 1, 2, 3 stream cases: 13, 19, 25 respectively */
static void createSeparationKernelCore(CALuint maxStreams, CALuint uavid)
{
    CALuint j;

/* Variables/kernel arguments etc. */
    named_variable<uint1> extra("cb1[9].x");
    named_variable<uint1> rSteps("cb1[10].x");
    named_variable<uint1> muSteps("cb1[11].x");
    named_variable<uint1> nuSteps("cb1[12].x");

    named_variable<double1> nuId("cb1[13].xy");
    named_variable<uint1> nuStep("cb1[14].x");

    named_variable<uint1> _nConvolve("cb3[1].x");
    named_variable<uint1> _nStream("cb3[1].y");

    named_variable<double1> m_sun_r0("cb3[2].xy");
    named_variable<double1> r0("cb3[2].zw");
    named_variable<double1> q_inv_sqr("cb3[3].xy");

/* Buffer base addresses */
    named_variable<uint1> bgOutBase("cb1[0].x");
    named_variable<uint1> streamsOutBase("cb1[1].x");
    named_variable<uint1> rConstsBase("cb1[2].x");
    named_variable<uint1> rPtsBase("cb1[3].x");
    named_variable<uint1> lTrigBase("cb1[4].x");
    named_variable<uint1> bSinBase("cb1[5].x");

    const uint1 nConvolve = _nConvolve;
    const uint1 nStream = _nStream;

    uav_raw<SumType> bgOut(uavid);
    uav_raw<SumType> streamsOut(uavid);

    uav_raw<double2> rPts(uavid);
    uav_raw<double2> rConsts(uavid);
    uav_raw<double2> lTrigBuf(uavid);
    uav_raw<double1> bTrigBuf(uavid);

    const uint1 gid = get_global_id0_with_offset() - extra;
    const uint1 muStep = gid % muSteps;
    const uint1 rStep = gid / muSteps;

    const uint1 sizeOfDouble2 = uint1(2 * sizeof(CALdouble));
    const uint1 sizeOfDouble1 = uint1(sizeof(CALdouble));

    uint1 inBounds = (rStep < rSteps) && (muStep < muSteps);

    il_if (inBounds)
    {
        /* 0 integrals and get stream constants */
        double1 bg_int = double1(0.0);
        std::vector<double1> streamIntegrals;

        emit_comment("Zero streams");
        for (j = 0; j < maxStreams; ++j)
            streamIntegrals.push_back(double1(0.0));


        const uint1 rPtsOffset = sizeOfDouble2 * nConvolve * rStep;
        uint1 i = rPtsOffset + rPtsBase;
        const uint1 rPtsMax = i + sizeOfDouble2 * nConvolve;

        uint1 trigIdx = nuStep * muSteps + muStep;

        uint1 lTrigIdx = sizeOfDouble2 * trigIdx + lTrigBase;
        uint1 bSinIdx = sizeOfDouble1 * trigIdx + bSinBase;

        emit_comment("Read Trig");
        const double2 lTrig = lTrigBuf[lTrigIdx];
        const double1 bSin = bTrigBuf[bSinIdx];

        il_whileloop
        {
            double2 rPt = rPts[i];

            emit_comment("coordinate conversion");
            double1 x = mad(rPt.x(), lTrig.x(), m_sun_r0);
            double1 y = rPt.x() * lTrig.y();
            double1 z = rPt.x() * bSin;

            emit_comment("x^2 + y^2 + z^2 / q^2");
            double1 tmp1 = x * x;
            double1 tmp3 = mad(y, y, tmp1);
            double1 tmp2 = z * z;
            double1 tmp4 = mad(q_inv_sqr, tmp2, tmp3);

            emit_comment("sqrt()");
            double1 rg = sqrt_custom(tmp4);
            double1 rs = rg + r0;

            emit_comment("qw_r3_N / (rg *rs^3)");
            bg_int += div_custom(rPt.y(), rg * rs * rs * rs);

            emit_comment("stream loops");
            for (j = 0; j < maxStreams; ++j)
            {
                emit_comment("Begin stream");

              #if FLEXIBLE_KERNEL
                il_if (nStream > uint1(j))
              #endif
                {
                    // a.x, c.x
                    named_variable<double2> streamX(str(format("cb4[%u]") % (4 * j + 0)));

                    // a.y, c.y
                    named_variable<double2> streamY(str(format("cb4[%u]") % (4 * j + 1)));

                    // a.z, c.z
                    named_variable<double2> streamZ(str(format("cb4[%u]") % (4 * j + 2)));

                    // sigma_sq2_inv, unused
                    named_variable<double1> sigma_sq2_inv(str(format("cb4[%u].xy") % (4 * j + 3)));

                    emit_comment("x - c");
                    double1 xs = x - streamX.y();
                    double1 ys = y - streamY.y();
                    double1 zs = z - streamZ.y();

                    emit_comment("a . x");
                    /* Dot product */
                    double1 dotted = streamX.x() * xs;
                    dotted = mad(streamY.x(), ys, dotted);
                    dotted = mad(streamZ.x(), zs, dotted);

                    emit_comment("xs += dotted * (-a)");
                    xs = mad(dotted, -streamX.x(), xs);
                    ys = mad(dotted, -streamY.x(), ys);
                    zs = mad(dotted, -streamZ.x(), zs);

                    emit_comment("sqrv");
                    double1 sqrv = xs * xs;
                    sqrv = mad(ys, ys, sqrv);
                    sqrv = mad(zs, zs, sqrv);

                    emit_comment("sqrv *= -sigma_sq2_inv");
                    sqrv *= -sigma_sq2_inv;

                    emit_comment("increment stream integral");
                    streamIntegrals[j] = mad(rPt.y(), exp_custom(sqrv), streamIntegrals[j]);
                }
              #if FLEXIBLE_KERNEL
                il_endif
              #endif
            }
            emit_comment("End streams");

            i += sizeOfDouble2;
            il_breakc(i >= rPtsMax);
        }
        il_endloop

        uint1 idx = muStep * rSteps + rStep;
        double1 V_reff_xr_rp3 = nuId * rConsts[sizeOfDouble2 * rStep + rConstsBase].x();

        emit_comment("read old values");

        SumType bgRead = bgOut[SIZEOF_SUM_TYPE * idx + bgOutBase];
        std::vector<SumType> streamRead;
        for (j = 0; j < maxStreams; ++j)
        {
            SumType val = SumType(0.0);

          #if FLEXIBLE_KERNEL
            il_if (nStream > uint1(j))
          #endif
            {
                val = streamsOut[SIZEOF_SUM_TYPE * (nStream * idx + uint1(j)) + streamsOutBase];
            }
          #if FLEXIBLE_KERNEL
            il_endif
          #endif

            streamRead.push_back(val);
        }

      #if USE_KAHAN
        emit_comment("multiply V_reff_xr_rp3 with Kahan summation");
        bg_int *= V_reff_xr_rp3
        bgOut = kahanAdd(bgRead, bg_int);
        for (j = 0; j < maxStreams; ++j)
        {
          #if FLEXIBLE_KERNEL
            il_if (nStream > uint1(j))
          #endif
            {
                streamIntegrals[j] *= V_reff_xr_rp3;
                streamIntegrajs[j] = kahanAdd(streamRead[j], streamIntegrals[j]);
            }
          #if FLEXIBLE_KERNEL
            il_endif
          #endif
        }
      #else
        emit_comment("multiply V_reff_xr_rp3");
        bg_int = mad(V_reff_xr_rp3, bg_int, bgRead);
        for (j = 0; j < maxStreams; ++j)
        {
            streamIntegrals[j] = mad(V_reff_xr_rp3, streamIntegrals[j], streamRead[j]);
        }
      #endif /* USE_KAHAN */

        emit_comment("Output");
        bgOut[SIZEOF_SUM_TYPE * idx + bgOutBase] = bg_int;

        for (j = 0; j < maxStreams; ++j)
        {
          #if FLEXIBLE_KERNEL
            il_if (nStream > uint1(j))
          #endif
            {
                streamsOut[SIZEOF_SUM_TYPE * (nStream * idx + uint1(j)) + streamsOutBase] = streamIntegrals[j];
            }
          #if FLEXIBLE_KERNEL
            il_endif
          #endif

        }
    }
    il_endif
}

static int separationKernelHeader(std::stringstream& code, CALuint maxStreams)
{
    if (4 * maxStreams > 1024)
    {
        std::cerr << format("Too many streams (%u)\n") % maxStreams << std::endl;
        return 1;
    }

    code << "il_cs_2_0\n";
  //code << format("dcl_max_thread_per_group %i\n") % 256;  // Only on Evergreen and later
    code << format("dcl_num_thread_per_group %i\n") % 64;  // R7XX and later
    code << format("; %i stream kernel\n") % maxStreams;
    code << "dcl_cb cb0[10]  ; Constant buffer that holds ABI data\n";
    code << "dcl_cb cb1[15]  ; Kernel arguments\n";
    code << "dcl_cb cb2[72]  ; I'm guessing the math constants AMD uses\n";
    code << "dcl_cb cb3[8]   ; ap constants\n";
    code << format("dcl_cb cb4[%i]  ; stream constants\n") % (4 * maxStreams);
    code << format("dcl_cb cb5[%i] ; sg_dx\n") % (IL_MAX_CONVOLVE / 2);

    return 0;
}

std::string createSeparationIntegralKernel(CALuint device, CALuint maxStreams)
{
    CALuint uavid = device >= CAL_TARGET_CYPRESS ? 11 : 1;
    std::stringstream code;

    if (separationKernelHeader(code, maxStreams))
        return "";

    Source::begin(device);

    createSeparationKernelCore(maxStreams, uavid);

    Source::end();

    Source::emitHeader(code);

    /* This must follow the UAV declaration, and CAL++ seems to be missing it. */
    code << "dcl_arena_uav_id(8)\n";

    Source::emitCode(code);

    return code.str();
}

static CALtargetEnum strToDeviceType(const char* str)
{
    if (!strcasecmp(str, "600"))
    {
        return CAL_TARGET_600;
    }
    else if (!strcasecmp(str, "610"))
    {
        return CAL_TARGET_610;
    }
    else if (!strcasecmp(str, "630"))
    {
        return CAL_TARGET_630;
    }
    else if (!strcasecmp(str, "670"))
    {
        return CAL_TARGET_670;
    }
    else if (!strcasecmp(str, "7XX"))
    {
        return CAL_TARGET_7XX;
    }
    else if (!strcasecmp(str, "770"))
    {
        return CAL_TARGET_770;
    }
    else if (!strcasecmp(str, "710"))
    {
        return CAL_TARGET_710;
    }
    else if (!strcasecmp(str, "730"))
    {
        return CAL_TARGET_730;
    }
    else if (!strcasecmp(str, "CYPRESS"))
    {
        return CAL_TARGET_CYPRESS;
    }
    else if (!strcasecmp(str, "JUNIPER"))
    {
        return CAL_TARGET_JUNIPER;
    }
    else if (!strcasecmp(str, "REDWOOD"))
    {
        return CAL_TARGET_REDWOOD;
    }
    else if (!strcasecmp(str, "CEDAR "))
    {
        return CAL_TARGET_CEDAR;
    }
    else if (!strcasecmp(str, "WRESTLER"))
    {
        return CAL_TARGET_WRESTLER;
    }
    else if (!strcasecmp(str, "CAYMAN"))
    {
        return CAL_TARGET_CAYMAN;
    }
    else if (!strcasecmp(str, "BARTS"))
    {
        return CAL_TARGET_BARTS;
    }
    else if (!strcasecmp(str, "RESERVED0"))
    {
        return CAL_TARGET_RESERVED0;
    }
    else if (!strcasecmp(str, "RESERVED1"))
    {
        return CAL_TARGET_RESERVED1;
    }
    else if (!strcasecmp(str, "RESERVED2"))
    {
        return CAL_TARGET_RESERVED2;
    }
    else
    {
        std::cerr << "Unknown device name " << str << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, const char* argv[])
{
    char* endp = NULL;

    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s <number of streams> <device type> <output file>\n", argc >= 1 ? argv[0] : "il_kernelgen");
        return 1;
    }

    CALuint nStream = (CALuint) strtoul(argv[1], &endp, 10);
    if (*endp != '\0')
    {
        perror("Reading number of streams");
        return 1;
    }

    CALuint deviceType = strToDeviceType(argv[2]);

    std::cout << "Creating a kernel with " << nStream << " for device type " << deviceType << std::endl;

    std::string src = createSeparationIntegralKernel(deviceType, nStream);

    std::ofstream outfile;
    outfile.open(argv[3], std:: ios::out);

    if (!outfile.is_open())
    {
        std::cerr << "Failed to open output file" << std::endl;
        return 1;
    }

    outfile << src;
    outfile.close();

    return 0;
}

