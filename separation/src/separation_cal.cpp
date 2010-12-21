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
#include <iostream>

#include "calculated_constants.h"
#include "r_points.h"
#include "separation_cal.h"
#include "separation_cal_kernelgen.h"


using namespace boost;
using namespace cal;
using namespace cal::il;

typedef struct
{
    Context      context;
    Program      program;
    Kernel       kernel;
    CommandQueue queue;
} SeparationCALInfo;

typedef struct
{
    Image2D outBg;
    Image2D outStreams;
    Image2D lbts;
    Image2D rPts;
    Image1D sg_dx;
    Image1D rc;
    NDRange range;
    Image2D         A0, A1, B0, B1, C;
} SeparationCALMem;


static void mapWriteBuffer(SeparationCALInfo& ci, Image& img, void* p, size_t size)
{
    void* mapP;
    CALuint pitch;

    mapP = ci.queue.mapMemObject(img, pitch);
    memcpy(mapP, p, size);
    ci.queue.unmapMemObject(img);
}

static void zeroBuffer(SeparationCALInfo& ci, Image& img, size_t size)
{
    void* p;
    CALuint pitch;

    p = ci.queue.mapMemObject(img, pitch);
    memset(p, 0, size);
    ci.queue.unmapMemObject(img);
}

/* Create buffers constant over r, and set kernel arguments */
static void createRBuffers(SeparationCALInfo& ci,
                           SeparationCALMem& cm,
                           const AstronomyParameters* ap,
                           const IntegralArea* ia,
                           const StreamGauss sg)
{

    RPoints* rPts;
    RConsts* rc;

    rPts = precalculateRPts(ap, ia, sg, &rc, TRUE);

    cm.rPts = Image2D(ci.context, ia->r_steps, ap->convolve, CAL_FORMAT_UNSIGNED_INT32_4, 0);
    mapWriteBuffer(ci, cm.rPts, (void*) rPts, ia->r_steps * ap->convolve * sizeof(RPoints));

    /* sg_dx, rc in constant buffers? */

    cm.rc = Image1D(ci.context, ia->r_steps, CAL_FORMAT_UNSIGNED_INT32_4, 0 );
    mapWriteBuffer(ci, cm.rc, (void*) rc, sizeof(RConsts) * ia->r_steps);

    cm.sg_dx = Image1D(ci.context, ia->r_steps, CAL_FORMAT_UNSIGNED_INT32_2, 0 );
    mapWriteBuffer(ci, cm.sg_dx, (void*) sg.dx, sizeof(real) * ia->r_steps);

    mwFreeA(rPts);
    mwFreeA(rc);
}


static void createLBTrigBuffer(SeparationCALInfo& ci,
                               SeparationCALMem& cm,
                               const AstronomyParameters* ap,
                               const IntegralArea* ia)
{
    LBTrig* lbts;
    CALuint width, height;

    width = ia->mu_steps;
    height = ia->nu_steps;

    lbts = precalculateLBTrig(ap, ia, TRUE);
    cm.lbts = Image2D(ci.context, width, 2 * height, CAL_FORMAT_UNSIGNED_INT32_4, 0);

    mapWriteBuffer(ci, cm.lbts, (void*) lbts, width * height * sizeof(LBTrig));

    mwFreeA(lbts);
}

static void createOutputBuffers(SeparationCALInfo& ci,
                                SeparationCALMem& cm,
                                const AstronomyParameters* ap,
                                const IntegralArea* ia)
{
    CALuint width, height;
    CALuint flags = 0;
    //CALuint flags = CAL_RESALLOC_GLOBAL_BUFFER;

    width = ia->mu_steps;
    height = ia->r_steps;

    /* Using CAL_RESALLOC_GLOBAL_BUFFER seems to have a minimum size or else it weirdly errors */
    cm.outBg = Image2D(ci.context, width, height, CAL_FORMAT_UNSIGNED_INT32_4, flags);
    zeroBuffer(ci, cm.outBg, width * height * sizeof(real));

    cm.outStreams = Image2D(ci.context, width, ap->number_streams * height, CAL_FORMAT_UNSIGNED_INT32_2, flags);
    zeroBuffer(ci, cm.outStreams, ap->number_streams * width * height * sizeof(real));
}


static void createProgram(SeparationCALInfo& ci,
                          const AstronomyParameters* ap,
                          const IntegralArea* ia,
                          const StreamConstants* sc)
{
    std::vector<Device> devices;
    std::string source;

    devices = ci.context.getInfo<CAL_CONTEXT_DEVICES>();
    source = createSeparationKernel(ap, ia, sc);

    std::cout << source << std::endl; // Uncomment to emit IL code
    ci.program = Program(ci.context, source.c_str(), source.length());
    ci.program.build(devices);

    ci.program.disassemble(std::cout);
}

static void createKernel(SeparationCALInfo& ci, const AstronomyParameters* ap, const IntegralArea* ia)
{
    ci.kernel = Kernel(ci.program, "main");

    /* bg out */
    ci.kernel.setArgBind(0, "o0");

    /* streams out */
    ci.kernel.setArgBind(1, "o1");

    /* ap */
    ci.kernel.setArgBind(2, "cb0");

    /* ia */
    ci.kernel.setArgBind(3, "cb1");

    /* stream constants */
    ci.kernel.setArgBind(4, "cb2");

    /* r consts */
    /* TODO: Move to image, not constant buffer */
    ci.kernel.setArgBind(5, "cb3");

    /* sg_dx */
    ci.kernel.setArgBind(6, "cb4");

    /* r_pts */
    ci.kernel.setArgBind(7, "i0");

    /* lb trig */
    ci.kernel.setArgBind(8, "i1");

    /* extra */
    ci.kernel.setArgBind(9, "cb6", 0, sizeof(CALuint));

    /* nu_id */
    ci.kernel.setArgBind(10, "cb7", 0, sizeof(real));


}

static CALint separationCreateCALBuffers(SeparationCALInfo& ci,
                                         SeparationCALMem& cm,
                                         const AstronomyParameters* ap,
                                         const IntegralArea* ia,
                                         const StreamGauss sg)
{
    try
    {
        createOutputBuffers(ci, cm, ap, ia);
        createLBTrigBuffer(ci, cm, ap, ia);
        createRBuffers(ci, cm, ap, ia, sg);
    }
    catch (Error err)
    {
        std::cerr << "Error creating buffers: " << err.what() << std::endl;
        return 1;
    }

    return 0;
}

static CALint separationSetKernelArgs(SeparationCALInfo& ci,
                                      SeparationCALMem& cm,
                                      const AstronomyParameters* ap,
                                      const IntegralArea* ia,
                                      const StreamGauss sg)
{
    const CALuint zero = 0;

    try
    {
        ci.kernel.setArg(0, cm.outBg);
        ci.kernel.setArg(1, cm.outStreams);
        //ci.kernel.setArg(2, ap);
        //ci.kernel.setArg(3, ia);
        //ci.kernel.setArg(4, sc);
        ci.kernel.setArg(5, cm.rc);
        ci.kernel.setArg(6, cm.sg_dx);
        ci.kernel.setArg(7, cm.rPts);
        ci.kernel.setArg(8, cm.lbts);
        ci.kernel.setArg(9, zero);
    } catch (Error err)
    {
        std::cerr << "Error setting kernel arguments: " << err.what() << std::endl;
        return 1;
    }

    return 0;
}

static void separationSetNuKernelArg(SeparationCALInfo& ci, CALuint nu_step)
{
    ci.kernel.setArg(10, nu_step);
}


static void createQueue(SeparationCALInfo& ci, SeparationCALMem& cm, int dev)
{
    std::vector<Device> devices = ci.context.getInfo<CAL_CONTEXT_DEVICES>();
    ci.queue = CommandQueue(ci.context, devices[dev]);
}

int separationCALInit(SeparationCALInfo& ci,
                      SeparationCALMem& cm,
                      const AstronomyParameters* ap,
                      const IntegralArea* ia,
                      const StreamConstants* sc,
                      const StreamGauss sg,
                      CALuint devNum)
{
    try
    {
        ci.context = Context(Device(devNum));
        createQueue(ci, cm, devNum);
    }
    catch (Error err)
    {
        std::cerr << format("Error getting context and queue for device %u: %s") % devNum % err.what()
                  << std::endl;
        return 1;
    }

    try
    {
        createProgram(ci, ap, ia, sc);
    }
    catch (Error err)
    {
        std::cerr << "Error creating program: " << err.what() << std::endl;
        return 1;
    }

    try
    {
        createKernel(ci, ap, ia);
    }
    catch (Error err)
    {
        std::cerr << "Error creating kernel: " << err.what() << std::endl;
        return 1;
    }

    if (separationCreateCALBuffers(ci, cm, ap, ia, sg))
        return 1;

    if (separationSetKernelArgs(ci, cm, ap, ia, sg))
        return 1;

    separationSetNuKernelArg(ci, 0);

    return 0;
}

void prepareRun(SeparationCALInfo& ci,
                SeparationCALMem& cm,
                const AstronomyParameters* ap,
                const IntegralArea* ia)
{
    cm.range = NDRange(ia->mu_steps, ia->r_steps);
}

static CALint runNuSteps(SeparationCALInfo& ci, SeparationCALMem& cm, const IntegralArea* ia)
{
    CALuint i;
    Event event;
    double t1, t2;

    try
    {
        for (i = 0; i < ia->nu_steps; ++i)
        {
            t1 = mwGetTime();

            ci.queue.enqueueNDRangeKernel(ci.kernel, cm.range, &event);
            ci.queue.flush();
            ci.queue.waitForEvent(event);

            t2 = mwGetTime();
            warn("Time[%u] = %f\n", i, t2 - t1);
        }

    }
    catch (Error err)
    {
        std::cerr << "Error running integral: " << err.what() << std::endl;
        return 1;
    }

    return 0;
}

void release_gpu_resources(SeparationCALInfo& ci)
{
    ci.queue   = CommandQueue();
    ci.kernel  = Kernel();
    ci.context = Context();
    ci.program = Program();
}

real integrateCAL_cpp(const AstronomyParameters* ap,
                      const IntegralArea* ia,
                      const StreamConstants* sc,
                      const StreamGauss sg,
                      real* st_probs,
                      EvaluationState* es,
                      const CLRequest* clr)
{
    SeparationCALInfo ci;
    SeparationCALMem cm;

    cal::Init();

    int devNum = 0;

    warn("Hello from C++ land\n");

    if (separationCALInit(ci, cm, ap, ia, sc, sg, devNum))
    {
        warn("Failed to init CAL run\n");
        return NAN;
    }

    warn("Setup complete\n");

    runNuSteps(ci, cm, ia);

    warn("Run complete\n");

    release_gpu_resources(ci);
    cal::Shutdown();

    return 0.0;
}

