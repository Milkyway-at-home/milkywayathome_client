/*
Copyright (C) 2011  Matthew Arsenault

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

#include "milkyway_util.h"
#include "separation_types.h"
#include "calculated_constants.h"
#include "r_points.h"
#include "show_cal_types.h"
#include "separation_cal_run.h"
#include "separation_cal_setup.h"
#include "separation_cal_types.h"
#include "separation_cal_kernelgen.h"
#include "cal_likelihood.h"
#include "integrals.h"

static CALresult createSGBuffers(const AstronomyParameters* ap, MWCALInfo* ci, SeparationCALMem* cm, StreamGauss sg)
{
    CALresult err;


    err = createConstantBuffer1D(&cm->sg_dx, ci, sg.dx, constantFormatReal1, ap->convolve);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create sg_dx buffer", err);

    err = createConstantBuffer1D(&cm->sg_qgauss_W, ci, sg.qgaus_W, constantFormatReal1, ap->convolve);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create sg_dx buffer", err);

    return CAL_RESULT_OK;

}


static CALresult createStarsBuffers(const StarPoints* sp, MWCALInfo* ci, SeparationCALMem* cm)
{
    real* starsLTrig;
    real* starsRBTrig;
    LB lb;
    LBTrig lbt;
    CALuint i;
    CALresult err = CAL_RESULT_OK;

    starsLTrig = mwMallocA(2 * sizeof(real) * sp->number_stars);
    starsRBTrig = mwMallocA(2 * sizeof(real) * sp->number_stars);

    for (i = 0; i < sp->number_stars; ++i)
    {
        LB_L(lb) = L(sp->stars[i]);
        LB_B(lb) = B(sp->stars[i]);
        lbt = lb_trig(lb);

        starsLTrig[2 * i + 0] = lbt.lCosBCos;
        starsLTrig[2 * i + 1] = lbt.lSinBCos;
        starsRBTrig[2 * i + 0] = R(sp->stars[i]);
        starsRBTrig[2 * i + 1] = lbt.bSin;
    }

    #if 0
    err |= createConstantBuffer1D(&cm->starsLTrig, ci, starsLTrig, formatReal2, sp->number_stars);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create starsLTrig buffer", err);

    err |= createConstantBuffer1D(&cm->starsRBTrig, ci, starsRBTrig, formatReal2, sp->number_stars);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create starsRBTrig buffer", err);
    #endif

    mwFreeA(starsLTrig);
    mwFreeA(starsRBTrig);

    return err;
}




static CALresult createLikelihoodBuffers(const AstronomyParameters* ap,
                                         const StreamGauss sg,
                                         const StarPoints* sp,
                                         MWCALInfo* ci,
                                         SeparationCALMem* cm)
{
    CALresult err = CAL_RESULT_OK;

    err |= createStarsBuffers(sp, ci, cm);
    err |= createSGBuffers(ap, ci, cm, sg);

    return err;
}

static CALresult getLikelihoodModuleNames(MWCALInfo* ci, SeparationCALNames* cn, CALuint numberStreams)
{
    CALresult err = CAL_RESULT_OK;


    if (err != CAL_RESULT_OK)
        cal_warn("Failed to get module names", err);

    return err;
}

static CALresult setLikelihoodKernelArguments(MWCALInfo* ci, SeparationCALMem* cm, SeparationCALNames* cn)
{
    CALresult err = CAL_RESULT_OK;


    if (err != CAL_RESULT_OK)
        cal_warn("Failed to set kernel arguments", err);

    return err;
}

CALresult likelihoodCAL(SeparationResults* results,
                        const AstronomyParameters* ap,
                        const StarPoints* sp,
                        const StreamConstants* sc,
                        const Streams* streams,
                        const StreamGauss sg,
                        const CLRequest* clr,
                        MWCALInfo* ci)
{
    CALresult err;
    SeparationCALMem cm;
    SeparationCALNames cn;

    memset(&cm, 0, sizeof(cm));
    if (separationLoadKernel(ci, ap, sc, CAL_TRUE) != CAL_RESULT_OK)
        fail("Failed to load likelihood kernel");



    err = createLikelihoodBuffers(ap, sg, sp, ci, &cm);
    if (err != CAL_RESULT_OK)
        return err;

    err = getLikelihoodModuleNames(ci, &cn, ap->number_streams);
    if (err != CAL_RESULT_OK)
        return err;

    err = setLikelihoodKernelArguments(ci, &cm, &cn);
    if (err != CAL_RESULT_OK)
        return err;

    /* Run kernel */

    err = releaseSeparationBuffers(ci, &cm);
    if (err != CAL_RESULT_OK)
        return err;

    mwUnloadKernel(ci);

    return CAL_RESULT_OK;
}


