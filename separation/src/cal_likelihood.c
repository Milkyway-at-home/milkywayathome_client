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
    real* starsXY;
    real* starsZ;
    CALuint i;
    CALresult err = CAL_RESULT_OK;

    starsXY = mwMallocA(2 * sizeof(real) * sp->number_stars);
    starsZ = mwMallocA(sizeof(real) * sp->number_stars);

    for (i = 0; i < sp->number_stars; ++i)
    {
        starsXY[2 * i + 0] = X(sp->stars[i]);
        starsXY[2 * i + 1] = Y(sp->stars[i]);
        starsZ[i] = Z(sp->stars[i]);
    }

    err |= createConstantBuffer1D(&cm->starsXY, ci, starsXY, formatReal2, sp->number_stars);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create starsXY buffer", err);

    err |= createConstantBuffer1D(&cm->starsZ, ci, starsZ, formatReal1, sp->number_stars);
    if (err != CAL_RESULT_OK)
        cal_warn("Failed to create starsZ buffer", err);

    mwFreeA(starsXY);
    mwFreeA(starsZ);

    return err;
}

static CALresult createLikelihoodBuffers(const AstronomyParameters* ap,
                                         const StreamGauss sg,
                                         const StarPoints* sp,
                                         MWCALInfo* ci,
                                         SeparationCALMem* cm)
{
    CALresult err = CAL_RESULT_OK;

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

    memset(&cm, 0, sizeof(cm));
    if (separationLoadKernel(ci, ap, sc, CAL_TRUE) != CAL_RESULT_OK)
        fail("Failed to load likelihood kernel");





    err = releaseSeparationBuffers(ci, &cm);
    if (err != CAL_RESULT_OK)
        return err;

    mwUnloadKernel(ci);

    return CAL_RESULT_OK;
}


