/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#include <string.h>
#include <assert.h>
#include "nbody_config.h"
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "calc_params.h"

static int processHalo(Halo* h)
{
    if (h->type == TriaxialHalo)
    {
        const real phi = h->triaxAngle;
        const real cp  = mw_cos(phi);
        const real cps = sqr(cp);
        const real sp  = mw_sin(phi);
        const real sps = sqr(sp);

        const real qxs = sqr(h->flattenX);
        const real qys = sqr(h->flattenY);

        h->c1 = (cps / qxs) + (sps / qys);
        h->c2 = (cps / qys) + (sps / qxs);

        /* 2 * sin(x) * cos(x) == sin(2 * x) */
        h->c3 = mw_sin(2.0 * phi) * ((qys - qxs) / (qxs * qys));
    }

    return 0;
}

static int processPotential(Potential* p)
{
    return processHalo(&p->halo);
}

static int processModel(DwarfModel* mod)
{
    const real r0 = mod->scale_radius;
    real eps;

    if (isnan(mod->time_dwarf) && isnan(mod->time_orbit))
    {
        warn("At least one of the evolution times must be specified for the dwarf model\n");
        return 1;
    }

    switch (mod->type)
    {
        case DwarfModelPlummer:
            /* If not set, and no default, it's calculated based on
             * other parameters. */
            if (isnan(mod->eps2))
            {
                eps = r0 / (10.0 * mw_sqrt((real) mod->nbody));
                mod->eps2 = sqr(eps);
            }

            if (isnan(mod->timestep))
                mod->timestep = sqr(1/10.0) * mw_sqrt((PI_4_3 * cube(r0)) / mod->mass);

            /* for the orbit, use dt = dtnbody/2 to make sure we get enough orbit precision. */
            if (isnan(mod->orbit_timestep))
                mod->orbit_timestep = mod->timestep / 2.0;

            break;

        case DwarfModelKing:
        case DwarfModelDehnen:
        default:
            warn("Unhandled model type: %d\n", mod->type);
            return 1;
    }

    if (isnan(mod->time_orbit))
        mod->time_orbit = mod->timestep / 2.0;

    /* if (isnan(mod->time_dwarf))
           ?????;
    */
    return 0;
}

static int processInitialConditions(const NBodyCtx* ctx, InitialConditions* ic)
{
    if (!ic->useGalC)
    {
        /* We aren't given galactic coordinates, so convert them */
        if (ic->useRadians)
            ic->position = lbrToCartesian_rad(ctx, ic->position);
        else
            ic->position = lbrToCartesian(ctx, ic->position);
        ic->useGalC = TRUE;
    }

    return 0;
}

/* Calculate needed parameters from whatever we read in */
/* TODO: I'm dissatisfied with having to throw these checks here at
 * the end */
static int postProcess(NBodyCtx* ctx)
{
    int rc;

    rc = processPotential(&ctx->pot);
    rc |= processModel(&ctx->model);

    /* These other pieces are dependent on the others being set up
     * first */

    ctx->freqout = inv(ctx->model.timestep);

    if (ctx->model.nbody < 1)
    {
        warn("nbody = %d is absurd\n", ctx->model.nbody);
        rc |= 1;
    }

    return rc;
}

int setCtxConsts(NBodyCtx* ctx,
                 const FitParams* fitParams, /* Hacked in overrides for using server's args */
                 InitialConditions* ic,
                 const long setSeed)
{
    int rc = 0;

    /* Hack: Ignore these parameters in the file if using the command
     * line arguments. */
    if (fitParams->useFitParams)
    {
        ctx->model.mass         = fitParams->modelMass;
        ctx->model.scale_radius = fitParams->modelRadius;
        ctx->model.time_dwarf   = fitParams->simulationTime;
        ctx->model.time_orbit   = fitParams->reverseOrbitTime;
        ctx->seed               = setSeed;
    }

    rc |= postProcess(ctx);
    rc |= processInitialConditions(ctx, ic);

    if (rc)
        warn("Failed to set context constants\n");

    return rc;
}

