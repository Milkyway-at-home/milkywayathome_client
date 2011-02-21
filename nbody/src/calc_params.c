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

static real calculateTimestep(real mass, real r0)
{
    return sqr(1/10.0) * mw_sqrt((PI_4_3 * cube(r0)) / mass);
}

static real calculateEps2(real nbody, real r0)
{
    real eps;

    eps = r0 / (10.0 * mw_sqrt(nbody));

    return sqr(eps);
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

static int processModel(NBodyCtx* ctx, DwarfModel* mod)
{
    int rc;

    rc = processInitialConditions(ctx, &mod->initialConditions);

    switch (mod->type)
    {
        case DwarfModelPlummer:
            if (isnan(mod->timestep))
                mod->timestep = calculateTimestep(mod->mass, mod->scale_radius);

            /* for the orbit, use dt = dtnbody/2 to make sure we get enough orbit precision. */
            if (isnan(mod->orbit_timestep))
                mod->orbit_timestep = mod->timestep / 2.0;

            break;

        case DwarfModelKing:
        case DwarfModelDehnen:
        default:
            warn("Unhandled model type: %s\n", showDwarfModelT(mod->type));
            return 1;
    }

    return rc;
}

static int processAllModels(NBodyCtx* ctx)
{
    real eps2;
    unsigned int i;
    int rc = 0;
    int findEps2 = isnan(ctx->eps2);
    int findTimestep = isnan(ctx->timestep);
    int findOrbitTimestep = isnan(ctx->orbit_timestep);

    /* Start out as high as possible */
    if (findTimestep)
        ctx->timestep = INFINITY;
    if (findOrbitTimestep)
        ctx->orbit_timestep = INFINITY;

    if (findEps2)      /* Setting this in the file overrides calculating it */
        ctx->eps2 = INFINITY;

    for (i = 0; i < ctx->modelNum; ++i)
    {
        rc |= processModel(ctx, &ctx->models[i]);

        /* Find the total number of bodies */
        ctx->nbody += ctx->models[i].nbody;

        /* Use a combination of light and dark models to generate timestep
        Assumes model 0 is light, model 1 is dark */

        /* Calculate the enclosed mass of the big sphere within the little sphere's scale length */
        double val, encmass = 0.0, step=1e-5;
        double smalla = ctx->models[0].scale_radius;
        double biga = ctx->models[1].scale_radius;
        double Ml = ctx->models[0].mass, Md = ctx->models[1].mass;
        double r;
	int intstep;

        if (findTimestep) {
           encmass = 0.0;
           for(r=0.0; r<=smalla; r+=step) {
                val = r*r / pow(r*r + biga*biga, 2.5);
                encmass += val*step;
           }
           encmass *= 3*Md*biga*biga;

           ctx->timestep = calculateTimestep(ctx->models[0].mass + encmass, ctx->models[0].scale_radius);
           warn("i = %i smalla = %f biga = %f Ml = %f Md = %f integral value = %f combined timestep = %f\n", i, smalla, biga, Ml, Md, encmass/(3*Md*biga*biga), ctx->timestep);
        }

        if (findOrbitTimestep) {
           encmass = 0.0;
           for(r=0.0; r<=smalla; r+=step) {
                val = r*r / pow(r*r + biga*biga, 2.5);
                encmass += val*step;
           }
           encmass *= 3*Md*biga*biga;

           /* Use 1/10th of the nbody timestep */ 
           ctx->orbit_timestep = calculateTimestep(ctx->models[0].mass + encmass, ctx->models[0].scale_radius)/10.0;
           warn("i = %i smalla = %f biga = %f Ml = %f Md = %f integral value = %f combined orbit timestep = %f\n", i, smalla, biga, Ml, Md, encmass/(3*Md*biga*biga), ctx->orbit_timestep);
        }

        if (findEps2)
        {
            /* Find the smallest eps2 */
            eps2 = calculateEps2((real) ctx->models[i].nbody, ctx->models[i].scale_radius);
            ctx->eps2 = mw_fmin(ctx->eps2, eps2);
        }

    }

    return rc;
}

/* Calculate needed parameters from whatever we read in */
static int postProcess(NBodyCtx* ctx)
{
    int rc = 0;

    rc |= processPotential(&ctx->pot);
    rc |= processAllModels(ctx);

    /* These other pieces are dependent on the others being set up
     * first */

    ctx->freqOut = 4;

    if (isnan(ctx->time_orbit))
        ctx->time_orbit = ctx->timestep / 2.0;

    /* if (isnan(ctx->time_dwarf))
           ?????;
    */

    return rc;
}

/* Make sure we aren't ignoring all of the models */
static int hasNonIgnorableModel(const NBodyCtx* ctx)
{
    unsigned int i;
    unsigned int totalNonIgnorableModels = 0;
    int rc;

    if (ctx->modelNum == 0)
    {
        warn("Context has no models\n");
        return 1;
    }

    if (ctx->modelNum > ((1 << 8 * sizeof(body_t)) - 1))
    {
        warn("Too many models: %u\n", ctx->modelNum);
        return 1;
    }

    for (i = 0; i < ctx->modelNum; ++i)
    {
        if (!ctx->models[i].ignoreFinal)
            ++totalNonIgnorableModels;
    }

    rc = (totalNonIgnorableModels == 0);
    if (rc)
        warn("Trying to ignore all of %u models\n", ctx->modelNum);

    return rc;
}

static int hasAcceptableEps2(const NBodyCtx* ctx)
{
    int rc = !isfinite(ctx->eps2) || ctx->eps2 <= 0.0;
    if (rc)
        warn("Got an absurd eps2\n");

    return rc;
}

static int hasAcceptableTimes(const NBodyCtx* ctx)
{
    int rc;

    rc =   !isnormal(ctx->time_evolve)
        || !isnormal(ctx->time_orbit)
        || ctx->time_evolve <= 0.0
        || ctx->time_orbit <= 0.0;

    if (rc)
        warn("Got an unacceptable orbit or evolution time\n");
    return rc;
}

static int hasAcceptableSteps(const NBodyCtx* ctx)
{
    int rc =    !isnormal(ctx->timestep)
             || !isnormal(ctx->orbit_timestep)
             || (ctx->timestep <= 0.0)
             || (ctx->orbit_timestep <= 0.0);
    if (rc)
        warn("Context has unacceptable timesteps\n");

    return rc;
}

static int hasAcceptableNbody(const NBodyCtx* ctx)
{
    int rc = ctx->nbody < 1;
    if (rc)
        warn("nbody = %d is absurd\n", ctx->nbody);
    return rc;
}

static int contextSanityCheck(const NBodyCtx* ctx)
{
    int rc = 0;

    rc |= hasNonIgnorableModel(ctx);
    rc |= hasAcceptableNbody(ctx);
    rc |= hasAcceptableTimes(ctx);
    rc |= hasAcceptableSteps(ctx);
    rc |= hasAcceptableEps2(ctx);

    return rc;
}

int setCtxConsts(NBodyCtx* ctx,
                 const FitParams* fitParams, /* Hacked in overrides for using server's args */
                 const long setSeed)
{
    int rc = 0;

    /* Hack: Ignore these parameters in the file if using the command
     * line arguments. */
    /* FIXME: Assumes using first model */
    /* willeb: assumes the first model is the smaller, light model */
    if (fitParams->useFitParams)
    {
        ctx->models[0].mass         = fitParams->modelMass;
        ctx->models[0].scale_radius = fitParams->modelRadius;

        ctx->time_orbit  = fitParams->reverseOrbitTime;
        ctx->time_evolve = fitParams->simulationTime;
        ctx->seed        = setSeed;
    }

    rc |= postProcess(ctx);
    rc |= contextSanityCheck(ctx);
    if (rc)
        warn("Failed to set context constants\n");

    return rc;
}

