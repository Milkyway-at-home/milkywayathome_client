/* Copyright 2010, 2011 Matthew Arsenault, Travis Desell, Dave Przybylo,
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
#include "nbody_check_params.h"

mwbool checkSphericalConstants(Spherical* s)
{
    if (mwCheckNormalPosNum(s->mass) || mwCheckNormalPosNum(s->scale))
    {
        warn("Invalid parameters for '%s': mass = %.15f, scale = %.15f\n",
             showSphericalT(s->type), s->mass, s->scale);
        return TRUE;
    }

    return FALSE;
}

mwbool checkDiskConstants(Disk* d)
{
    mwbool badDisk = FALSE;

    if (mwCheckNormalPosNum(d->mass))
    {
        warn("Invalid disk mass (%.15f)\n", d->mass);
        badDisk = TRUE;
    }

    switch (d->type)
    {
        case MiyamotoNagaiDisk:
            if (mwCheckNormalPosNum(d->scaleLength) || mwCheckNormalPosNum(d->scaleHeight))
            {
                warn("Invalid parameters for disk type '%s': scaleLength = %.15f, scaleHeight = %.15f\n",
                     showDiskT(d->type),
                     d->scaleLength,
                     d->scaleHeight);
                badDisk = TRUE;
            }
            break;

        case ExponentialDisk:
            if (isnormal(d->scaleHeight))
            {
                warn("Scale height unused for disk type '%s'\n", showDiskT(d->type));
                badDisk = TRUE;
            }

            if (mwCheckNormalPosNum(d->scaleLength))
            {
                warn("Invalid parameter for disk type '%s': scaleLength = %.15f\n",
                     showDiskT(d->type),
                     d->scaleLength);
                badDisk = TRUE;
            }

            break;

        case InvalidDisk:
        default:
            return warn1("Invalid disk type: %s (%d)\n", showDiskT(d->type), d->type);
    }

    return badDisk;
}

static mwbool invalidHaloWarning(halo_t type)
{
    warn("Got non-finite required field for halo type '%s'\n", showHaloT(type));
    return TRUE;
}

/* Check for valid halo values and calculate constants. Return true if error. */
mwbool checkHaloConstants(Halo* h)
{
    real phi, cp, cps, sp, sps;
    real qxs, qys;

    /* Common to all 3 models */
    if (!isfinite(h->vhalo) || mwCheckNormalPosNum(h->scaleLength))
        return invalidHaloWarning(h->type);

    /* TODO: Check for set things which don't match the model? */
    switch (h->type)
    {
        case LogarithmicHalo:
            if (!isfinite(h->flattenZ))
                return invalidHaloWarning(h->type);
            break;

        case NFWHalo:
            break;

        case TriaxialHalo:
            if (   !isfinite(h->triaxAngle)
                || !isfinite(h->flattenX)
                || !isfinite(h->flattenY)
                || !isfinite(h->flattenZ))
            {
                return invalidHaloWarning(h->type);
            }

            phi = d2r(h->triaxAngle);
            cp  = mw_cos(phi);
            cps = sqr(cp);
            sp  = mw_sin(phi);
            sps = sqr(sp);

            qxs = sqr(h->flattenX);
            qys = sqr(h->flattenY);

            h->c1 = (cps / qxs) + (sps / qys);
            h->c2 = (cps / qys) + (sps / qxs);

            /* 2 * sin(x) * cos(x) == sin(2 * x) */
            h->c3 = mw_sin(2.0 * phi) * ((qys - qxs) / (qxs * qys));
            break;

        case InvalidHalo:
        default:
            return warn1("Trying to use invalid halo type: %s (%d)\n", showHaloT(h->type), h->type);
    }

    return FALSE;
}

mwbool checkPotentialConstants(Potential* p)
{
    return checkSphericalConstants(&p->sphere[0]) || checkDiskConstants(&p->disk) || checkHaloConstants(&p->halo);
}

static int hasAcceptableEps2(const NBodyCtx* ctx)
{
    int rc = mwCheckNormalPosNumEps(ctx->eps2);
    if (rc)
        warn("Got an absurd eps2 (%.15f)\n", ctx->eps2);

    return rc;
}

static int hasAcceptableTimes(const NBodyCtx* ctx)
{
    int rc = mwCheckNormalPosNumEps(ctx->timeEvolve);
    if (rc)
        warn("Got an unacceptable evolution time (%.15f)\n", ctx->timeEvolve);
    return rc;
}

static int hasAcceptableSteps(const NBodyCtx* ctx)
{
    int rc = mwCheckNormalPosNumEps(ctx->timestep);
    if (rc)
        warn("Got an unacceptable timestep (%.15f)\n", ctx->timestep);

    return rc;
}

mwbool checkNBodyCtxConstants(const NBodyCtx* ctx)
{
    return hasAcceptableTimes(ctx) || hasAcceptableSteps(ctx) || hasAcceptableEps2(ctx);
}

