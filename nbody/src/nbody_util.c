/*
 * Copyright (c) 2011-2012 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
 * Copyright (c) 2016-2018 Siddhartha Shelton
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_util.h"
#include "milkyway_math.h"

/* Correct timestep so an integer number of steps covers the exact
 * evolution time */
real_0 nbCorrectTimestep(real_0 time, real_0 dt)
{
    real_0 nStep = mw_ceil_0(time / dt);
    return time / nStep;
}

mwvector nbCenterOfMass(const NBodyState* st)
{
    int i;
    const Body* b;
    int nbody = st->nbody;
    mwvector cm = ZERO_VECTOR;
    mwvector tmp;
    Kahan mass;
    Kahan pos[3];

    CLEAR_KAHAN(&mass);
    CLEAR_KAHAN(&pos[0]);
    CLEAR_KAHAN(&pos[1]);
    CLEAR_KAHAN(&pos[2]);

    for (i = 0; i < nbody; ++i)
    {
        b = &st->bodytab[i];

        tmp.x = mw_mul(&Pos(b).x, &Mass(b));
        tmp.y = mw_mul(&Pos(b).y, &Mass(b));
        tmp.z = mw_mul(&Pos(b).z, &Mass(b));

        KAHAN_ADD(&pos[0], &tmp.x);
        KAHAN_ADD(&pos[1], &tmp.y);
        KAHAN_ADD(&pos[2], &tmp.z);
        KAHAN_ADD(&mass, &Mass(b));
    }

    X(&cm) = mw_div(&pos[0].sum, &mass.sum);
    Y(&cm) = mw_div(&pos[1].sum, &mass.sum);
    Z(&cm) = mw_div(&pos[2].sum, &mass.sum);
    W(&cm) = mass.sum;

    return cm;
}

mwvector nbCenterOfMass_Best(const NBodyState* st)
{
    int i;
    const Body* b;
    int nbody = st->nbody;
    mwvector cm = ZERO_VECTOR;
    mwvector tmp;
    Kahan mass;
    Kahan pos[3];

    CLEAR_KAHAN(&mass);
    CLEAR_KAHAN(&pos[0]);
    CLEAR_KAHAN(&pos[1]);
    CLEAR_KAHAN(&pos[2]);

    for (i = 0; i < nbody; ++i)
    {
        b = &st->bestLikelihoodBodyTab[i];

        tmp.x = mw_mul(&Pos(b).x, &Mass(b));
        tmp.y = mw_mul(&Pos(b).y, &Mass(b));
        tmp.z = mw_mul(&Pos(b).z, &Mass(b));

        KAHAN_ADD(&pos[0], &tmp.x);
        KAHAN_ADD(&pos[1], &tmp.y);
        KAHAN_ADD(&pos[2], &tmp.z);
        KAHAN_ADD(&mass, &Mass(b));
    }

    X(&cm) = mw_div(&pos[0].sum, &mass.sum);
    Y(&cm) = mw_div(&pos[1].sum, &mass.sum);
    Z(&cm) = mw_div(&pos[2].sum, &mass.sum);
    W(&cm) = mass.sum;

    return cm;
}

mwvector nbCenterOfMom(const NBodyState* st)
{
    int i;
    const Body* b;
    int nbody = st->nbody;
    mwvector cm = ZERO_VECTOR;
    mwvector tmp;
    Kahan mass;
    Kahan pos[3];

    CLEAR_KAHAN(&mass);
    CLEAR_KAHAN(&pos[0]);
    CLEAR_KAHAN(&pos[1]);
    CLEAR_KAHAN(&pos[2]);

    for (i = 0; i < nbody; ++i)
    {
        b = &st->bodytab[i];
        tmp.x = mw_mul(&Vel(b).x, &Mass(b));
        tmp.y = mw_mul(&Vel(b).y, &Mass(b));
        tmp.z = mw_mul(&Vel(b).z, &Mass(b));

        KAHAN_ADD(&pos[0], &tmp.x);
        KAHAN_ADD(&pos[1], &tmp.y);
        KAHAN_ADD(&pos[2], &tmp.z);
        KAHAN_ADD(&mass, &Mass(b));
    }

    X(&cm) = mw_div(&pos[0].sum, &mass.sum);
    Y(&cm) = mw_div(&pos[1].sum, &mass.sum);
    Z(&cm) = mw_div(&pos[2].sum, &mass.sum);
    W(&cm) = mass.sum;

    return cm;
}

mwvector nbCenterOfMom_Best(const NBodyState* st)
{
    int i;
    const Body* b;
    int nbody = st->nbody;
    mwvector cm = ZERO_VECTOR;
    mwvector tmp;
    Kahan mass;
    Kahan pos[3];

    CLEAR_KAHAN(&mass);
    CLEAR_KAHAN(&pos[0]);
    CLEAR_KAHAN(&pos[1]);
    CLEAR_KAHAN(&pos[2]);

    for (i = 0; i < nbody; ++i)
    {
        b = &st->bestLikelihoodBodyTab[i];
        tmp.x = mw_mul(&Vel(b).x, &Mass(b));
        tmp.y = mw_mul(&Vel(b).y, &Mass(b));
        tmp.z = mw_mul(&Vel(b).z, &Mass(b));

        KAHAN_ADD(&pos[0], &tmp.x);
        KAHAN_ADD(&pos[1], &tmp.y);
        KAHAN_ADD(&pos[2], &tmp.z);
        KAHAN_ADD(&mass, &Mass(b));
    }

    X(&cm) = mw_div(&pos[0].sum, &mass.sum);
    Y(&cm) = mw_div(&pos[1].sum, &mass.sum);
    Z(&cm) = mw_div(&pos[2].sum, &mass.sum);
    W(&cm) = mass.sum;

    return cm;
}

static inline real_0 log8(real_0 x)
{
    return mw_log_0(x) / mw_log_0(8.0);
}

/* The estimate formula has the unfortunate property of being negative
   for small n.  This will be the most negative. Add this as an extra
   boost to prevent negative flops estimates.
 */
static real_0 worstFlops(real_0 cQ, real_0 d, real_0 f)
{
    real_0 a = mw_pow_0(2.0, 3.0 - 3.0 * d / cQ);
    real_0 b = (cQ - d) * mw_log_0(8.0);
    real_0 c = cQ * mw_log_0(mw_pow_0(8.0, 1.0 - d / cQ));

    return -a * sqr_0(f) * (cQ + b - c) / (M_E * mw_log_0(8.0));
}

/* Estimate number of operations based on formula derived in
   "A Practical Comparison of N-Body Algorithms" (Blelloc, Narlikar 1995)

   Should be more accurate for more uniform distributions.  Does not
   include the flops from the external potential. However, the effect
   of the potential actually reduces the total number of flops by
   tearing apart the system in general.

   Does not account for newer opening criteria.
 */
real_0 nbEstimateNumberFlops(const NBodyCtx* ctx, int nbody)
{
    real_0 quadTerm, baseTerm;

    real_0 n = (real_0) nbody;
    real_0 nSteps = ctx->timeEvolve / ctx->timestep;

    /* Cost of interaction for a cell using a quadrupole moment. */
    const real_0 cQ = ctx->useQuad ? 50.0 : 0;

    /* Cost of a direct interaction */
    const real_0 d = 13;

    /* Based on BH86 opening criterion. */
    real_0 f = 28.0 * M_PI / (3.0 * cube_0(ctx->theta));

    /* FIXME: Don't be lazy and try rederiving for these. It should be
     * some number larger than for BH86. Somewhere I remember
     * something saying about 3x more operations for SW93
     */
    if (ctx->criterion != BH86)
    {
        f *= 3.0;
    }

    quadTerm = cQ * n * f * (log8(n / f) - 1.0);
    baseTerm = d * n * f;

    /* Total flops is then this times the number of timesteps */
    return nSteps * (quadTerm + baseTerm - worstFlops(cQ, d, f));
}

/* These estimates seem to sometimes work OK but very often not */
real_0 nbEstimateTime(const NBodyCtx* ctx, int nbody, real_0 flops)
{
    /* Spends < ~5% of the time in tree construction. Spends about
     * half the time in tree traversal / memory access as actually
     * calculating forces. */
    const real_0 factor = 2.05;

    /* Not 100% efficient. Bullshit number */
    const real_0 efficiency = 0.95;

    real_0 nflop = nbEstimateNumberFlops(ctx, nbody);

    return factor * nflop / (efficiency * flops);
}

void nbReportTreeIncest(const NBodyCtx* ctx, NBodyState* st)
{
    if (!st->treeIncest)   /* don't repeat warning */
    {
        st->treeIncest = TRUE;

        if (!ctx->quietErrors) /* Avoid massive printing of tests causing incest */
        {
            if (ctx->allowIncest)
            {
                mw_printf("[tree-incest detected at step %u / %u (%f%%)]\n",
                          st->step,
                          ctx->nStep,
                          100.0 * (real_0) st->step / (real_0) ctx->nStep
                    );
            }
            else
            {
                mw_printf("tree-incest detected (fatal) at step %u / %u (%f%%)\n",
                          st->step,
                          ctx->nStep,
                          100.0 * (real_0) st->step / (real_0) ctx->nStep
                    );
            }
        }
    }
}

