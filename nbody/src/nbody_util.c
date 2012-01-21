/*
Copyright (C) 2011  Matthew Arsenault
Copyright (c) 2011 Rensselaer Polytechnic Institute.

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

#include "nbody_util.h"
#include "milkyway_math.h"

mwvector nbCenterOfMass(const NBodyState* st)
{
    int i;
    const Body* b;
    mwvector cm = ZERO_VECTOR;
    mwvector tmp;
    real mass = 0.0;

    for (i = 0; i < st->nbody; ++i)
    {
        b = &st->bodytab[i];

        tmp = mw_mulvs(Pos(b), Mass(b));
        mass += Mass(b);
        mw_incaddv(cm, tmp);
    }

    mw_incdivs(cm, mass);

    return cm;
}

static inline double log8(double x)
{
    return log(x) / log(8.0);
}

/* The estimate formula has the unfortunate property of being negative
   for small n.  This will be the most negative. Add this as an extra
   boost to prevent negative flops estimates.
 */
static double worstFlops(double cQ, double d, double f)
{
    double a = pow(2.0, 3.0 - 3.0 * d / cQ);
    double b = (cQ - d) * log(8.0);
    double c = cQ * log(pow(8.0, 1.0 - d / cQ));

    return -a * sqr(f) * (cQ + b - c) / (M_E * log(8.0));
}

/* Estimate number of operations based on formula derived in
   "A Practical Comparison of N-Body Algorithms" (Blelloc, Narlikar 1995)

   Should be more accurate for more uniform distributions.  Does not
   include the flops from the external potential. However, the effect
   of the potential actually reduces the total number of flops by
   tearing apart the system in general.

   Does not account for newer opening criteria.
 */
double nbEstimateNumberFlops(const NBodyCtx* ctx, int nbody)
{
    double quadTerm, baseTerm;

    double n = (double) nbody;
    double nSteps = ctx->timeEvolve / ctx->timestep;

    /* Cost of interaction for a cell using a quadrupole moment. */
    const double cQ = ctx->useQuad ? 50.0 : 0;

    /* Cost of a direct interaction */
    const double d = 13;

    /* Based on BH86 opening criterion. */
    double f = 28.0 * M_PI / (3.0 * cube(ctx->theta));

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
double nbEstimateTime(const NBodyCtx* ctx, int nbody, double flops)
{
    /* Spends < ~5% of the time in tree construction. Spends about
     * half the time in tree traversal / memory access as actually
     * calculating forces. */
    const double factor = 2.05;

    /* Not 100% efficient. Bullshit number */
    const double efficiency = 0.95;

    double nflop = nbEstimateNumberFlops(ctx, nbody);

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
                          100.0 * (double) st->step / (double) ctx->nStep
                    );
            }
            else
            {
                mw_printf("tree-incest detected (fatal) at step %u / %u (%f%%)\n",
                          st->step,
                          ctx->nStep,
                          100.0 * (double) st->step / (double) ctx->nStep
                    );
            }
        }
    }
}

