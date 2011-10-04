/*
 *  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_priv.h"
#include "milkyway_util.h"
#include "nbody_grav.h"

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */


static void reportTreeIncest(const NBodyCtx* ctx, NBodyState* st)
{
    if (!st->treeIncest)   /* don't repeat warning */
    {
        st->treeIncest = TRUE;

        if (!ctx->quietErrors) /* Avoid massive printing of tests causing incest */
        {
            if (ctx->allowIncest)
            {
                mw_printf("[hackGrav: tree-incest detected at time %f / %f (%f%%)]\n",
                          st->tnow,
                          ctx->timeEvolve,
                          100.0 * st->tnow / ctx->timeEvolve
                    );
            }
            else
            {
                mw_printf("hackGrav: tree-incest detected (fatal) at time %f / %f (%f%%)\n",
                          st->tnow,
                          ctx->timeEvolve,
                          100.0 * st->tnow / ctx->timeEvolve
                    );
            }
        }
    }
}


/*
 * nbodyGravity: Walk the tree starting at the root to do force
 * calculations.
 *
 *
 * Random notes:
 *   - Not inlined without inline from multiple calls in
 *     mapForceBody(). Measurably better with the inline, but only
 *     slightly.
 */
static inline mwvector nbodyGravity(const NBodyCtx* ctx, NBodyState* st, const Body* p)
{
    mwbool skipSelf = FALSE;

    mwvector pos0 = Pos(p);
    mwvector acc0 = ZERO_VECTOR;

    const NBodyNode* q = (const NBodyNode*) st->tree.root; /* Start at the root */

    while (q != NULL)               /* while not at end of scan */
    {
        mwvector dr = mw_subv(Pos(q), pos0);   /* Then compute distance */
        real drsq = mw_sqrv(dr);               /* and distance squared */

        if (isBody(q) || (drsq >= Rcrit2(q)))      /* If is a body or far enough away to approximate */
        {
            if (mw_likely((const Body*) q != p))   /* self-interaction? */
            {
                real drab, phii, mor3;

                /* Compute gravity */

                drsq += ctx->eps2;   /* use standard softening */
                drab = mw_sqrt(drsq);
                phii = Mass(q) / drab;
                mor3 = phii / drsq;

                mw_incaddv(acc0, mw_mulvs(dr, mor3));   /* ... and to total accel. */

                if (ctx->useQuad && isCell(q))          /* if cell, add quad term */
                {
                    real dr5inv, drquaddr, phiquad;
                    mwvector quaddr;
                    mwvector ai;

                    dr5inv = 1.0 / (sqr(drsq) * drab);  /* form dr^-5 */
                    quaddr = mw_mulmv(Quad(q), dr);     /* form Q * dr */
                    drquaddr = mw_dotv(dr, quaddr);     /* form dr * Q * dr */
                    phiquad = -0.5 * dr5inv * drquaddr; /* get quad. part of phi */
                    phiquad = 5.0 * phiquad / drsq;     /* save for acceleration */
                    ai = mw_mulvs(dr, phiquad);         /* components of acc. */
                    mw_incsubv(acc0, ai);               /* increment */
                    mw_incmulvs(quaddr, dr5inv);
                    mw_incsubv(acc0, quaddr);           /* acceleration */
                }
            }
            else
            {
                skipSelf = TRUE;   /* Encountered self */
            }

            q = Next(q);  /* Follow next link */
        }
        else
        {
             q = More(q); /* Follow to the next level if need to go deeper */
        }
    }

    if (!skipSelf)
    {
        /* If a body does not encounter itself in its traversal of the
         * tree, it is "tree incest" */

        reportTreeIncest(ctx, st);
    }

    return acc0;
}

static inline void mapForceBody(const NBodyCtx* ctx, NBodyState* st)
{
    int i;
    const int nbody = st->nbody;  /* Prevent reload on each loop */
    mwvector a, externAcc;
    const Body* b;

  #ifdef _OPENMP
    #pragma omp parallel for private(i, b, a, externAcc) schedule(dynamic)
  #endif
    for (i = 0; i < nbody; ++i)      /* get force on each body */
    {
        /* Repeat the base hackGrav part in each case or else GCC's
         * -funswitch-loops doesn't happen. Without that this constant
         * gets checked on every body on every step which is dumb.  */
        switch (ctx->potentialType)
        {
            case EXTERNAL_POTENTIAL_DEFAULT:
                /* Include the external potential */
                b = &st->bodytab[i];
                a = nbodyGravity(ctx, st, b);

                externAcc = acceleration(&ctx->pot, Pos(b));
                st->acctab[i] = mw_addv(a, externAcc);
                break;

            case EXTERNAL_POTENTIAL_NONE:
                st->acctab[i] = nbodyGravity(ctx, st, &st->bodytab[i]);
                break;

            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                mw_panic("Implement me!\n");

            default:
                mw_fail("Bad external potential type: %d\n", ctx->potentialType);
        }
    }
}

static inline mwvector nbodyGravity_Exact(const NBodyCtx* ctx, NBodyState* st, const Body* p)
{
    int i;
    const int nbody = st->nbody;
    mwvector a = ZERO_VECTOR;
    const real eps2 = ctx->eps2;

    for (i = 0; i < nbody; ++i)
    {
        const Body* b = &st->bodytab[i];

        mwvector dr = mw_subv(Pos(b), Pos(p));
        real drsq = mw_sqrv(dr) + eps2;

        real drab = mw_sqrt(drsq);
        real phii = Mass(b) / drab;
        real mor3 = phii / drsq;

        mw_incaddv(a, mw_mulvs(dr, mor3));
    }

    return a;
}

static inline void mapForceBody_Exact(const NBodyCtx* ctx, NBodyState* st)
{
    int i;
    const int nbody = st->nbody;  /* Prevent reload on each loop */
    mwvector a, externAcc;
    const Body* b;

  #ifdef _OPENMP
    #pragma omp parallel for private(i, b, a, externAcc) schedule(dynamic)
  #endif
    for (i = 0; i < nbody; ++i)      /* get force on each body */
    {
        switch (ctx->potentialType)
        {
            case EXTERNAL_POTENTIAL_DEFAULT:
                b = &st->bodytab[i];
                a = nbodyGravity_Exact(ctx, st, b);

                externAcc = acceleration(&ctx->pot, Pos(b));
                st->acctab[i] = mw_addv(a, externAcc);
                break;

            case EXTERNAL_POTENTIAL_NONE:
                st->acctab[i] = nbodyGravity_Exact(ctx, st, &st->bodytab[i]);
                break;

            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                mw_panic("Implement me!\n");

            default:
                mw_fail("Bad external potential type: %d\n", ctx->potentialType);
        }
    }

}

static inline NBodyStatus incestStatusCheck(const NBodyCtx* ctx, const NBodyState* st)
{
    if (st->treeIncest)
        return ctx->allowIncest ? NBODY_TREE_INCEST_NONFATAL : NBODY_TREE_INCEST_FATAL;

    return NBODY_SUCCESS;
}

NBodyStatus gravMap(const NBodyCtx* ctx, NBodyState* st)
{
    NBodyStatus rc;

    if (mw_likely(ctx->criterion != Exact))
    {
        rc = makeTree(ctx, st);
        if (nbodyStatusIsFatal(rc))
            return rc;

        mapForceBody(ctx, st);
    }
    else
    {
        mapForceBody_Exact(ctx, st);
    }

    return incestStatusCheck(ctx, st); /* Check if incest occured during step */
}

