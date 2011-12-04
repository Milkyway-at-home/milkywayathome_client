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
#include "nbody_util.h"
#include "nbody_grav.h"

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */


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
static inline mwvector nbGravity(const NBodyCtx* ctx, NBodyState* st, const Body* p)
{
    mwbool skipSelf = FALSE;

    mwvector pos0 = Pos(p);
    mwvector acc0 = ZERO_VECTOR;

    const NBodyNode* q = (const NBodyNode*) st->tree.root; /* Start at the root */

    while (q != NULL)               /* while not at end of scan */
    {
        mwvector dr = mw_subv(Pos(q), pos0);   /* Then compute distance */
        real drSq = mw_sqrv(dr);               /* and distance squared */

        if (isBody(q) || (drSq >= Rcrit2(q)))      /* If is a body or far enough away to approximate */
        {
            if (mw_likely((const Body*) q != p))   /* self-interaction? */
            {
                real drab, phii, mor3;

                /* Compute gravity */

                drSq += ctx->eps2;   /* use standard softening */
                drab = mw_sqrt(drSq);
                phii = Mass(q) / drab;
                mor3 = phii / drSq;

                acc0.x += mor3 * dr.x;
                acc0.y += mor3 * dr.y;
                acc0.z += mor3 * dr.z;

                if (ctx->useQuad && isCell(q))          /* if cell, add quad term */
                {
                    real dr5inv, drQdr, phiQ;
                    mwvector Qdr;

                    /* form Q * dr */
                    Qdr.x = Quad(q).xx * dr.x + Quad(q).xy * dr.y + Quad(q).xz * dr.z;
                    Qdr.y = Quad(q).xy * dr.x + Quad(q).yy * dr.y + Quad(q).yz * dr.z;
                    Qdr.z = Quad(q).xz * dr.x + Quad(q).yz * dr.y + Quad(q).zz * dr.z;


                    /* form dr * Q * dr */
                    drQdr = Qdr.x * dr.x + Qdr.y * dr.y + Qdr.z * dr.z;

                    dr5inv = 1.0 / (sqr(drSq) * drab);  /* form dr^-5 */

                    /* get quad. part of phi */
                    phiQ = 2.5 * (dr5inv * drQdr) / drSq;

                    acc0.x += phiQ * dr.x;
                    acc0.y += phiQ * dr.y;
                    acc0.z += phiQ * dr.z;

                    /* acceleration */
                    acc0.x -= dr5inv * Qdr.x;
                    acc0.y -= dr5inv * Qdr.y;
                    acc0.z -= dr5inv * Qdr.z;
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

        nbReportTreeIncest(ctx, st);
    }

    return acc0;
}

static inline void nbMapForceBody(const NBodyCtx* ctx, NBodyState* st)
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
                a = nbGravity(ctx, st, b);

                externAcc = nbExtAcceleration(&ctx->pot, Pos(b));
                mw_incaddv(a, externAcc);
                st->acctab[i] = a;
                break;

            case EXTERNAL_POTENTIAL_NONE:
                st->acctab[i] = nbGravity(ctx, st, &st->bodytab[i]);
                break;

            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                a = nbGravity(ctx, st, &st->bodytab[i]);
                nbEvalPotentialClosure(st, Pos(&st->bodytab[i]), &externAcc);
                mw_incaddv(a, externAcc)
                st->acctab[i] = a;
                break;

            default:
                mw_fail("Bad external potential type: %d\n", ctx->potentialType);
        }
    }
}

static mwvector nbGravity_Exact(const NBodyCtx* ctx, NBodyState* st, const Body* p)
{
    int i;
    const int nbody = st->nbody;
    mwvector a = ZERO_VECTOR;
    const real eps2 = ctx->eps2;

    for (i = 0; i < nbody; ++i)
    {
        const Body* b = &st->bodytab[i];

        mwvector dr = mw_subv(Pos(b), Pos(p));
        real drSq = mw_sqrv(dr) + eps2;

        real drab = mw_sqrt(drSq);
        real phii = Mass(b) / drab;
        real mor3 = phii / drSq;

        mw_incaddv(a, mw_mulvs(dr, mor3));
    }

    return a;
}

static inline void nbMapForceBody_Exact(const NBodyCtx* ctx, NBodyState* st)
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
                a = nbGravity_Exact(ctx, st, b);
                mw_incaddv(a, nbExtAcceleration(&ctx->pot, Pos(b)));
                st->acctab[i] = a;
                break;

            case EXTERNAL_POTENTIAL_NONE:
                st->acctab[i] = nbGravity_Exact(ctx, st, &st->bodytab[i]);
                break;

            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                a = nbGravity_Exact(ctx, st, &st->bodytab[i]);
                nbEvalPotentialClosure(st, Pos(&st->bodytab[i]), &externAcc);
                mw_incaddv(a, externAcc);
                st->acctab[i] = a;
                break;

            default:
                mw_fail("Bad external potential type: %d\n", ctx->potentialType);
        }
    }
}

static inline NBodyStatus nbIncestStatusCheck(const NBodyCtx* ctx, const NBodyState* st)
{
    if (st->treeIncest)
    {
        return ctx->allowIncest ? NBODY_TREE_INCEST_NONFATAL : NBODY_TREE_INCEST_FATAL;
    }

    return NBODY_SUCCESS;
}

NBodyStatus nbGravMap(const NBodyCtx* ctx, NBodyState* st)
{
    NBodyStatus rc;

    if (mw_likely(ctx->criterion != Exact))
    {
        rc = makeTree(ctx, st);
        if (nbStatusIsFatal(rc))
            return rc;

        nbMapForceBody(ctx, st);
    }
    else
    {
        nbMapForceBody_Exact(ctx, st);
    }

    if (st->potentialEvalError)
    {
        return NBODY_LUA_POTENTIAL_ERROR;
    }

    return nbIncestStatusCheck(ctx, st); /* Check if incest occured during step */
}

