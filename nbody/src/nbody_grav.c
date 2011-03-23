/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
   Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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
#include "nbody_grav.h"

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */

typedef struct
{
    mwvector pos0;      /* point to evaluate field */
    mwvector acc0;      /* resulting acceleration */
    mwvector dr;        /* vector from q to pos0 */
    const Body* pskip;  /* skip in force evaluation */
    NBodyCell* qmem;         /* data shared with gravSub */
    real drsq;          /* squared distance to pos0 */
} ForceEvalState;

#define EMPTY_FORCE_EVAL_STATE(p) { Pos(p), ZERO_VECTOR, ZERO_VECTOR, p, NULL, 0.0 }


/* subDivP: decide if cell q is too close to accept as a single
 * term. Also sets qmem, dr, and drsq for use by gravsub.
 */
HOT
static inline mwbool subDivP(ForceEvalState* fest, NBodyCell* q)
{
    fest->dr = mw_subv(Pos(q), fest->pos0);   /* compute displacement */
    fest->drsq = mw_sqrv(fest->dr);           /* and find dist squared */
    fest->qmem = q;                           /* remember we know them */
    return (fest->drsq < Rcrit2(q));          /* apply standard rule */
}

HOT
static inline void cellQuadTerm(ForceEvalState* fest, const NBodyNode* q, const real drab)
{
    real dr5inv, drquaddr, phiquad;
    mwvector ai, quaddr;

    dr5inv = 1.0 / (sqr(fest->drsq) * drab); /* form dr^-5 */
    quaddr = mw_mulmv(Quad(q), fest->dr);    /* form Q * dr */
    drquaddr = mw_dotv(fest->dr, quaddr);    /* form dr * Q * dr */
    phiquad = -0.5 * dr5inv * drquaddr;      /* get quad. part of phi */
    phiquad = 5.0 * phiquad / fest->drsq;    /* save for acceleration */
    ai = mw_mulvs(fest->dr, phiquad);        /* components of acc. */
    mw_incsubv(fest->acc0, ai);              /* increment */
    mw_incmulvs(quaddr, dr5inv);
    mw_incsubv(fest->acc0, quaddr);          /* acceleration */
}

/* gravSub: compute contribution of node q to gravitational field at
 * point pos0, and add to running totals phi0 and acc0.
 */
HOT
static inline void gravSub(const NBodyCtx* ctx, ForceEvalState* fest, const NBodyNode* q)
{
    real drab, phii, mor3;
    mwvector ai;

    if (q != (NBodyNode*) fest->qmem)                    /* cant use memorized data? */
    {
        fest->dr = mw_subv(Pos(q), fest->pos0);       /* then compute sep. */
        fest->drsq = mw_sqrv(fest->dr);               /* and sep. squared */
    }

    fest->drsq += ctx->eps2;   /* use standard softening */
    drab = mw_sqrt(fest->drsq);
    phii = Mass(q) / drab;
    mor3 = phii / fest->drsq;
    ai = mw_mulvs(fest->dr, mor3);
    mw_incaddv(fest->acc0, ai);         /* ... and to total accel. */

    if (ctx->useQuad && isCell(q))      /* if cell, add quad term */
        cellQuadTerm(fest, q, drab);
}

/* treeScan: iterative routine to do force calculation, starting with
 * node q, which is typically the t.root cell. Watches for tree
 * incest.
 */
static inline mwbool treeScan(const NBodyCtx* ctx,
                              ForceEvalState* fest,
                              const NBodyNode* q)
{
    mwbool skipself = FALSE;

    while (q != NULL)               /* while not at end of scan */
    {
        if (   isCell(q)                       /* is node a cell and... */
            && subDivP(fest, (NBodyCell*) q))       /* too close to accept? */
        {
            q = More(q);            /* follow to next level */
        }
        else                    /* else accept this term */
        {
            if (q == (NBodyNode*) fest->pskip)    /* self-interaction? */
                skipself = TRUE;             /* then just skip it */
            else                             /* not self-interaction */
                gravSub(ctx, fest, q);       /* so compute gravity */
            q = Next(q);                     /* follow next link */
        }
    }

    return skipself;
}

static void reportTreeIncest(const NBodyCtx* ctx, NBodyState* st)
{
    if (!st->treeIncest)   /* don't repeat warning */
    {
        st->treeIncest = TRUE;

        if (!ctx->quietErrors) /* Avoid massive printing of tests causing incest */
        {
            if (ctx->allowIncest)
                warn("[hackGrav: tree-incest detected]\n");
            else
                warn("hackGrav: tree-incest detected (fatal)\n");
        }
    }
}

/* hackGrav: evaluate gravitational field on body p; checks to be
 * sure self-interaction was handled correctly if intree is true.
 */

/* Not inlined without inline from multiple calls in
 * mapForceBody(). Measurably better with the inline, but only
 * slightly. */
static inline mwvector hackGrav(const NBodyCtx* ctx, NBodyState* st, const NBodyNode* root, const Body* p)
{
    ForceEvalState fest = EMPTY_FORCE_EVAL_STATE(p);
    mwbool intree = Mass(p) > 0.0;

    /* scan tree from root, watch for tree incest */
    if (intree && !treeScan(ctx, &fest, root))
        reportTreeIncest(ctx, st);

    return fest.acc0;         /* and acceleration */
}

static inline void mapForceBody(const NBodyCtx* ctx, NBodyState* st)
{
    unsigned int i;
    const unsigned int nbody = st->nbody;  /* Prevent reload on each loop */
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
                a = hackGrav(ctx, st, (NBodyNode*) st->tree.root, b);

                externAcc = acceleration(&ctx->pot, Pos(b));
                st->acctab[i] = mw_addv(a, externAcc);
                break;

            case EXTERNAL_POTENTIAL_NONE:
                st->acctab[i] = hackGrav(ctx, st, (NBodyNode*) st->tree.root, &st->bodytab[i]);
                break;

            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                mw_panic("Implement me!\n");

            default:
                fail("Bad external potential type: %d\n", ctx->potentialType);
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

    rc = makeTree(ctx, st);
    if (nbodyStatusIsFatal(rc))
        return rc;

    mapForceBody(ctx, st);

    return incestStatusCheck(ctx, st); /* Check if incest occured during step */
}

