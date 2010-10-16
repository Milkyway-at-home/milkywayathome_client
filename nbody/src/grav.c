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
#include "grav.h"

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */

typedef struct
{
    mwvector pos0;    /* point to evaluate field */
    mwvector acc0;    /* resulting acceleration */
    bodyptr pskip;    /* skip in force evaluation */
    cellptr qmem;     /* data shared with gravsub */
    mwvector dr;      /* vector from q to pos0 */
    real drsq;        /* squared distance to pos0 */
} ForceEvalState;

#define EMPTY_FORCE_EVAL_STATE { ZERO_VECTOR, ZERO_VECTOR, NULL, NULL, ZERO_VECTOR, 0.0 }


/* subdivp: decide if cell q is too close to accept as a single
 * term. Also sets qmem, dr, and drsq for use by gravsub.
 */
static inline bool subdivp(ForceEvalState* fest, cellptr q)
{
    fest->dr = mw_subv(Pos(q), fest->pos0);   /* compute displacement */
    fest->drsq = mw_sqrv(fest->dr);           /* and find dist squared */
    fest->qmem = q;                           /* remember we know them */
    return (fest->drsq < Rcrit2(q));          /* apply standard rule */
}


static inline void cellQuadTerm(ForceEvalState* fest, const nodeptr q, const real drab)
{
    real dr5inv, drquaddr, phiquad;
    mwvector ai, quaddr;

    dr5inv = 1.0 / (sqr(fest->drsq) * drab); /* form dr^-5 */
    quaddr = mw_mulmv(Quad(q), fest->dr);    /* form Q * dr */
    drquaddr = mw_dotv(fest->dr, quaddr);    /* form dr * Q * dr */
    phiquad = -0.5 * dr5inv * drquaddr;      /* get quad. part of phi */
    phiquad = 5.0 * phiquad / fest->drsq;    /* save for acceleration */
    ai = mw_mulvs(phiquad, fest->dr);        /* components of acc. */
    mw_incsubv(fest->acc0, ai);              /* increment */
    mw_incmulvs(quaddr, dr5inv);
    mw_incsubv(fest->acc0, quaddr);          /* acceleration */
}

/* gravsub: compute contribution of node q to gravitational field at
 * point pos0, and add to running totals phi0 and acc0.
 */
static inline void gravsub(const NBodyCtx* ctx, ForceEvalState* fest, const nodeptr q)
{
    real drab, phii, mor3;
    mwvector ai;

    if (q != (nodeptr) fest->qmem)                    /* cant use memorized data? */
    {
        fest->dr = mw_subv(Pos(q), fest->pos0);       /* then compute sep. */
        fest->drsq = mw_sqrv(fest->dr);               /* and sep. squared */
    }

    fest->drsq += ctx->model.eps2;   /* use standard softening */
    drab = mw_sqrt(fest->drsq);
    phii = Mass(q) / drab;
    mor3 = phii / fest->drsq;
    ai = mw_mulvs(mor3, fest->dr);
    mw_incaddv(fest->acc0, ai);         /* ... and to total accel. */

    if (ctx->usequad && Type(q) == CELL)             /* if cell, add quad term */
        cellQuadTerm(fest, q, drab);
}

/* treescan: iterative routine to do force calculation, starting with
 * node q, which is typically the t.root cell. Watches for tree
 * incest.
 */
static inline bool treescan(const NBodyCtx* ctx,
                            ForceEvalState* fest,
                            nodeptr q)
{
    bool skipself = FALSE;

    while (q != NULL)               /* while not at end of scan */
    {
        if (   Type(q) == CELL                   /* is node a cell and... */
            && subdivp(fest, (cellptr) q))       /* too close to accept? */
        {
            q = More(q);            /* follow to next level */
        }
        else                    /* else accept this term */
        {
            if (q == (nodeptr) fest->pskip)    /* self-interaction? */
                skipself = TRUE;               /* then just skip it */
            else                               /* not self-interaction */
                gravsub(ctx, fest, q);         /* so compute gravity */
            q = Next(q);            /* follow next link */
        }
    }

    return skipself;
}

/* hackGrav: evaluate gravitational field on body p; checks to be
 * sure self-interaction was handled correctly if intree is true.
 */
static inline mwvector hackGrav(const NBodyCtx* ctx, nodeptr root, bodyptr p)
{
    static bool treeincest = FALSE;     /* tree-incest occured */
    bool skipself          = FALSE;     /* self-interaction skipped */
    ForceEvalState fest = EMPTY_FORCE_EVAL_STATE;
    mwvector externalacc;
    bool intree = Mass(p) > 0.0;

    fest.pskip = p;           /* exclude p from f.c. */
    fest.pos0 = Pos(p);       /* set field point */

                  /* watch for tree-incest */
    skipself = treescan(ctx, &fest, root);         /* scan tree from t.root */
    if (intree && !skipself)            /* did tree-incest occur? */
    {
        if (!ctx->allowIncest) /* treat as catastrophic? */
            fail("hackgrav: tree-incest detected\n");
        if (!treeincest)           /* for the first time? */
            warn("\n[hackgrav: tree-incest detected]\n");
        treeincest = TRUE;          /* don't repeat warning */
    }

    /* Adding the external potential */
    externalacc = acceleration(ctx, Pos(p));

    mw_incaddv(fest.acc0, externalacc);

    /* TODO: Sharing */
    return fest.acc0;         /* and acceleration */
}

#ifndef _OPENMP

static inline void mapForceBody(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    mwvector* a;

    const bodyptr endp = st->bodytab + ctx->model.nbody;

    for (p = st->bodytab, a = st->acctab; p < endp; ++p, ++a)      /* get force on each body */
        *a = hackGrav(ctx, (nodeptr) st->tree.root, p);
}

#else

static inline void mapForceBody(const NBodyCtx* ctx, NBodyState* st)
{
    unsigned int i;
    const unsigned int nbody = ctx->model.nbody;

    #pragma omp parallel for private(i) schedule(dynamic)
    for (i = 0; i < nbody; ++i)      /* get force on each body */
    {
        st->acctab[i] = hackGrav(ctx,
                                 (nodeptr) st->tree.root,
                                 &st->bodytab[i]);
    }
}

#endif /* _OPENMP */

void gravMap(const NBodyCtx* ctx, NBodyState* st)
{
    makeTree(ctx, st);
    mapForceBody(ctx, st);
}

