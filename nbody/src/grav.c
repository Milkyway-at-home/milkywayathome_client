/* ************************************************************************** */
/* grav.c: routines to compute gravity. Public routines: hackgrav(). */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "nbody_priv.h"

typedef struct
{
    vector pos0;      /* point to evaluate field */
    vector acc0;      /* resulting acceleration */
    bodyptr pskip;    /* skip in force evaluation */
    real phi0;        /* resulting potential */
    cellptr qmem;     /* data shared with gravsub */
    vector dr;        /* vector from q to pos0 */
    real drsq;        /* squared distance to pos0 */
} ForceEvalState;

#define EMPTY_FORCE_EVAL_STATE { ZERO_VECTOR, ZERO_VECTOR, NULL, 0.0, NULL, ZERO_VECTOR, 0.0 }


/* subdivp: decide if cell q is too close to accept as a single
 * term. Also sets qmem, dr, and drsq for use by gravsub.
 */
inline static bool subdivp(ForceEvalState* fest, cellptr q)
{
    SUBV(fest->dr, Pos(q), fest->pos0);   /* compute displacement */
    SQRV(fest->drsq, fest->dr);           /* and find dist squared */
    fest->qmem = q;                       /* remember we know them */
    return (fest->drsq < Rcrit2(q));      /* apply standard rule */
}

/* gravsub: compute contribution of node q to gravitational field at
 * point pos0, and add to running totals phi0 and acc0.
 */
inline static void gravsub(const NBodyCtx* ctx, ForceEvalState* fest, nodeptr q)
{
    real drab, phii, mor3;
    vector ai, quaddr;
    real dr5inv, phiquad, drquaddr;

    if (q != (nodeptr) fest->qmem)                    /* cant use memorized data? */
    {
        SUBV(fest->dr, Pos(q), fest->pos0);           /* then compute sep. */
        SQRV(fest->drsq, fest->dr);                   /* and sep. squared */
    }
    fest->drsq += sqr(ctx->model.eps);              /* use standard softening */
    drab = rsqrt(fest->drsq);
    phii = Mass(q) / drab;
    mor3 = phii / fest->drsq;
    MULVS(ai, fest->dr, mor3);
    fest->phi0 -= phii;                         /* add to total grav. pot. */
    INCADDV(fest->acc0, ai);                   /* ... and to total accel. */

    if (ctx->usequad && Type(q) == CELL)             /* if cell, add quad term */
    {
        dr5inv = 1.0 / (sqr(fest->drsq) * drab); /* form dr^-5 */
        MULMV(quaddr, Quad(q), fest->dr);        /* form Q * dr */
        DOTVP(drquaddr, fest->dr, quaddr);       /* form dr * Q * dr */
        phiquad = -0.5 * dr5inv * drquaddr;      /* get quad. part of phi */
        fest->phi0 += phiquad;                   /* increment potential */
        phiquad = 5.0 * phiquad / fest->drsq;    /* save for acceleration */
        MULVS(ai, fest->dr, phiquad);            /* components of acc. */
        INCSUBV(fest->acc0, ai);                 /* increment */
        INCMULVS(quaddr, dr5inv);
        INCSUBV(fest->acc0, quaddr);             /* acceleration */
    }
}

/* treescan: iterative routine to do force calculation, starting with
 * node q, which is typically the t.root cell. Watches for tree
 * incest.
 */
inline static bool treescan(const NBodyCtx* ctx,
                            NBodyState* st,
                            ForceEvalState* fest,
                            nodeptr q)
{
    bool skipself = FALSE;

    while (q != NULL)               /* while not at end of scan */
    {
        if (Type(q) == CELL &&                /* is node a cell and... */
            subdivp(fest, (cellptr) q))       /* too close to accept? */
            q = More(q);            /* follow to next level */
        else                    /* else accept this term */
        {
            if (q == (nodeptr) fest->pskip)    /* self-interaction? */
                skipself = TRUE;               /* then just skip it */
            else                               /* not self-interaction */
            {
                gravsub(ctx, fest, q);         /* so compute gravity */
                if (Type(q) == BODY)
                    st->n2bterm++;          /* count body-body */
                else
                    st->nbcterm++;          /* count body-cell */
            }
            q = Next(q);            /* follow next link */
        }
    }

    return skipself;
}

/* hackgrav: evaluate gravitational field on body p; checks to be
 * sure self-interaction was handled correctly if intree is true.
 */
void hackgrav(const NBodyCtx* ctx, NBodyState* st, bodyptr p, bool intree)
{
    vector externalacc;
    int n2bterm;    /* number 2-body of terms evaluated */
    int nbcterm;    /* num of body-cell terms evaluated */
    static bool treeincest = FALSE;     /* tree-incest occured */
    bool skipself          = FALSE;     /* self-interaction skipped */

    ForceEvalState fest = EMPTY_FORCE_EVAL_STATE;


    fest.pskip = p;                /* exclude p from f.c. */
    SETV(fest.pos0, Pos(p));       /* set field point */

    n2bterm = nbcterm = 0;          /* count body & cell terms */
                  /* watch for tree-incest */
    skipself = treescan(ctx, st, &fest, (nodeptr) st->tree.root);         /* scan tree from t.root */
    if (intree && !skipself)            /* did tree-incest occur? */
    {
        if (!ctx->allowIncest) /* treat as catastrophic? */
            error("hackgrav: tree-incest detected\n");
        if (!treeincest)           /* for the first time? */
            eprintf("\n[hackgrav: tree-incest detected]\n");
        treeincest = TRUE;          /* don't repeat warning */
    }

    /* Adding the external potential */
    acceleration(externalacc, ctx, Pos(p));

    INCADDV(fest.acc0, externalacc);

    /* TODO: Sharing */
    /* CHECKME: Only the acceleration here seems to have an effect on
     * the results */

    fest.phi0 -=   miyamotoNagaiPhi(&ctx->pot.disk, Pos(p))
                 + sphericalPhi(&ctx->pot.sphere[0], Pos(p))
                 + logHaloPhi(&ctx->pot.halo, Pos(p));

    Phi(p) = fest.phi0;              /* store total potential */
    SETV(Acc(p), fest.acc0);         /* and acceleration */

}

