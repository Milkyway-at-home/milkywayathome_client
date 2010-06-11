/* ************************************************************************** */
/* GRAV.C: routines to compute gravity. Public routines: hackgrav(). */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "defs.h"
#include "code.h"
#include <stdio.h>

static void treescan(nodeptr);           /* does force calculation */
static bool subdivp(cellptr);            /* can cell be accepted? */
static void gravsub(nodeptr);            /* compute grav interaction */

/* HACKGRAV: evaluate gravitational field on body p; checks to be
 * sure self-interaction was handled correctly if intree is true.
 */

static bodyptr pskip;                /* skip in force evaluation */
static vector pos0;              /* point to evaluate field */
static real phi0;                /* resulting potential */
static vector acc0;              /* resulting acceleration */
static bool skipself;                /* self-interaction skipped */
static bool treeincest = FALSE;          /* tree-incest occured */

void hackgrav(bodyptr p, bool intree)
{
    vector externalacc;
    pskip = p;                  /* exclude p from f.c. */
    SETV(pos0, Pos(p));             /* set field point */
    phi0 = 0.0;                 /* init total potential */
    CLRV(acc0);                 /* and total acceleration */
    ps.n2bterm = ps.nbcterm = 0;          /* count body & cell terms */
    skipself = FALSE;               /* watch for tree-incest */
    treescan((nodeptr) t.root);           /* scan tree from t.root */
    if (intree && !skipself)            /* did tree-incest occur? */
    {
        if (!ctx.allowIncest) /* treat as catastrophic? */
            error("hackgrav: tree-incest detected\n");
        if (!treeincest)           /* for the first time? */
            eprintf("\n[hackgrav: tree-incest detected]\n");
        treeincest = TRUE;          /* don't repeat warning */
    }

// Adding the external potential

    static real miya_ascal = 6.5;
    static real miya_bscal = 0.26;
    static real miya_mass = 4.45865888E5;
    static real plu_rc = 0.7;
    static real plu_mass = 1.52954402E5;
    static real vhalo = 73;
    static real haloq = 1.0;
    static real halod = 12.0;

    static real apar, qpar, spar, ppar, lpar, rpar, rcyl;

    rcyl = sqrt(Pos(p)[0] * Pos(p)[0] + Pos(p)[1] * Pos(p)[1]);
    qpar = sqrt(Pos(p)[2] * Pos(p)[2] + miya_bscal * miya_bscal);
    apar = miya_ascal + qpar;
    spar = (Pos(p)[0] * Pos(p)[0]) + (Pos(p)[1] * Pos(p)[1]) + ((miya_ascal + qpar) * (miya_ascal + qpar));
    ppar = sqrt ((Pos(p)[0] * Pos(p)[0]) + (Pos(p)[1] * Pos(p)[1]) + (Pos(p)[2] * Pos(p)[2])) + plu_rc;
    rpar = ppar - plu_rc;
    lpar = (rcyl * rcyl) + ((Pos(p)[2] / haloq) * (Pos(p)[2] / haloq)) + (halod * halod);

    phi0 -= (-(miya_mass) / sqrt(spar)) + (-(plu_mass) / ppar) + (vhalo * vhalo * log(lpar));

    externalacc[0] = - ( ( (2.0 * vhalo * vhalo * Pos(p)[0]) / (lpar) ) + ((plu_mass * Pos(p)[0]) / (rpar * ppar * ppar) ) + ((miya_mass * Pos(p)[0]) / (pow(spar, 1.5)) ) );
    externalacc[1] = - ( ( (2.0 * vhalo * vhalo * Pos(p)[1]) / (lpar) ) + ((plu_mass * Pos(p)[1]) / (rpar * ppar * ppar) ) + ((miya_mass * Pos(p)[1]) / (pow(spar, 1.5)) ) );
    externalacc[2] = - ( ( (2.0 * vhalo * vhalo * Pos(p)[2]) / (haloq * haloq * lpar) ) + ((plu_mass * Pos(p)[2]) / (rpar * ppar * ppar) ) + ( (miya_mass * Pos(p)[2] * apar) / (qpar * pow(spar, 1.5)) ) );

    ADDV(acc0, acc0, externalacc);

    Phi(p) = phi0;              /* store total potential */
    SETV(Acc(p), acc0);             /* and acceleration */
}

/*  * TREESCAN: iterative routine to do force calculation, starting
 * with node q, which is typically the t.root cell.
 */

static void treescan(nodeptr q)
{
    while (q != NULL)               /* while not at end of scan */
    {
        if (Type(q) == CELL &&          /* is node a cell and... */
                subdivp((cellptr) q))       /* too close to accept? */
            q = More(q);            /* follow to next level */
        else                    /* else accept this term */
        {
            if (q == (nodeptr) pskip)       /* self-interaction? */
                skipself = TRUE;        /* then just skip it */
            else                /* not self-interaction */
            {
                gravsub(q);                     /* so compute gravity */
                if (Type(q) == BODY)
                    ps.n2bterm++;          /* count body-body */
                else
                    ps.nbcterm++;          /* count body-cell */
            }
            q = Next(q);            /* follow next link */
        }
    }
}

/*  * SUBDIVP: decide if cell q is too close to accept as a single
 * term.  Also sets qmem, dr, and drsq for use by gravsub.
 */

static cellptr qmem;                         /* data shared with gravsub */
static vector dr;                /* vector from q to pos0 */
static real drsq;                /* squared distance to pos0 */

static bool subdivp(cellptr q)
{
    SUBV(dr, Pos(q), pos0);         /* compute displacement */
    DOTVP(drsq, dr, dr);            /* and find dist squared */
    qmem = q;                   /* remember we know them */
    return (drsq < Rcrit2(q));          /* apply standard rule */
}

/*  * GRAVSUB: compute contribution of node q to gravitational field at
 * point pos0, and add to running totals phi0 and acc0.
 */

static void gravsub(nodeptr q)
{
    real drab, phii, mor3;
    vector ai, quaddr;
    real dr5inv, phiquad, drquaddr;

    if (q != (nodeptr) qmem)                    /* cant use memorized data? */
    {
        SUBV(dr, Pos(q), pos0);                 /* then compute sep. */
        DOTVP(drsq, dr, dr);            /* and sep. squared */
    }
    drsq += ps.eps * ps.eps;                          /* use standard softening */
    drab = rsqrt(drsq);
    phii = Mass(q) / drab;
    mor3 = phii / drsq;
    MULVS(ai, dr, mor3);
    phi0 -= phii;                               /* add to total grav. pot. */
    ADDV(acc0, acc0, ai);                       /* ... and to total accel. */

    if (ctx.usequad && Type(q) == CELL)             /* if cell, add quad term */
    {
        dr5inv = 1.0 / (drsq * drsq * drab);    /* form dr^-5 */
        MULMV(quaddr, Quad(q), dr);             /* form Q * dr */
        DOTVP(drquaddr, dr, quaddr);            /* form dr * Q * dr */
        phiquad = -0.5 * dr5inv * drquaddr;     /* get quad. part of phi */
        phi0 = phi0 + phiquad;                  /* increment potential */
        phiquad = 5.0 * phiquad / drsq;         /* save for acceleration */
        MULVS(ai, dr, phiquad);                 /* components of acc. */
        SUBV(acc0, acc0, ai);                   /* increment */
        MULVS(quaddr, quaddr, dr5inv);
        SUBV(acc0, acc0, quaddr);               /* acceleration */
    }
}

