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

#ifndef _WIN32
  #include <sys/stat.h>
  #include <sys/mman.h>
  #include <unistd.h>
#else
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#endif /* _WIN32 */

#include "nbody_priv.h"
#include "load.h"
#include "milkyway_util.h"

/* subIndex: compute subcell index for body p in cell q. */
static int subIndex(bodyptr p, cellptr q)
{
    int ind = 0;

    /* accumulate subcell index */
    /* loop over dimensions */
    if (X(Pos(q)) <= X(Pos(p)))     /* if beyond midpoint */
        ind += NSUB >> (0 + 1);     /* skip over subcells */

    if (Y(Pos(q)) <= Y(Pos(p)))
        ind += NSUB >> (1 + 1);

    if (Z(Pos(q)) <= Z(Pos(p)))
        ind += NSUB >> (2 + 1);

    return ind;
}

/* hackQuad: descend tree, evaluating quadrupole moments.  Note that this
 * routine is coded so that the Subp() and Quad() components of a cell can
 * share the same memory locations.
 */

/* TODO: Incremental matrix operations */
static inline void hackQuad(cellptr p)
{
    unsigned int ndesc, i;
    nodeptr desc[NSUB], q;
    mwvector dr;
    real drsq;
    mwmatrix drdr, Idrsq, tmpm;

    ndesc = 0;                                  /* count occupied subnodes  */
    for (i = 0; i < NSUB; ++i)                  /* loop over all subnodes   */
    {
        if (Subp(p)[i] != NULL)                 /* if this one's occupied   */
            desc[ndesc++] = Subp(p)[i];         /* copy it to safety        */
    }

    mw_set_matrix_zero(Quad(p));                /* init quadrupole moment   */
    for (i = 0; i < ndesc; ++i)                 /* loop over real subnodes  */
    {
        q = desc[i];                            /* access each one in turn  */
        if (Type(q) == CELL)                    /* if it's also a cell      */
            hackQuad((cellptr) q);              /* then process it first    */
        dr = mw_subv(Pos(q), Pos(p));           /* find displacement vect.  */
        mw_outv(drdr, dr, dr);                  /* form outer prod. of dr   */
        drsq = mw_sqrv(dr);                     /* and dot prod. dr * dr    */
        mw_set_matrix_identity(Idrsq);          /* init unit matrix         */
        mw_incmulms(Idrsq, drsq);               /* and scale by dr * dr     */
        mw_mulms(tmpm, drdr, 3.0);              /* scale drdr by 3          */
        mw_incsubm(tmpm, Idrsq);                /* now form quad. moment    */
        mw_incmulms(tmpm, Mass(q));             /* from cm of subnode       */
        if (Type(q) == CELL)                    /* if subnode is cell       */
            mw_incaddm(tmpm, Quad(q));          /* then include its moment  */
        mw_incaddm(Quad(p), tmpm);              /* increment moment of cell */
    }
}


/* threadTree: do a recursive treewalk starting from node p,
 * with next stop n, installing Next and More links.
 */
static void threadTree(nodeptr p, nodeptr n)
{
    unsigned int ndesc, i;
    nodeptr desc[NSUB+1];

    Next(p) = n;                                /* link to next node */
    if (Type(p) == CELL)                        /* any children to thread? */
    {
        ndesc = 0;                              /* count extant children */
        for (i = 0; i < NSUB; ++i)              /* loop over subnodes */
        {
            if (Subp(p)[i] != NULL)             /* found a live one? */
                desc[ndesc++] = Subp(p)[i];     /* store in table */
        }
        More(p) = desc[0];                      /* link to first child */
        desc[ndesc] = n;                        /* end table with next */
        for (i = 0; i < ndesc; i++)             /* loop over children */
            threadTree(desc[i], desc[i+1]);     /* thread each w/ next */
    }
}

/* expandBox: find range of coordinate values (with respect to t.root)
 * and expand t.root cell to fit.  The size is doubled at each step to
 * take advantage of exact representation of powers of two.
 */
static void expandBox(Tree* t, bodyptr btab, unsigned int nbody)
{
    real xyzmax;
    bodyptr p;

    const cellptr root = t->root;

    xyzmax = 0.0;
    for (p = btab; p < btab + nbody; ++p)
    {
        xyzmax = mw_max(xyzmax, mw_abs(X(Pos(p)) - X(Pos(root))));
        xyzmax = mw_max(xyzmax, mw_abs(Y(Pos(p)) - Y(Pos(root))));
        xyzmax = mw_max(xyzmax, mw_abs(Z(Pos(p)) - Z(Pos(root))));
    }

    while (t->rsize < 2.0 * xyzmax)
        t->rsize *= 2.0;
}

/* newTree: reclaim cells in tree, prepare to build new one. */
static void newTree(NBodyState* st, Tree* t)
{
    static bool firstcall = TRUE;
    nodeptr p;

    if (!firstcall)                             /* tree data to reclaim? */
    {
        p = (nodeptr) t->root;                  /* start with the t.root */
        while (p != NULL)                       /* loop scanning tree */
        {
            if (Type(p) == CELL)                /* found cell to free? */
            {
                Next(p) = st->freecell;         /* link to front of */
                st->freecell = p;               /* ...existing list */
                p = More(p);                    /* scan down tree */
            }
            else                                /* skip over bodies */
                p = Next(p);                    /* go on to next */
        }
    }
    else                                        /* first time through */
        firstcall = FALSE;                      /* so just note it */
    t->root = NULL;                             /* flush existing tree */
    t->cellused = 0;                            /* reset cell count */
}

/* makecell: return pointer to free cell. */
static cellptr makeCell(NBodyState* st, Tree* t)
{
    cellptr c;
    size_t i;

    if (st->freecell == NULL)                    /* no free cells left? */
        c = (cellptr) mallocSafe(sizeof(cell));  /* allocate a new one */
    else                                         /* use existing free cell */
    {
        c = (cellptr) st->freecell;             /* take one on front */
        st->freecell = Next(c);                 /* go on to next one */
    }
    Type(c) = CELL;                             /* initialize cell type */
    for (i = 0; i < NSUB; i++)                  /* loop over subcells */
        Subp(c)[i] = NULL;                      /* and empty each one */
    t->cellused++;                              /* count one more cell */
    return c;
}

ALWAYS_INLINE
static inline real calcOffset(real pPos, real qPos, real qsize)
{
    /* offset from parent */
    return qPos + 0.25 * (pPos < qPos ? -qsize : qsize);
}

ALWAYS_INLINE
static inline void initMidpoint(cellptr c, const bodyptr p, const cellptr q, real qsize)
{
    X(Pos(c)) = calcOffset(X(Pos(p)), X(Pos(q)), qsize);
    Y(Pos(c)) = calcOffset(Y(Pos(p)), Y(Pos(q)), qsize);
    Z(Pos(c)) = calcOffset(Z(Pos(p)), Z(Pos(q)), qsize);
}

/* loadBody: descend tree and insert body p in appropriate place. */
static void loadBody(NBodyState* st, Tree* t, bodyptr p)
{
    cellptr q, c;
    size_t qind;
    unsigned int lev;
    real qsize;

    q = t->root;                                /* start with tree t.root */
    qind = subIndex(p, q);                      /* get index of subcell */
    qsize = t->rsize;                           /* keep track of cell size */
    lev = 0;                                    /* count levels descended */
    while (Subp(q)[qind] != NULL)               /* loop descending tree */
    {
        if (Type(Subp(q)[qind]) == BODY)        /* reached a "leaf"? */
        {
            c = makeCell(st, t);               /* allocate new cell */
            initMidpoint(c, p, q, qsize);      /* initialize midpoint */

            Subp(c)[subIndex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
            /* put body in cell */
            Subp(q)[qind] = (nodeptr) c;        /* link cell in tree */
        }
        q = (cellptr) Subp(q)[qind];        /* advance to next level */
        qind = subIndex(p, q);              /* get index to examine */
        qsize /= 2;                         /* shrink current cell */
        ++lev;                              /* count another level */
    }
    Subp(q)[qind] = (nodeptr) p;            /* found place, store p */
    t->maxlevel = MAX(t->maxlevel, lev);    /* remember maximum level */
}

ALWAYS_INLINE
static inline real bmax2Inc(real cmPos, real pPos, real psize)
{
    real dmin;
    dmin = cmPos - (pPos - 0.5 * psize);         /* dist from 1st corner */
    return sqr(mw_max(dmin, psize - dmin));      /* sum max distance^2 */
}

ALWAYS_INLINE
static inline real calcSW93MaxDist2(const cellptr p, const mwvector cmpos, real psize)
{
    real bmax2;

    /* compute max distance^2 */
    /* loop over dimensions */
    bmax2 = bmax2Inc(X(cmpos), X(Pos(p)), psize);
    bmax2 += bmax2Inc(Y(cmpos), Y(Pos(p)), psize);
    bmax2 += bmax2Inc(Z(cmpos), Z(Pos(p)), psize);

    return bmax2;
}

/* setRCrit: assign critical radius for cell p, using center-of-mass
 * position cmpos and cell size psize. */
static void setRCrit(const NBodyCtx* ctx, NBodyState* st, cellptr p, mwvector cmpos, real psize)
{
    real rc, bmax2;

    switch (ctx->criterion)
    {
        case NEWCRITERION:
            rc = psize / ctx->theta + mw_distv(cmpos, Pos(p));
            /* use size plus offset */
            break;
        case EXACT:                         /* exact force calculation? */
            rc = 2.0 * st->tree.rsize;      /* always open cells */
            break;
        case BH86:                          /* use old BH criterion? */
            rc = psize / ctx->theta;        /* using size of cell */
            break;
        case SW93:                           /* use S&W's criterion? */
            /* compute max distance^2 */
            bmax2 = calcSW93MaxDist2(p, cmpos, psize);
            rc = mw_sqrt(bmax2) / ctx->theta;      /* using max dist from cm */
            break;
        default:
            rc = 0.0; /* Stop clang static analysis warning */
            fail("Got unknown criterion: %d\n", ctx->criterion);
    }

    Rcrit2(p) = sqr(rc);           /* store square of radius */
}

static inline void checkTreeDim(const real pPos, const real cmPos, const real halfPsize)
{
    /* CHECKME: Precision: This gets angry as N gets big, and the divisions get small */
    if (   cmPos < pPos - halfPsize       /* if out of bounds */
        || cmPos > pPos + halfPsize)      /* in either direction */
    {
        warn("hackCofM: tree structure error.\n"
             "\tcmpos out of bounds\n"
             "\tPos(p)           = %e\n"
             "\tpsize/2          = %e\n"
             "\tPos(p) + psize/2 = %e\n"
             "\tcmpos            = %e\n"
             "\tPos(p) - psize/2 = %e\n",
             pPos,
             halfPsize,
             pPos + halfPsize,
             cmPos,
             pPos - halfPsize);
    }
}

static inline void checkTreeStructure(const mwvector pPos, const mwvector cmPos, const real psize)
{
    real halfPsize = 0.5 * psize;

    checkTreeDim(X(pPos), X(cmPos), halfPsize);
    checkTreeDim(Y(pPos), Y(cmPos), halfPsize);
    checkTreeDim(Z(pPos), Z(cmPos), halfPsize);
}


/* hackCofM: descend tree finding center-of-mass coordinates and
 * setting critical cell radii.
 */
static void hackCofM(const NBodyCtx* ctx, NBodyState* st, cellptr p, real psize)
{
    int i;
    nodeptr q;
    mwvector cmpos = ZERO_VECTOR;                 /* init center of mass */

    Mass(p) = 0.0;                              /* init total mass... */
    for (i = 0; i < NSUB; ++i)                  /* loop over subnodes */
    {
        if ((q = Subp(p)[i]) != NULL)           /* does subnode exist? */
        {
            if (Type(q) == CELL)                /* and is it a cell? */
                hackCofM(ctx, st, (cellptr) q, 0.5 * psize); /* find subcell cm */
            Mass(p) += Mass(q);                       /* sum total mass */
                                                      /* weight pos by mass */
            mw_incaddv_s(cmpos, Pos(q), Mass(q));     /* sum c-of-m position */
        }
    }

    if (Mass(p) > 0.0)                          /* usually, cell has mass   */
    {
        mw_incdivs(cmpos, Mass(p));            /* so find c-of-m position  */
    }
    else                                        /* but if no mass inside    */
    {
        warn("Found massless cell\n"); /* Debugging */
        cmpos = Pos(p);                /* use geo. center for now  */
    }

    checkTreeStructure(Pos(p), cmpos, psize);
    setRCrit(ctx, st, p, cmpos, psize);            /* set critical radius */
    Pos(p) = cmpos;             /* and center-of-mass pos */
}

/* makeTree: initialize tree structure for hierarchical force calculation
 * from body array btab, which contains ctx.nbody bodies.
 */
void makeTree(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    const bodyptr endp = st->bodytab + ctx->nbody;
    Tree* t = &st->tree;

    newTree(st, t);                                  /* flush existing tree, etc */
    t->root = makeCell(st, t);                       /* allocate the t.root cell */
    mw_zerov(Pos(t->root));                          /* initialize the midpoint */
    expandBox(t, st->bodytab, ctx->nbody);           /* and expand cell to fit */
    t->maxlevel = 0;                                 /* init count of levels */
    for (p = st->bodytab; p < endp; p++)             /* loop over bodies... */
    {
        if (Mass(p) != 0.0)                  /* exclude test particles */
            loadBody(st, t, p);              /* and insert into tree */
    }

    hackCofM(ctx, st, t->root, t->rsize);       /* find c-of-m coordinates */
    threadTree((nodeptr) t->root, NULL);        /* add Next and More links */
    if (ctx->usequad)                           /* including quad moments? */
        hackQuad(t->root);                      /* assign Quad moments */
}

