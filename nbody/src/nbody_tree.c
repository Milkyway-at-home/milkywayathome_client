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
#include "nbody_tree.h"

#include <lua.h>
#include <lauxlib.h>
#include "nbody_lua_types.h"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif


/* subIndex: compute subcell index for body p in cell q. */
static int subIndex(Body* p, NBodyCell* q)
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
static void hackQuad(NBodyCell* p)
{
    unsigned int ndesc, i;
    NBodyNode* desc[NSUB];
    NBodyNode* q;
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
        if (isCell(q))                          /* if it's also a cell      */
            hackQuad((NBodyCell*) q);                /* then process it first    */
        dr = mw_subv(Pos(q), Pos(p));           /* find displacement vect.  */
        mw_outv(drdr, dr, dr);                  /* form outer prod. of dr   */
        drsq = mw_sqrv(dr);                     /* and dot prod. dr * dr    */
        mw_set_matrix_identity(Idrsq);          /* init unit matrix         */
        mw_incmulms(Idrsq, drsq);               /* and scale by dr * dr     */
        mw_mulms(tmpm, drdr, 3.0);              /* scale drdr by 3          */
        mw_incsubm(tmpm, Idrsq);                /* now form quad. moment    */
        mw_incmulms(tmpm, Mass(q));             /* from cm of subnode       */
        if (isCell(q))                          /* if subnode is cell       */
            mw_incaddm(tmpm, Quad(q));          /* then include its moment  */
        mw_incaddm(Quad(p), tmpm);              /* increment moment of cell */
    }
}


/* threadTree: do a recursive treewalk starting from node p,
 * with next stop n, installing Next and More links.
 */
static void threadTree(NBodyNode* p, NBodyNode* n)
{
    unsigned int ndesc, i;
    NBodyNode* desc[NSUB+1];

    Next(p) = n;                                /* link to next node */
    if (isCell(p))                              /* any children to thread? */
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

/* expandBox: find range of coordinate values (with respect to root)
 * and expand root cell to fit. The size is doubled at each step to
 * take advantage of exact representation of powers of two.
 */
static void expandBox(NBodyTree* t, const Body* btab, unsigned int nbody)
{
    real xyzmax;
    const Body* p;

    const NBodyCell* root = t->root;

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

/* makecell: return pointer to free cell. */
static NBodyCell* makeCell(NBodyState* st, NBodyTree* t)
{
    NBodyCell* c;
    size_t i;

    if (st->freecell == NULL)                  /* no free cells left? */
        c = (NBodyCell*) mwMalloc(sizeof(*c)); /* allocate a new one */
    else                                       /* use existing free cell */
    {
        c = (NBodyCell*) st->freecell;             /* take one on front */
        st->freecell = Next(c);               /* go on to next one */
    }
    Type(c) = CELL(0);                          /* initialize cell type */
    for (i = 0; i < NSUB; i++)                  /* loop over subcells */
        Subp(c)[i] = NULL;                      /* and empty each one */
    t->cellused++;                              /* count one more cell */
    return c;
}

/* newTree: reclaim cells in tree, prepare to build new one. */
static void newTree(NBodyState* st, NBodyTree* t)
{
    NBodyNode* p = (NBodyNode*) t->root;              /* start with the root */

    while (p != NULL)                       /* loop scanning tree */
    {
        if (isCell(p))                      /* found cell to free? */
        {
            Next(p) = st->freecell;         /* link to front of */
            st->freecell = p;               /* ...existing list */
            p = More(p);                    /* scan down tree */
        }
        else                                /* skip over bodies */
            p = Next(p);                    /* go on to next */
    }

    t->cellused = 0;   /* init count of cells, levels */
    t->maxlevel = 0;

    t->root = makeCell(st, t);    /* allocate the root cell */
    mw_zerov(Pos(t->root));                          /* initialize the midpoint */
}


ALWAYS_INLINE
static inline real calcOffset(real pPos, real qPos, real qsize)
{
    /* offset from parent */
    return qPos + 0.25 * (pPos < qPos ? -qsize : qsize);
}

ALWAYS_INLINE
static inline void initMidpoint(NBodyCell* c, const Body* p, const NBodyCell* q, real qsize)
{
    X(Pos(c)) = calcOffset(X(Pos(p)), X(Pos(q)), qsize);
    Y(Pos(c)) = calcOffset(Y(Pos(p)), Y(Pos(q)), qsize);
    Z(Pos(c)) = calcOffset(Z(Pos(p)), Z(Pos(q)), qsize);
}

/* loadBody: descend tree and insert body p in appropriate place. */
static void loadBody(NBodyState* st, NBodyTree* t, Body* p)
{
    NBodyCell* q;
    NBodyCell* c;
    size_t qind;
    unsigned int lev;
    real qsize;

    q = t->root;                                /* start with tree t.root */
    qind = subIndex(p, q);                      /* get index of subcell */
    qsize = t->rsize;                           /* keep track of cell size */
    lev = 0;                                    /* count levels descended */
    while (Subp(q)[qind] != NULL)               /* loop descending tree */
    {
        if (isBody(Subp(q)[qind]))              /* reached a "leaf"? */
        {
            c = makeCell(st, t);               /* allocate new cell */
            initMidpoint(c, p, q, qsize);      /* initialize midpoint */

            Subp(c)[subIndex((Body*) Subp(q)[qind], c)] = Subp(q)[qind];
            /* put body in cell */
            Subp(q)[qind] = (NBodyNode*) c;    /* link cell in tree */
        }
        q = (NBodyCell*) Subp(q)[qind];        /* advance to next level */
        qind = subIndex(p, q);            /* get index to examine */
        qsize *= 0.5;                     /* shrink current cell */
        ++lev;                            /* count another level */
    }
    Subp(q)[qind] = (NBodyNode*) p;            /* found place, store p */
    t->maxlevel = MAX(t->maxlevel, lev);  /* remember maximum level */
}

ALWAYS_INLINE
static inline real bmax2Inc(real cmPos, real pPos, real psize)
{
    real dmin;
    dmin = cmPos - (pPos - 0.5 * psize);         /* dist from 1st corner */
    return sqr(mw_max(dmin, psize - dmin));      /* sum max distance^2 */
}

ALWAYS_INLINE
static inline real calcSW93MaxDist2(const NBodyCell* p, const mwvector cmpos, real psize)
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
static inline real findRCrit(const NBodyCtx* ctx, const NBodyCell* p, real treeRSize, mwvector cmpos, real psize)
{
    real rc, bmax2;

    switch (ctx->criterion)
    {
        case NewCriterion:
            rc = psize / ctx->theta + mw_distv(cmpos, Pos(p));
            /* use size plus offset */
            break;

        case SW93:                           /* use S&W's criterion? */
            /* compute max distance^2 */
            bmax2 = calcSW93MaxDist2(p, cmpos, psize);
            rc = mw_sqrt(bmax2) / ctx->theta;      /* using max dist from cm */
            break;

        case BH86:                          /* use old BH criterion? */
            rc = psize / ctx->theta;        /* using size of cell */
            break;

        case Exact:                       /* exact force calculation? */
            rc = 2.0 * treeRSize;         /* always open cells */
            break;

        case InvalidCriterion:
        default:
            rc = 0.0; /* Stop clang static analysis warning */
            fail("Bad criterion: %d\n", ctx->criterion);
    }

    return sqr(rc);          /* store square of radius */
}

static inline void checkTreeDim(NBodyTree* tree, const real pPos, const real cmPos, const real halfPsize)
{
    /* CHECKME: Precision: This gets angry as N gets big, and the divisions get small */
    if (   cmPos < pPos - halfPsize       /* if out of bounds */
        || cmPos > pPos + halfPsize)      /* in either direction */
    {
        if (!tree->structureError)
        {
            /* Only print if we don't know about the error
             * already. The error will be caught later and this
             * otherwise could print a lot of noise */
            tree->structureError = TRUE;
            warn("hackCofM: tree structure error.\n"
                 "    cmPos out of bounds\n"
                 "    Pos(p)           = %.15e\n"
                 "    tpsize/2          = %.15e\n"
                 "    Pos(p) + psize/2 = %.15e\n"
                 "    cmpos            = %.15e\n"
                 "    Pos(p) - psize/2 = %.15e\n",
                 pPos,
                 halfPsize,
                 pPos + halfPsize,
                 cmPos,
                 pPos - halfPsize);
        }
    }
}

static inline void checkTreeStructure(NBodyTree* tree, const mwvector pPos, const mwvector cmPos, const real psize)
{
    real halfPsize = 0.5 * psize;

    checkTreeDim(tree, X(pPos), X(cmPos), halfPsize);
    checkTreeDim(tree, Y(pPos), Y(cmPos), halfPsize);
    checkTreeDim(tree, Z(pPos), Z(cmPos), halfPsize);
}


/* hackCofM: descend tree finding center-of-mass coordinates and
 * setting critical cell radii.
 */
static void hackCofM(const NBodyCtx* ctx, NBodyTree* tree, NBodyCell* p, real psize)
{
    int i;
    NBodyNode* q;
    mwvector cmpos = ZERO_VECTOR;                /* init center of mass */

    assert(psize >= REAL_EPSILON);

    Mass(p) = 0.0;                              /* init total mass... */
    for (i = 0; i < NSUB; ++i)                  /* loop over subnodes */
    {
        if ((q = Subp(p)[i]) != NULL)           /* does subnode exist? */
        {
            if (isCell(q))                     /* and is it a cell? */
            {
                hackCofM(ctx, tree, (NBodyCell*) q, 0.5 * psize); /* find subcell cm */
            }

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

    checkTreeStructure(tree, Pos(p), cmpos, psize);

    Rcrit2(p) = findRCrit(ctx, p, tree->rsize, cmpos, psize);            /* set critical radius */
    Pos(p) = cmpos;             /* and center-of-mass pos */
}

/* makeTree: initialize tree structure for hierarchical force calculation
 * from body array btab, which contains ctx.nbody bodies.
 */
NBodyStatus makeTree(const NBodyCtx* ctx, NBodyState* st)
{
    Body* p;
    const Body* endp = st->bodytab + st->nbody;
    NBodyTree* t = &st->tree;

    newTree(st, t);                                  /* flush existing tree, etc */

    expandBox(t, st->bodytab, st->nbody);            /* and expand cell to fit */
    for (p = st->bodytab; p < endp; p++)             /* loop over bodies... */
    {
        if (Mass(p) != 0.0)                  /* exclude test particles */
            loadBody(st, t, p);              /* and insert into tree */
    }

    hackCofM(ctx, &st->tree, t->root, t->rsize);   /* find c-of-m coordinates */

    /* Check if tree structure error occured */
    if (st->tree.structureError)
        return NBODY_TREE_STRUCTURE_ERROR;

    threadTree((NBodyNode*) t->root, NULL);        /* add Next and More links */
    if (ctx->useQuad)                           /* including quad moments? */
        hackQuad(t->root);                      /* assign Quad moments */

    return NBODY_SUCCESS;
}

#if 0
/* For testing */
static int luaFindRCrit(lua_State* luaSt)
{
    const NBodyCtx* ctx;
    NBodyCell p;  /* Test cell, just need a set position */
    real rSize, pSize;
    mwvector cmPos;

    ctx = checkNBodyCtx(luaSt, 1);
    Pos(&p) = *checkVector(luaSt, 2);
    rSize = luaL_checknumber(luaSt, 3);
    cmPos = *checkVector(luaSt, 4);
    pSize = luaL_checknumber(luaSt, 5);

    lua_pushnumber(luaSt, findRCrit(ctx, &p, rSize, cmPos, pSize));

    return 1;
}

void registerFindRCrit(lua_State* luaSt)
{
    lua_register(luaSt, "findRCrit", luaFindRCrit);
}
#endif

