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
#include "nbody_tree.h"

#include <lua.h>
#include <lauxlib.h>
#include "nbody_lua_types.h"
#include "milkyway_util.h"

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif


/* subIndex: compute subcell index for body p in cell q. */
static inline int nbSubIndex(Body* p, NBodyCell* q)
{
    int ind = 0;

    /* accumulate subcell index */
    /* loop over dimensions */
    if (showRealValue(&X(&Pos(q))) <= showRealValue(&X(&Pos(p))))     /* if beyond midpoint */
        ind += NSUB >> (0 + 1);     /* skip over subcells */

    if (showRealValue(&Y(&Pos(q))) <= showRealValue(&Y(&Pos(p))))
        ind += NSUB >> (1 + 1);

    if (showRealValue(&Z(&Pos(q))) <= showRealValue(&Z(&Pos(p))))
        ind += NSUB >> (2 + 1);

    return ind;
}

static inline void nbIncAddNBodyQuadMatrix(NBodyQuadMatrix* RESTRICT a, NBodyQuadMatrix* RESTRICT b)
{
    a->xx = mw_add(&a->xx, &b->xx);
    a->xy = mw_add(&a->xy, &b->xy);
    a->xz = mw_add(&a->xz, &b->xz);

    a->yy = mw_add(&a->yy, &b->yy);
    a->yz = mw_add(&a->yz, &b->yz);

    a->zz = mw_add(&a->zz, &b->zz);
}

/* hackQuad: descend tree, evaluating quadrupole moments.  Note that this
 * routine is coded so that the Subp() and Quad() components of a cell can
 * share the same memory locations.
 */
static void hackQuad(NBodyCell* p)
{
    real tmp1, tmp2;
    unsigned int ndesc, i;
    NBodyNode* desc[NSUB];
    NBodyNode* q;
    mwvector dr;
    real drsq;
    NBodyQuadMatrix quad = { ZERO_REAL, ZERO_REAL, ZERO_REAL, ZERO_REAL, ZERO_REAL, ZERO_REAL };

    ndesc = 0;                                  /* count occupied subnodes  */
    for (i = 0; i < NSUB; ++i)                  /* loop over all subnodes   */
    {
        if (Subp(p)[i] != NULL)                 /* if this one's occupied   */
        {
            desc[ndesc++] = Subp(p)[i];         /* copy it to safety        */
        }
    }

    for (i = 0; i < ndesc; ++i)                 /* loop over real subnodes  */
    {
        q = desc[i];                            /* access each one in turn  */
        if (isCell(q))                          /* if it's also a cell      */
        {
            hackQuad((NBodyCell*) q);           /* then process it first    */
        }

        dr = mw_subv(&Pos(q), &Pos(p));           /* find displacement vect.  */
        drsq = mw_sqrv(&dr);                     /* and dot prod. (dr . dr)  */

        /* Outer product scaled by 3, then subtract drsq off the
         * diagonal to form quad moment*/
        {
            real m = Mass(q);   /* from CM of subnode */

            tmp1 = mw_mul(&X(&dr), &X(&dr));
            tmp1 = mw_mul_s(&tmp1, 3.0);
            tmp1 = mw_sub(&tmp1, &drsq);
            quad.xx = mw_mul(&m, &tmp1);

            tmp1 = mw_mul(&X(&dr), &Y(&dr));
            tmp1 = mw_mul_s(&tmp1, 3.0);
            quad.xy = mw_mul(&m, &tmp1);

            tmp1 = mw_mul(&X(&dr), &Z(&dr));
            tmp1 = mw_mul_s(&tmp1, 3.0);
            quad.xz = mw_mul(&m, &tmp1);

            tmp1 = mw_mul(&Y(&dr), &Y(&dr));
            tmp1 = mw_mul_s(&tmp1, 3.0);
            tmp1 = mw_sub(&tmp1, &drsq);
            quad.yy = mw_mul(&m, &tmp1);

            tmp1 = mw_mul(&Y(&dr), &Z(&dr));
            tmp1 = mw_mul_s(&tmp1, 3.0);
            quad.yz = mw_mul(&m, &tmp1);

            tmp1 = mw_mul(&Z(&dr), &Z(&dr));
            tmp1 = mw_mul_s(&tmp1, 3.0);
            tmp1 = mw_sub(&tmp1, &drsq);
            quad.zz = mw_mul(&m, &tmp1);
        }

        if (isCell(q)) /* if subnode is cell       */
        {
            nbIncAddNBodyQuadMatrix(&quad, &Quad(q));     /* then include its moment  */
        }

        nbIncAddNBodyQuadMatrix(&Quad(p), &quad); /* increment moment of cell */
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
            {
                desc[ndesc++] = Subp(p)[i];     /* store in table */
            }
        }
        More(p) = desc[0];                      /* link to first child */
        desc[ndesc] = n;                        /* end table with next */
        for (i = 0; i < ndesc; i++)             /* loop over children */
        {
            threadTree(desc[i], desc[i + 1]);     /* thread each w/ next */
        }
    }
}

/* expandBox: find range of coordinate values (with respect to root)
 * and expand root cell to fit. The size is doubled at each step to
 * take advantage of exact representation of powers of two.
 */
static void expandBox(NBodyTree* t, const Body* btab, int nbody)
{
    real_0 xyzmax;
    const Body* p;
    const NBodyCell* root = t->root;

    assert(t->rsize > 0.0);

    xyzmax = 0.0;
    for (p = btab; p < btab + nbody; ++p)
    {
        xyzmax = mw_fmax_0(xyzmax, mw_abs_0(showRealValue(&X(&Pos(p))) - showRealValue(&X(&Pos(root)))));
        xyzmax = mw_fmax_0(xyzmax, mw_abs_0(showRealValue(&Y(&Pos(p))) - showRealValue(&Y(&Pos(root)))));
        xyzmax = mw_fmax_0(xyzmax, mw_abs_0(showRealValue(&Z(&Pos(p))) - showRealValue(&Z(&Pos(root)))));
    }

    while (t->rsize < 2.0 * xyzmax)
    {
        t->rsize *= 2.0;
    }
}

/* makecell: return pointer to free cell. */
static NBodyCell* nbMakeCell(NBodyState* st, NBodyTree* t)
{
    NBodyCell* c;

    if (st->freeCell == NULL)                   /* no free cells left? */
    {
        c = (NBodyCell*) mwMallocA(sizeof(*c)); /* allocate a new one */
    }
    else                                        /* use existing free cell */
    {
        c = (NBodyCell*) st->freeCell;          /* take one on front */
        st->freeCell = Next(c);                 /* go on to next one */
    }
    Type(c) = CELL(0);                          /* initialize cell type */
    More(c) = NULL;
    memset(&c->stuff, 0, sizeof(c->stuff));     /* empty sub cells */
    t->cellUsed++;                              /* count one more cell */
    return c;
}

/* reclaim cells in tree, prepare to build new one. */
static void nbNewTree(NBodyState* st, NBodyTree* t)
{
    NBodyNode* p = (NBodyNode*) t->root;              /* start with the root */

    while (p != NULL)                       /* loop scanning tree */
    {
        if (isCell(p))                      /* found cell to free? */
        {
            Next(p) = st->freeCell;         /* link to front of */
            st->freeCell = p;               /* ...existing list */
            p = More(p);                    /* scan down tree */
        }
        else                                /* skip over bodies */
        {
            p = Next(p);                    /* go on to next */
        }
    }

    t->cellUsed = 0;   /* init count of cells, levels */
    t->maxDepth = 0;

    t->root = nbMakeCell(st, t);      /* allocate the root cell */
    mw_zerov(&Pos(t->root));           /* initialize the midpoint */
}


ALWAYS_INLINE
static inline real calcOffset(real* pPos, real* qPos, real_0 qsize)
{
    /* offset from parent */
    real_0 tmp_val = 0.25 * (showRealValue(pPos) < showRealValue(qPos) ? -qsize : qsize);
    return mw_add_s(qPos, tmp_val);
}

ALWAYS_INLINE
static inline void nbInitMidpoint(NBodyCell* c, const Body* p, const NBodyCell* q, real_0 qsize)
{
    Pos(c).x = calcOffset(&X(&Pos(p)), &X(&Pos(q)), qsize);
    Pos(c).y = calcOffset(&Y(&Pos(p)), &Y(&Pos(q)), qsize);
    Pos(c).z = calcOffset(&Z(&Pos(p)), &Z(&Pos(q)), qsize);
}

/* loadBody: descend tree and insert body p in appropriate place. */
static void nbLoadBody(NBodyState* st, NBodyTree* t, Body* p)
{
    NBodyCell* q;
    NBodyCell* c;
    size_t qind;
    size_t cind;
    unsigned int lev;
    real_0 qsize;
    //mw_printf("Body POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))) );
    //mw_printf("Body VEL = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&Vel(p))), showRealValue(&Y(&Vel(p))), showRealValue(&Z(&Vel(p))) );

    q = (NBodyCell*) t->root;                                /* start with tree t.root */
    qind = nbSubIndex(p, q);                    /* get index of subcell */
    qsize = t->rsize;                           /* keep track of cell size */
    lev = 0;                                    /* count levels descended */
    while (Subp(q)[qind] != NULL)               /* loop descending tree */
    {
        if (qsize <= REAL_EPSILON)
        {
            if (!t->structureError)
            {
                mw_printf("qsize (= %.15f) <= epsilon at level %u (initial root = %.15f)\n", qsize, lev, t->rsize);
                t->structureError = TRUE; /* FIXME: Not quite the same as the other structure error */
            }
            return;
        }

        if (isBody(Subp(q)[qind]))              /* reached a "leaf"? */
        {
            c = nbMakeCell(st, t);             /* allocate new cell */
            nbInitMidpoint(c, p, q, qsize);    /* initialize midpoint */

            Subp(c)[nbSubIndex((Body*) Subp(q)[qind], c)] = Subp(q)[qind];
            /* put body in cell */
            Subp(q)[qind] = (NBodyNode*) c;    /* link cell in tree */
        }
        q = (NBodyCell*) Subp(q)[qind];   /* advance to next level */
        qind = nbSubIndex(p, q);          /* get index to examine */
        qsize *= 0.5;                     /* shrink current cell */
        ++lev;                            /* count another level */
    }
    Subp(q)[qind] = (NBodyNode*) p;            /* found place, store p */
    t->maxDepth = MAX(t->maxDepth, lev);  /* remember maximum level */
}

ALWAYS_INLINE
static inline real_0 bmax2Inc(real_0 cmPos, real_0 pPos, real_0 psize)
{
    real_0 dmin;
    dmin = cmPos - (pPos - 0.5 * psize);         /* dist from 1st corner */
    return sqr_0(mw_fmax_0(dmin, psize - dmin));      /* sum max distance^2 */
}

ALWAYS_INLINE
static inline real_0 calcSW93MaxDist2(const NBodyCell* p, const mwvector* cmpos, real_0 psize)
{
    real_0 bmax2;

    /* compute max distance^2 */
    /* loop over dimensions */
    bmax2 = bmax2Inc(showRealValue(&X(cmpos)), showRealValue(&X(&Pos(p))), psize);
    bmax2 += bmax2Inc(showRealValue(&Y(cmpos)), showRealValue(&Y(&Pos(p))), psize);
    bmax2 += bmax2Inc(showRealValue(&Z(cmpos)), showRealValue(&Z(&Pos(p))), psize);

    return bmax2;
}

/* assign critical radius for cell p, using center-of-mass position
 * cmpos and cell size psize. */
static inline real_0 findRCrit(const NBodyCtx* ctx, const NBodyCell* p, real_0 treeRSize, mwvector* cmpos, real_0 psize)
{
    real_0 rc, bmax2;
    real tmp;
    mwvector diff;

    if (mw_unlikely(ctx->theta == 0.0))
    {
        /* Do an exact force calculation by always opening cells */
        rc = 2.0 * treeRSize;
        return sqr_0(rc);
    }

    /* return square of radius */
    switch (ctx->criterion)
    {
        case TreeCode:
            /* use size plus offset */
            diff = mw_subv(cmpos, &Pos(p));
            tmp = mw_length(&diff);
            if (showRealValue(&tmp) < 0.000001)
            {
                mw_printf("ZERO DISTANCE CALCULATED!\n");
                mw_printf("    P POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))) );
                mw_printf("    CMPOS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(cmpos)), showRealValue(&Y(cmpos)), showRealValue(&Z(cmpos)) );
            }
            rc = (psize / ctx->theta) + showRealValue(&tmp);
            return sqr_0(rc);

        case SW93:                           /* use S&W's criterion? */
            /* compute max distance^2 */
            bmax2 = calcSW93MaxDist2(p, cmpos, psize);
            return bmax2 / sqr_0(ctx->theta);      /* using max dist from cm */

        case BH86:                          /* use old BH criterion? */
            rc = psize / ctx->theta;        /* using size of cell */
            return sqr_0(rc);

        case InvalidCriterion:
        case Exact: /* Uses separate path */
        default:
            rc = 0.0; /* Stop clang static analysis warning */
            mw_fail("Invalid criterion: %s (%d)\n", showCriterionT(ctx->criterion), ctx->criterion);
	    return 0.0;
    }
}

static inline void nbCheckTreeDim(NBodyTree* tree, real_0 pPos, real_0 cmPos, real_0 halfPsize)
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
            mw_printf("hackCofM: tree structure error.\n"
                      "    cmPos out of bounds\n"
                      "    Pos(p)           = %.15e\n"
                      "    psize/2          = %.15e\n"
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

static inline void nbCheckTreeStructure(NBodyTree* tree, const mwvector* pPos, const mwvector* cmPos, const real_0 psize)
{
    real_0 halfPsize = 0.5 * psize;

    //mw_printf("P POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(pPos)), showRealValue(&Y(pPos)), showRealValue(&Z(pPos)) );
    //mw_printf("CMPOS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(cmPos)), showRealValue(&Y(cmPos)), showRealValue(&Z(cmPos)) );
    nbCheckTreeDim(tree, showRealValue(&X(pPos)), showRealValue(&X(cmPos)), halfPsize);
    nbCheckTreeDim(tree, showRealValue(&Y(pPos)), showRealValue(&Y(cmPos)), halfPsize);
    nbCheckTreeDim(tree, showRealValue(&Z(pPos)), showRealValue(&Z(cmPos)), halfPsize);
}


/* hackCofM: descend tree finding center-of-mass coordinates and
 * setting critical cell radii.
 */
static void hackCofM(const NBodyCtx* ctx, NBodyTree* tree, NBodyCell* p, real_0 psize)
{
    int i;
    NBodyNode* q;
    mwvector cmpos = ZERO_VECTOR;                /* init center of mass */
    mwvector tmp;

    assert(psize >= REAL_EPSILON);

    Mass(p) = ZERO_REAL;                              /* init total mass... */
    for (i = 0; i < NSUB; ++i)                  /* loop over subnodes */
    {
        if ((q = Subp(p)[i]) != NULL)           /* does subnode exist? */
        {
            if (isCell(q))                     /* and is it a cell? */
            {
                hackCofM(ctx, tree, (NBodyCell*) q, 0.5 * psize); /* find subcell cm */
            }

            Mass(p) = mw_add(&Mass(p), &Mass(q));                       /* sum total mass */
                                                      /* weight pos by mass */
            tmp = mw_mulvs(&Pos(q), &Mass(q));
            mw_incaddv(&cmpos, &tmp);     /* sum c-of-m position */
            //mw_printf("Mass(q) = %.15f\n", showRealValue(&Mass(q)));
            //mw_printf("Q POS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&Pos(q))), showRealValue(&Y(&Pos(q))), showRealValue(&Z(&Pos(q))) );
            //mw_printf("CMPOS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&cmpos)), showRealValue(&Y(&cmpos)), showRealValue(&Z(&cmpos)) );
        }
    }

    if (showRealValue(&Mass(p)) > 0.0)                          /* usually, cell has mass   */
    {
        cmpos = mw_divvs(&cmpos, &Mass(p));            /* so find c-of-m position  */
    }
    else                                        /* but if no mass inside    */
    {
        mw_printf("Found massless cell\n"); /* Debugging */
        cmpos = Pos(p);                /* use geo. center for now  */
    }

    //mw_printf("CMPOS = [%.15f, %.15f, %.15f]\n", showRealValue(&X(&cmpos)), showRealValue(&Y(&cmpos)), showRealValue(&Z(&cmpos)) );
    nbCheckTreeStructure(tree, &Pos(p), &cmpos, psize);

    Rcrit2(p) = findRCrit(ctx, p, tree->rsize, &cmpos, psize);            /* set critical radius */
    Pos(p) = cmpos;             /* and center-of-mass pos */
}

/* nbMakeTree: initialize tree structure for hierarchical force calculation
 * from body array btab, which contains ctx.nbody bodies.
 */
NBodyStatus nbMakeTree(const NBodyCtx* ctx, NBodyState* st)
{
    Body* p;
    const Body* endp = st->bodytab + st->nbody;
    NBodyTree* t = &st->tree;
    NBodyCell* q;

    nbNewTree(st, t);                                /* flush existing tree, etc */

    expandBox(t, st->bodytab, st->nbody);            /* and expand cell to fit */
    for (p = st->bodytab; p < endp; p++)             /* loop over bodies... */
    {
        if (showRealValue(&Mass(p)) != 0.0)                  /* exclude test particles */
            nbLoadBody(st, t, p);              /* and insert into tree */
    }

    /* Check if tree structure error occured */
    if (st->tree.structureError)
        return NBODY_TREE_STRUCTURE_ERROR;


    hackCofM(ctx, &st->tree, t->root, t->rsize);   /* find c-of-m coordinates */

    /* Check if tree structure error occured */
    if (st->tree.structureError)
        return NBODY_TREE_STRUCTURE_ERROR;

    threadTree((NBodyNode*) t->root, NULL);        /* add Next and More links */

    if (ctx->useQuad)                           /* including quad moments? */
        hackQuad(t->root);                      /* assign Quad moments */

    /*Code here is for debugging tree-incest errors*/
    mwbool checkIncest = FALSE;
    if (checkIncest)
    {
        for (p = st->bodytab; p < endp; p++)             /* loop over bodies... */
        {
            mwbool skipSelf = FALSE;
            q = (NBodyCell*) t->root; 
            while (q != NULL)               /* while not at end of scan */
            {
                if (isBody(q))      /* If is a body or far enough away to approximate */
                {
                    if (mw_likely((Body*) q == p))   /* self-interaction? */
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

                mw_printf("Incest found in nbody_tree.c at %p = [%.15f, %.15f, %.15f]\n", p, showRealValue(&X(&Pos(p))), showRealValue(&Y(&Pos(p))), showRealValue(&Z(&Pos(p))) );
            }
        }
    }
    /*----------------------------------------------*/

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

