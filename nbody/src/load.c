/* ************************************************************************** */
/* load.C: routines to create tree.  Public routines: makeTree(). */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _WIN32
  #include <sys/stat.h>
  #include <sys/mman.h>
  #include <unistd.h>
#else
  #include <windows.h>
#endif /* _WIN32 */

#include "nbody_priv.h"
#include "load.h"

/* subIndex: compute subcell index for body p in cell q. */
static int subIndex(bodyptr p, cellptr q)
{
    size_t k, ind = 0;

    /* accumulate subcell index */
    for (k = 0; k < NDIM; ++k)          /* loop over dimensions */
    {
        if (Pos(q)[k] <= Pos(p)[k])     /* if beyond midpoint */
            ind += NSUB >> (k + 1);     /* skip over subcells */
    }
    return ind;
}

/* hackQuad: descend tree, evaluating quadrupole moments.  Note that this
 * routine is coded so that the Subp() and Quad() components of a cell can
 * share the same memory locations.
 */

/* TODO: Incremental matrix operations */
inline static void hackQuad(cellptr p)
{
    unsigned int ndesc, i;
    nodeptr desc[NSUB], q;
    vector dr;
    real drsq;
    matrix drdr, Idrsq, tmpm;

    ndesc = 0;                                  /* count occupied subnodes  */
    for (i = 0; i < NSUB; ++i)                  /* loop over all subnodes   */
    {
        if (Subp(p)[i] != NULL)                 /* if this one's occupied   */
            desc[ndesc++] = Subp(p)[i];         /* copy it to safety        */
    }
    CLRM(Quad(p));                              /* init quadrupole moment   */
    for (i = 0; i < ndesc; ++i)                 /* loop over real subnodes  */
    {
        q = desc[i];                            /* access ech one in turn   */
        if (Type(q) == CELL)                    /* if it's also a cell      */
            hackQuad((cellptr) q);              /* then process it first    */
        SUBV(dr, Pos(q), Pos(p));               /* find displacement vect.  */
        OUTVP(drdr, dr, dr);                    /* form outer prod. of dr   */
        SQRV(drsq, dr);                         /* and dot prod. dr * dr    */
        SETMI(Idrsq);                           /* init unit matrix         */
        MULMS(Idrsq, Idrsq, drsq);              /* and scale by dr * dr     */
        MULMS(tmpm, drdr, 3.0);                 /* scale drdr by 3          */
        SUBM(tmpm, tmpm, Idrsq);                /* now form quad. moment    */
        MULMS(tmpm, tmpm, Mass(q));             /* from cm of subnode       */
        if (Type(q) == CELL)                    /* if subnode is cell       */
            ADDM(tmpm, tmpm, Quad(q));          /* then include its moment  */
        ADDM(Quad(p), Quad(p), tmpm);           /* increment moment of cell */
    }
}


/* threadTree: do a recursive treewalk starting from node p,
 * with next stop n, installing Next and More links.
 */
static void threadTree(nodeptr p, nodeptr n)
{
    int ndesc, i;
    nodeptr desc[NSUB+1];

    Next(p) = n;                                /* link to next node */
    if (Type(p) == CELL)                        /* any children to thread? */
    {
        ndesc = 0;                              /* count extant children */
        for (i = 0; i < NSUB; i++)              /* loop over subnodes */
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
static void expandBox(Tree* t, bodyptr btab, int nbody)
{
    real xyzmax;
    bodyptr p;
    size_t k;

    const cellptr root = t->root;

    xyzmax = 0.0;
    for (p = btab; p < btab + nbody; ++p)
    {
        for (k = 0; k < NDIM; ++k)
            xyzmax = MAX(xyzmax, rabs(Pos(p)[k] - Pos(root)[k]));
    }

    while (t->rsize < 2 * xyzmax)
        t->rsize *= 2;
}

/* newTree: reclaim cells in tree, prepare to build new one. */
static nodeptr freecell = NULL;              /* list of free cells */

static void newTree(Tree* t)
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
                Next(p) = freecell;             /* link to front of */
                freecell = p;                   /* ...existing list */
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
static cellptr makeCell(Tree* t)
{
    cellptr c;
    size_t i;

    if (freecell == NULL)                        /* no free cells left? */
        c = (cellptr) mallocSafe(sizeof(cell));  /* allocate a new one */
    else                                         /* use existing free cell */
    {
        c = (cellptr) freecell;                 /* take one on front */
        freecell = Next(c);                     /* go on to next one */
    }
    Type(c) = CELL;                             /* initialize cell type */
    for (i = 0; i < NSUB; i++)                  /* loop over subcells */
        Subp(c)[i] = NULL;                      /* and empty each one */
    t->cellused++;                              /* count one more cell */
    return c;
}

/* loadBody: descend tree and insert body p in appropriate place. */
static void loadBody(Tree* t, bodyptr p)
{
    cellptr q, c;
    size_t qind, k;
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
            c = makeCell(t);                    /* allocate new cell */
            for (k = 0; k < NDIM; k++)          /* initialize midpoint */
            {
                Pos(c)[k] = Pos(q)[k] +         /* offset from parent */
                            (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;
            }
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

/* setRCrit: assign critical radius for cell p, using center-of-mass
 * position cmpos and cell size psize. */
static void setRCrit(const NBodyCtx* ctx, NBodyState* st, cellptr p, vector cmpos, real psize)
{
    real rc, bmax2, dmin, tmp;
    size_t k;

    switch (ctx->criterion)
    {
        case NEWCRITERION:
            DISTV(tmp, cmpos, Pos(p));
            rc = psize / ctx->theta + tmp;
            /* use size plus offset */
            break;
        case EXACT:                         /* exact force calculation? */
            rc = 2 * st->tree.rsize;        /* always open cells */
            break;
        case BH86:                          /* use old BH criterion? */
            rc = psize / ctx->theta;        /* using size of cell */
            break;
        case SW93:                           /* use S&W's criterion? */
            bmax2 = 0.0;                     /* compute max distance^2 */
            for (k = 0; k < NDIM; ++k)       /* loop over dimensions */
            {
                dmin = cmpos[k] - (Pos(p)[k] - psize / 2);
                /* dist from 1st corner */
                bmax2 += sqr(MAX(dmin, psize - dmin));
                /* sum max distance^2 */
            }
            rc = rsqrt(bmax2) / ctx->theta;      /* using max dist from cm */
            break;
        default:
            fail("Got unknown criterion: %d\n", ctx->criterion);
    }

    Rcrit2(p) = sqr(rc);           /* store square of radius */
}

/* hackCofM: descend tree finding center-of-mass coordinates and
 * setting critical cell radii.
 */
static void hackCofM(const NBodyCtx* ctx, NBodyState* st, cellptr p, real psize)
{
    int i, k;
    nodeptr q;
    vector tmpv;
    vector cmpos = ZERO_VECTOR;                 /* init center of mass */

    Mass(p) = 0.0;                              /* init total mass... */
    for (i = 0; i < NSUB; ++i)                  /* loop over subnodes */
    {
        if ((q = Subp(p)[i]) != NULL)           /* does subnode exist? */
        {
            if (Type(q) == CELL)                /* and is it a cell? */
                hackCofM(ctx, st, (cellptr) q, psize / 2); /* find subcell cm */
            Mass(p) += Mass(q);                 /* sum total mass */
            MULVS(tmpv, Pos(q), Mass(q));       /* weight pos by mass */
            INCADDV(cmpos, tmpv);               /* sum c-of-m position */
        }
    }

    if (Mass(p) > 0.0)                          /* usually, cell has mass   */
    {
        INCDIVVS(cmpos, Mass(p));               /* so find c-of-m position  */
    }
    else                                        /* but if no mass inside    */
    {
        warn("Found massless cell\n"); /* Debugging */
        SETV(cmpos, Pos(p));                    /* use geo. center for now  */
    }

    for (k = 0; k < NDIM; k++)          /* check tree structure... */
    {
        /* CHECKME: Precision: This gets angry as N gets big, and the divisions get small */
        if (cmpos[k] < Pos(p)[k] - psize / 2 ||    /* if out of bounds */
                Pos(p)[k] + psize / 2 < cmpos[k])  /* in either direction */
        {
            warn("hackCofM: tree structure error.\n"
                 "\tcmpos out of bounds\n"
                 "\tPos(p)[%d]           = %e\n"
                 "\tpsize               = %e\n"
                 "\tPos(p)[%d] + psize/2 = %e\n"
                 "\tcmpos[%d]            = %e\n"
                 "\tPos(p)[%d] - psize/2 = %e\n",
                 k, Pos(p)[k],
                 psize,
                 k, Pos(p)[k] + psize / 2,
                 k, cmpos[k],
                 k, Pos(p)[k] - psize / 2);
        }
    }
    setRCrit(ctx, st, p, cmpos, psize);            /* set critical radius */
    SETV(Pos(p), cmpos);            /* and center-of-mass pos */
}

/* makeTree: initialize tree structure for hierarchical force calculation
 * from body array btab, which contains ctx.nbody bodies.
 */
void makeTree(const NBodyCtx* ctx, NBodyState* st)
{
    bodyptr p;
    const bodyptr endp = st->bodytab + ctx->model.nbody;
    Tree* t = &st->tree;

    newTree(t);                                      /* flush existing tree, etc */
    t->root = makeCell(t);                           /* allocate the t.root cell */
    CLRV(Pos(t->root));                              /* initialize the midpoint */
    expandBox(t, st->bodytab, ctx->model.nbody);     /* and expand cell to fit */
    t->maxlevel = 0;                                 /* init count of levels */
    for (p = st->bodytab; p < endp; p++)             /* loop over bodies... */
    {
        if (Mass(p) != 0.0)                     /* exclude test particles */
            loadBody(t, p);                     /* and insert into tree */
    }

    hackCofM(ctx, st, t->root, t->rsize);       /* find c-of-m coordinates */
    threadTree((nodeptr) t->root, NULL);        /* add Next and More links */
    if (ctx->usequad)                           /* including quad moments? */
        hackQuad(t->root);                      /* assign Quad moments */

    //printf("Cells used: %d, max height = %d, rsize = %g\n", t->cellused, t->maxlevel, t->rsize);
}

