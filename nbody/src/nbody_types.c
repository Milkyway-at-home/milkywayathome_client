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

#include "milkyway_util.h"
#include "nbody_types.h"
#include "nbody_show.h"

int destroyNBodyCtx(NBodyCtx* ctx)
{
    if (ctx->outfile && ctx->outfile != DEFAULT_OUTPUT_FILE)
    {
        if (fclose(ctx->outfile))
        {
            perror("closing output\n");
            return TRUE;
        }
    }

    return FALSE;
}

static void freeTree(Tree* t)
{
    Node* p;
    Node* tmp;

    p = (Node*) t->root;
    while (p != NULL)
    {
        if (isCell(p))
        {
            tmp = More(p);
            free(p);
            p = tmp;
        }
        else                        /* skip over bodies */
            p = Next(p);
    }

    t->root = NULL;
    t->cellused = 0;
    t->maxlevel = 0;
}

static void freeFreeCells(Node* freecell)
{
    Node* p;
    Node* tmp;

    p = freecell;
    while (p)
    {
        tmp = Next(p);
        free(p);
        p = tmp;
    }
}

void destroyNBodyState(NBodyState* st)
{
    freeTree(&st->tree);
    freeFreeCells(st->freecell);
    mwFreeA(st->bodytab);
    mwFreeA(st->acctab);

  #if NBODY_OPENCL
    cleanupNBodyCL(st);
  #endif /* NBODY_OPENCL */
}

void setInitialNBodyState(NBodyState* st, const NBodyCtx* ctx, Body* bodies, unsigned int nbody)
{
    static const Tree emptyTree = EMPTY_TREE;

    st->tree = emptyTree;
    st->freecell = NULL;

    st->tree.rsize = ctx->treeRSize;
    st->tnow = 0.0;
    st->nbody = nbody;
    st->bodytab = bodies;

    /* The tests may step the system from an arbitrary place, so make sure this is 0'ed */
    st->acctab = (mwvector*) mwCallocA(nbody, sizeof(mwvector));
}

NBodyState* newNBodyState()
{
    return mwCalloc(1, sizeof(NBodyState));
}

static int equalMaybeArray(const void* a, const void* b, size_t n)
{
    if (!a && !b)  /* Both not set, equal */
        return 1;

    if (!a || !b)  /* One is not set, not equal */
        return 0;

    assert(a && b);
    return memcmp(a, b, n);  /* Compare actual values */
}


#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wfloat-equal".
#endif

/* TODO: Doesn't handle tree */
/* Returns nonzero if states are equal, 0 otherwise */
int equalNBodyState(const NBodyState* st1, const NBodyState* st2)
{
    assert(st1 && st2);

    if (   st1->tnow != st2->tnow
        || st1->outputTime != st2->outputTime
        || st1->lastCheckpoint != st2->lastCheckpoint
        || st1->nbody != st2->nbody)
    {
        return 0;
    }

    if (!equalMaybeArray(st1->bodytab, st2->bodytab, st1->nbody * sizeof(Body)))
        return 1;

    return equalMaybeArray(st1->acctab, st2->acctab, st1->nbody * sizeof(mwvector));
}

/* TODO: Doesn't clone tree */
void cloneNBodyState(NBodyState* st, const NBodyState* oldSt, const unsigned int nbody)
{
    static const Tree emptyTree = EMPTY_TREE;

    st->tree = emptyTree;
    st->freecell = NULL;

    st->outputTime = oldSt->outputTime;
    st->lastCheckpoint = oldSt->lastCheckpoint;
    st->tnow = oldSt->tnow;

    st->bodytab = (Body*) mwMallocA(sizeof(Body) * nbody);
    memcpy(st->bodytab, oldSt->bodytab, sizeof(Body) * nbody);

    st->acctab = (mwvector*) mwMallocA(sizeof(mwvector) * nbody);
    memcpy(st->acctab, oldSt->acctab, sizeof(mwvector) * nbody);

    st->treeIncest = oldSt->treeIncest;
    st->incestReported = oldSt->incestReported;
}


static inline int compareComponents(real a, real b)
{
    if (a > b)
        return 1;
    if (a < b)
        return -1;

    return 0;
}

static int compareVectors(mwvector a, mwvector b)
{
    int rc;
    real ar, br;

    ar = mw_absv(a);
    br = mw_absv(b);

    if (ar > br)
        return 1;
    else if (ar < br)
        return -1;
    else
    {
        /* Resort to comparing by each component */
        if ((rc = compareComponents(X(a), X(b))))
            return rc;

        if ((rc = compareComponents(Y(a), Y(b))))
            return rc;

        if ((rc = compareComponents(Z(a), Z(b))))
            return rc;
    }

    return 0;  /* Equal */
}

/* Function for sorting bodies */
static int compareBodies(const void* _a, const void* _b)
{
    const Body* a = (const Body*) _a;
    const Body* b = (const Body*) _b;
    int rc;
    char* bufA;
    char* bufB;

    if ((rc = compareComponents(Mass(a), Mass(b))))
        return rc;

    /* Masses equal, compare positions */
    rc = compareVectors(Pos(a), Pos(b));
    if (rc == 0)
    {
        bufA = showBody(a);
        bufB = showBody(b);
        mw_panic("Comparing bodies with equal positions: %s, %s\n", bufA, bufB);
        free(bufA);  /* Never reached */
        free(bufB);
    }

    return rc;
}

/* Sort the bodies. The actual order doesn't matter, it just needs to
 * be consistent when we hash. This is so when if we end up shifting
 * bodies around for the GPU, the tests will still work as
 * expected. */
void sortBodies(Body* bodies, unsigned int nbody)
{
    qsort(bodies, (size_t) nbody, sizeof(Body), compareBodies);
}

