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

static void freeNBodyTree(NBodyTree* t)
{
    NBodyNode* p;
    NBodyNode* tmp;

    p = (NBodyNode*) t->root;
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

static void freeFreeCells(NBodyNode* freecell)
{
    NBodyNode* p;
    NBodyNode* tmp;

    p = freecell;
    while (p)
    {
        tmp = Next(p);
        free(p);
        p = tmp;
    }
}

int detachSharedScene(NBodyState* st)
{
 #if USE_SHMEM
    if (st->scene)
    {
        if (shmdt(st->scene) < 0)
        {
            perror("Closing shared memory");
            return 1;
        }

        st->shmId = -1;
        st->scene = NULL;
    }
  #endif /* _WIN32 */

    return 0;
}

int destroyNBodyState(NBodyState* st)
{
    int failed = FALSE;

    freeNBodyTree(&st->tree);
    freeFreeCells(st->freecell);
    mwFreeA(st->bodytab);
    mwFreeA(st->acctab);

    free(st->checkpointResolved);

    if (st->outFile && st->outFile != DEFAULT_OUTPUT_FILE)
    {
        if (fclose(st->outFile))
        {
            perror("closing output\n");
            failed = TRUE;
        }
    }

    failed |= detachSharedScene(st);

    return failed;
}

void setInitialNBodyState(NBodyState* st, const NBodyCtx* ctx, Body* bodies, unsigned int nbody)
{
    static const NBodyTree emptyTree = EMPTY_TREE;

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
    return mwCallocA(1, sizeof(NBodyState));
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
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

/* TODO: Doesn't handle tree */
/* Returns nonzero if states are equal, 0 otherwise */
int equalNBodyState(const NBodyState* st1, const NBodyState* st2)
{
    assert(st1 && st2);

    if (   st1->tnow != st2->tnow
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
void cloneNBodyState(NBodyState* st, const NBodyState* oldSt)
{
    static const NBodyTree emptyTree = EMPTY_TREE;
    unsigned int nbody = oldSt->nbody;

    st->tree = emptyTree;
    st->tree.rsize = oldSt->tree.rsize;

    st->freecell = NULL;

    st->lastCheckpoint = oldSt->lastCheckpoint;
    st->tnow           = oldSt->tnow;
    st->nbody          = oldSt->nbody;

    st->treeIncest = oldSt->treeIncest;
    st->tree.structureError = oldSt->tree.structureError;

    assert(nbody > 0);
    assert(st->bodytab == NULL && st->acctab == NULL);

    st->bodytab = (Body*) mwMallocA(sizeof(Body) * nbody);
    memcpy(st->bodytab, oldSt->bodytab, sizeof(Body) * nbody);

    st->acctab = (mwvector*) mwMallocA(sizeof(mwvector) * nbody);
    memcpy(st->acctab, oldSt->acctab, sizeof(mwvector) * nbody);
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

/* Floating point comparison where nan compares equal */
static int feqWithNan(real a, real b)
{
    return (isnan(a) && isnan(b)) ? TRUE : (a == b);
}


int equalDisk(const Disk* d1, const Disk* d2)
{
    return (d1->type == d2->type)
        && feqWithNan(d1->mass,d2->mass)
        && feqWithNan(d1->scaleLength, d1->scaleLength)
        && feqWithNan(d1->scaleHeight, d1->scaleHeight);
}

int equalHalo(const Halo* h1, const Halo* h2)
{
    return (h1->type == h2->type)
        && feqWithNan(h1->vhalo, h2->vhalo)
        && feqWithNan(h1->scaleLength, h2->scaleLength)
        && feqWithNan(h1->flattenZ, h2->flattenZ)
        && feqWithNan(h1->flattenY, h2->flattenY)
        && feqWithNan(h1->flattenX, h2->flattenX)
        && feqWithNan(h1->triaxAngle, h2->triaxAngle)
        && feqWithNan(h1->c1, h2->c1)
        && feqWithNan(h1->c2, h2->c2)
        && feqWithNan(h1->c3, h2->c3);
}

int equalSpherical(const Spherical* s1, const Spherical* s2)
{
    return (s1->type == s2->type)
        && feqWithNan(s1->mass, s2->mass)
        && feqWithNan(s1->scale, s2->scale);
}

int equalPotential(const Potential* p1, const Potential* p2)
{
    return equalSpherical(&p1->sphere[0], &p2->sphere[0])
        && equalDisk(&p1->disk, &p2->disk)
        && equalHalo(&p1->halo, &p2->halo);
}

int equalHistogramParams(const HistogramParams* hp1, const HistogramParams* hp2)
{
    return feqWithNan(hp1->phi, hp2->phi)
        && feqWithNan(hp1->theta, hp2->theta)
        && feqWithNan(hp1->psi, hp2->psi)
        && feqWithNan(hp1->startRaw, hp2->startRaw)
        && feqWithNan(hp1->endRaw, hp2->endRaw)
        && feqWithNan(hp1->binSize, hp2->binSize)
        && feqWithNan(hp1->center, hp2->center);
}

int equalNBodyCtx(const NBodyCtx* ctx1, const NBodyCtx* ctx2)
{
    return (ctx1->potentialType == ctx2->potentialType)
        && feqWithNan(ctx1->timestep, ctx2->timestep)
        && feqWithNan(ctx1->timeEvolve, ctx2->timeEvolve)
        && feqWithNan(ctx1->theta, ctx2->theta)
        && feqWithNan(ctx1->eps2, ctx2->eps2)
        && feqWithNan(ctx1->treeRSize, ctx2->treeRSize)
        && feqWithNan(ctx1->sunGCDist, ctx2->sunGCDist)
        && feqWithNan(ctx1->criterion, ctx2->criterion)
        && feqWithNan(ctx1->useQuad, ctx2->useQuad)
        && feqWithNan(ctx1->allowIncest, ctx2->allowIncest)
        && feqWithNan(ctx1->quietErrors, ctx2->quietErrors)
        && feqWithNan(ctx1->checkpointT, ctx2->checkpointT)
        && feqWithNan(ctx1->freqOut, ctx2->freqOut)
        && equalHistogramParams(&ctx1->histogramParams, &ctx2->histogramParams)
        && equalPotential(&ctx1->pot, &ctx2->pot);
}

