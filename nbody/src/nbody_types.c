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
#include "nbody_util.h"
#include "nbody_types.h"
#include "nbody_show.h"
#include "nbody_defaults.h"

#if NBODY_OPENCL
  #include "nbody_cl.h"
#endif /* NBODY_OPENCL */

#if USE_SHMEM
  #include <sys/mman.h>
#endif

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
            mwFreeA(p);
            p = tmp;
        }
        else                        /* skip over bodies */
        {
            p = Next(p);
        }
    }

    t->root = NULL;
    t->cellUsed = 0;
    t->maxDepth = 0;
}

static void freeFreeCells(NBodyNode* freeCell)
{
    NBodyNode* p;
    NBodyNode* tmp;

    p = freeCell;
    while (p)
    {
        tmp = Next(p);
        mwFreeA(p);
        p = tmp;
    }
}

int nbDetachSharedScene(NBodyState* st)
{
    (void) st;

  #if USE_SHMEM
    if (st->scene)
    {
        if (shm_unlink(st->scene->shmemName) < 0)
        {
            mwPerror("Closing shared scene memory '%s'", st->scene->shmemName);
            return 1;
        }

        st->shmId = -1;
        st->scene = NULL;
    }
  #endif /* USE_SHMEM */

    return 0;
}



int destroyNBodyState(NBodyState* st)
{
    int failed = FALSE;
    int nThread = nbGetMaxThreads();
    int i;

    freeNBodyTree(&st->tree);
    freeFreeCells(st->freeCell);
    mwFreeA(st->bodytab);
    mwFreeA(st->acctab);
    mwFreeA(st->orbitTrace);

    free(st->checkpointResolved);

    if (st->potEvalStates)
    {
        for (i = 0; i < nThread; ++i)
        {
            lua_close(st->potEvalStates[i]);
        }
        free(st->potEvalClosures);
        free(st->potEvalStates);
    }

  #if NBODY_OPENCL

    if (st->ci)
    {
        mwDestroyCLInfo(st->ci);
        free(st->ci);
        st->ci = NULL;
    }

    if (st->kernels)
    {
        cl_int err;

        err = nbReleaseKernels(st);
        free(st->kernels);
        failed |= (err != CL_SUCCESS);
    }

    if (st->workSizes)
    {
        free(st->workSizes);
    }

    if (st->nbb)
    {
        cl_int err;

        err = nbReleaseBuffers(st);
        free(st->nbb);
        st->nbb = NULL;
        failed |= (err != CL_SUCCESS);
    }

  #endif /* NBODY_OPENCL */


    failed |= nbDetachSharedScene(st);

    return failed;
}

void setInitialNBodyState(NBodyState* st, const NBodyCtx* ctx, Body* bodies, int nbody)
{
    static const NBodyTree emptyTree = EMPTY_TREE;
    static const mwvector maxV = mw_vec(REAL_MAX, REAL_MAX, REAL_MAX);
    int i;

    st->tree = emptyTree;
    st->freeCell = NULL;
    st->usesQuad = ctx->useQuad;
    st->usesExact = (ctx->criterion == Exact);

    st->tree.rsize = ctx->treeRSize;
    st->step = 0;
    st->nbody = nbody;
    st->bodytab = bodies;


    st->orbitTrace = (mwvector*) mwMallocA(N_ORBIT_TRACE_POINTS * sizeof(mwvector));
    for (i = 0; i < N_ORBIT_TRACE_POINTS; ++i)
    {
        st->orbitTrace[i] = maxV;
    }

    /* The tests may step the system from an arbitrary place, so make sure this is 0'ed */
    st->acctab = (mwvector*) mwCallocA(nbody, sizeof(mwvector));
}

NBodyState* newNBodyState()
{
    return mwCallocA(1, sizeof(NBodyState));
}

#if NBODY_OPENCL

NBodyStatus nbInitNBodyStateCL(NBodyState* st, const NBodyCtx* ctx, const CLRequest* clr)
{
    cl_int err;

    /* Bodies must be set before trying to use this */
    if (!st->bodytab)
    {
        return NBODY_CONSISTENCY_ERROR;
    }

    if (ctx->potentialType == EXTERNAL_POTENTIAL_CUSTOM_LUA)
    {
        mw_printf("Cannot use Lua potential with OpenCL\n");
        return NBODY_UNSUPPORTED;
    }

    st->usesQuad = ctx->useQuad;
    st->usesExact = (ctx->criterion == Exact);
    st->usesCL = TRUE;

    st->ci = mwCalloc(1, sizeof(CLInfo));
    st->nbb = mwCalloc(1, sizeof(NBodyBuffers));
    st->workSizes = mwCalloc(1, sizeof(NBodyWorkSizes));
    st->kernels = mwCalloc(1, sizeof(NBodyKernels));

    err = mwSetupCL(st->ci, clr);
    if (err != CL_SUCCESS)
        return NBODY_CL_ERROR;

    if (!nbCheckDevCapabilities(&st->ci->di, ctx, st->nbody))
        return NBODY_CAPABILITY_ERROR;

    if (nbSetThreadCounts(st->workSizes, &st->ci->di, ctx) || nbSetWorkSizes(st->workSizes, &st->ci->di))
        return NBODY_ERROR;

    st->effNBody = nbFindEffectiveNBody(st->workSizes, st->usesExact, st->nbody);

    if (nbLoadKernels(ctx, st))
        return NBODY_CL_ERROR;

    err = nbCreateBuffers(ctx, st);
    if (err != CL_SUCCESS)
        return NBODY_CL_ERROR;

    err = nbSetInitialTreeStatus(st);
    if (err != CL_SUCCESS)
        return NBODY_CL_ERROR;

    err = nbSetAllKernelArguments(st);
    if (err != CL_SUCCESS)
        return NBODY_CL_ERROR;

    err = nbMarshalBodies(st, CL_TRUE);
    if (err != CL_SUCCESS)
        return NBODY_CL_ERROR;

    return NBODY_SUCCESS;
}

#endif /* NBODY_OPENCL */


#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

static int equalVector(const mwvector* a, const mwvector* b)
{
    return (a->x == b->x && a->y == b->y && a->z == b->z);
}

int equalBody(const Body* a, const Body* b)
{
    if (Mass(a) != Mass(b))
    {
        mw_printf("mass differ\n");
        return FALSE;
    }
    if (Type(a) != Type(b))
    {
        mw_printf("type ndiffer\n");
        return FALSE;
    }
    if (!equalVector(&Pos(a), &Pos(b)))
    {
        mw_printf("pos differ\n");
        return FALSE;
    }
    if (!equalVector(&Vel(a), &Vel(b)))
    {
        mw_printf("VElocity differ\n");
        return FALSE;
    }

    return TRUE;
}

static int equalVectorArray(const mwvector* a, const mwvector* b, size_t n)
{
    size_t i;

    if (!a && !b)
        return TRUE;

    if (!a || !b)
        return FALSE;

    for (i = 0; i < n; ++i)
    {
        if (!equalVector(&a[i], &b[i]))
        {
            return FALSE;
        }
    }

    return TRUE;
}

static int equalBodyArray(const Body* a, const Body* b, size_t n)
{
    size_t i;

    if (!a && !b)
        return TRUE;

    if (!a || !b)
        return FALSE;

    for (i = 0; i < n; ++i)
    {
        if (!equalBody(&a[i], &b[i]))
        {
            return FALSE;
        }
    }

    return TRUE;
}

/* TODO: Doesn't handle tree or other parts */
/* Returns nonzero if states are equal, 0 otherwise */
int equalNBodyState(const NBodyState* st1, const NBodyState* st2)
{
    assert(st1 && st2);

    if (   st1->step != st2->step
        || st1->lastCheckpoint != st2->lastCheckpoint
        || st1->nbody != st2->nbody)
    {
        return FALSE;
    }

    if (!equalBodyArray(st1->bodytab, st2->bodytab, st1->nbody))
    {
        return FALSE;
    }

    return equalVectorArray(st1->acctab, st2->acctab, st1->nbody);
}

/* TODO: Doesn't clone tree or CL stuffs */
void cloneNBodyState(NBodyState* st, const NBodyState* oldSt)
{
    static const NBodyTree emptyTree = EMPTY_TREE;
    unsigned int nbody = oldSt->nbody;

    st->tree = emptyTree;
    st->tree.rsize = oldSt->tree.rsize;

    st->freeCell = NULL;

    st->lastCheckpoint = oldSt->lastCheckpoint;
    st->step           = oldSt->step;
    st->nbody          = oldSt->nbody;
    st->effNBody       = oldSt->effNBody;

    st->ignoreResponsive = oldSt->ignoreResponsive;
    st->usesExact = oldSt->usesExact;
    st->usesQuad = oldSt->usesQuad,
    st->dirty = oldSt->dirty;
    st->usesCL = oldSt->usesCL;
    st->reportProgress = oldSt->reportProgress;

    st->treeIncest = oldSt->treeIncest;
    st->tree.structureError = oldSt->tree.structureError;

    assert(nbody > 0);
    assert(st->bodytab == NULL && st->acctab == NULL);

    st->bodytab = (Body*) mwMallocA(nbody * sizeof(Body));
    memcpy(st->bodytab, oldSt->bodytab, nbody * sizeof(Body));

    st->acctab = (mwvector*) mwMallocA(nbody * sizeof(mwvector));
    memcpy(st->acctab, oldSt->acctab, nbody * sizeof(mwvector));

    st->orbitTrace = (mwvector*) mwMallocA(N_ORBIT_TRACE_POINTS * sizeof(mwvector));
    memcpy(st->orbitTrace, oldSt->orbitTrace, N_ORBIT_TRACE_POINTS * sizeof(mwvector));

    if (st->ci)
    {
        mw_panic("OpenCL NBodyState cloning not implemented\n");

        /*
        st->ci = (CLInfo*) mwCalloc(1, sizeof(CLInfo));
        st->nbb = (NBodyBuffers*) mwCalloc(1, sizeof(NBodyBuffers));

        memcpy(st->ci, oldSt->ci, sizeof(CLInfo));

        clRetainContext(oldSt->ci->clctx);
        clRetainProgram(oldSt->ci->prog);
        clRetainCommandQueue(oldSt->ci->queue);

        // copy buffers
        mwDuplicateBuffer(st->ci, oldSt->nbb.blah)
        */
    }

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
void sortBodies(Body* bodies, int nbody)
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
        && feqWithNan(d1->mass, d2->mass)
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
    return feqWithNan(ctx1->eps2, ctx2->eps2)
        && feqWithNan(ctx1->theta, ctx2->theta)
        && feqWithNan(ctx1->timestep, ctx2->timestep)
        && feqWithNan(ctx1->timeEvolve, ctx2->timeEvolve)
        && feqWithNan(ctx1->treeRSize, ctx2->treeRSize)
        && feqWithNan(ctx1->sunGCDist, ctx2->sunGCDist)
        && feqWithNan(ctx1->criterion, ctx2->criterion)
        && (ctx1->potentialType == ctx2->potentialType)
        && feqWithNan(ctx1->useQuad, ctx2->useQuad)
        && feqWithNan(ctx1->allowIncest, ctx2->allowIncest)
        && feqWithNan(ctx1->quietErrors, ctx2->quietErrors)
        && feqWithNan(ctx1->checkpointT, ctx2->checkpointT)
        && feqWithNan(ctx1->nStep, ctx2->nStep)
        && equalPotential(&ctx1->pot, &ctx2->pot);
}

