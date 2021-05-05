/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
   Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
   Copyright (c) 2016-2018 Siddhartha Shelton
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkyway@Home.

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

#if USE_POSIX_SHMEM
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
  #if USE_POSIX_SHMEM
    if (st->scene)
    {
        if (shm_unlink(st->scene->shmemName) < 0)
        {
            mwPerror("Closing shared scene memory '%s'", st->scene->shmemName);
            return 1;
        }
    }
  #elif USE_WIN32_SHARED_MAP
    if (st->scene)
    {
        if (!UnmapViewOfFile((LPCVOID) st->scene))
        {
            mwPerrorW32("Error unmapping shared scene memory");
            return 1;
        }
    }
  #endif /* USE_POSIX_SHMEM */

    st->scene = NULL;
    return 0;
}


int destroyNBodyState(NBodyState* st)
{
    int failed = FALSE;
    int nThread = nbGetMaxThreads();
    int i;
    //mw_printf("After Initial\n");

    freeNBodyTree(&st->tree);
    //mw_printf("After Free Tree\n");
    freeFreeCells(st->freeCell);
    //mw_printf("After Free Cells\n");
    mwFreeA(st->bodytab);
    //mw_printf("After Free bodytab\n");
    mwFreeA(st->bestLikelihoodBodyTab);
    //mw_printf("After Free bestLikeTabl\n");
    mwFreeA(st->acctab);
    //mw_printf("After Free acctab\n");
    mwFreeA(st->orbitTrace);
    //mw_printf("After Free orbitTrace\n");
    
    if(st->shiftByLMC) {   
        mwFreeA(st->shiftByLMC);
    }
    //mw_printf("After Free LMCShift\n");
    
    free(st->checkpointResolved);
    //mw_printf("After Free checkpointResolved\n");

    if (st->potEvalStates)
    {
        for (i = 0; i < nThread; ++i)
        {
            lua_close(st->potEvalStates[i]);
        }
        free(st->potEvalClosures);
        free(st->potEvalStates);
    }
    //mw_printf("After Free potEvalStates\n");

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
    //mw_printf("After OpenCL\n");


    failed |= nbDetachSharedScene(st);
    //mw_printf("After nbDetachSharedScene\n");

    return failed;
}

void setInitialNBodyState(NBodyState* st, const NBodyCtx* ctx, Body* bodies, int nbody)
{
    static const NBodyTree emptyTree = EMPTY_TREE;

    st->tree = emptyTree;
    st->freeCell = NULL;
    st->usesQuad = ctx->useQuad;
    st->usesExact = (ctx->criterion == Exact);

    st->tree.rsize = ctx->treeRSize;
    st->step = 0;
    st->nbody = nbody;
    st->bodytab = bodies;
    st->bestLikelihoodBodyTab = (Body*) mwMallocA(nbody * sizeof(Body));
    memcpy(st->bestLikelihoodBodyTab, st->bodytab, nbody * sizeof(Body));
    st->bestLikelihood         = DEFAULT_WORST_CASE;
    st->bestLikelihood_EMD     = DEFAULT_WORST_CASE;
    st->bestLikelihood_Mass    = DEFAULT_WORST_CASE;
    st->bestLikelihood_Beta    = DEFAULT_WORST_CASE;
    st->bestLikelihood_Vel     = DEFAULT_WORST_CASE;
    st->bestLikelihood_BetaAvg = DEFAULT_WORST_CASE;
    st->bestLikelihood_VelAvg  = DEFAULT_WORST_CASE;
    st->bestLikelihood_Dist    = DEFAULT_WORST_CASE;
    st->bestLikelihood_time    = 0.0;
    st->bestLikelihood_count   = 0;
    
    /* We'll report the center of mass for each step + the initial one */
    st->nOrbitTrace = ctx->nStep + 1;
    st->orbitTrace = (mwvector*) mwCallocA(st->nOrbitTrace, sizeof(mwvector));

    /* The tests may step the system from an arbitrary place, so make sure this is 0'ed */
    st->acctab = (mwvector*) mwCallocA(nbody, sizeof(mwvector));

}

void setRandomLMCNBodyState(NBodyState* st, int nShift, dsfmt_t* dsfmtState)
{
    int j;

    st->shiftByLMC = (mwvector*)mwCallocA(nShift, sizeof(mwvector));
    for(j = 0; j < nShift; j++) {
        st->shiftByLMC[j] = mwRandomVector(dsfmtState, mwXrandom(dsfmtState,0.0,1.0));
        //SET_VECTOR(st->shiftByLMC[j],0.0,0.0,0.0);
    }

    st->LMCpos = mwRandomVector(dsfmtState, mwXrandom(dsfmtState,0.01,200.0));
    st->LMCvel = mwRandomVector(dsfmtState, mwXrandom(dsfmtState,0.01,200.0));
    //SET_VECTOR(*(st->LMCpos),0.0,0.0,0.0);
    //SET_VECTOR(*(st->LMCvel),0.0,0.0,0.0);
    st->nShiftLMC = nShift;
}

void setLMCShiftArray(NBodyState* st, mwvector* shiftArray, size_t shiftSize) {
    //Set the state variable for the LMC shift array
    st->shiftByLMC = shiftArray;
    st->nShiftLMC = shiftSize;
}

void setLMCPosVel(NBodyState* st, mwvector Pos, mwvector Vel) {
    //Set the state variable for the LMC position and velocity
    st->LMCpos = Pos;
    st->LMCvel = Vel;
}

NBodyState* newNBodyState()
{
    return mwCallocA(1, sizeof(NBodyState));
}

#if NBODY_OPENCL

NBodyStatus nbInitCL(NBodyState* st, const NBodyCtx* ctx, const CLRequest* clr)
{
    cl_int err;

    st->usesQuad = ctx->useQuad;
    st->usesExact = (ctx->criterion == Exact);
    st->usesCL = TRUE;
    st->useCLCheckpointing = clr->enableCheckpointing;

    st->ci = mwCalloc(1, sizeof(CLInfo));
    st->nbb = mwCalloc(1, sizeof(NBodyBuffers));
    st->workSizes = mwCalloc(1, sizeof(NBodyWorkSizes));
    st->kernels = mwCalloc(1, sizeof(NBodyKernels));

    err = mwSetupCL(st->ci, clr);
    if (err != CL_SUCCESS)
        return NBODY_CL_ERROR;

    return NBODY_SUCCESS;
}

NBodyStatus nbInitNBodyStateCL(NBodyState* st, const NBodyCtx* ctx)
{
    cl_int err;
    const DevInfo* devInfo;

    if (!st->usesCL)
    {
        mw_printf("CL not setup for CL state initialization\n");
        return NBODY_CONSISTENCY_ERROR;
    }

    /* Bodies must be set before trying to use this */
    if (!st->bodytab)
    {
        mw_printf("Bodies not set for CL state initialization\n");
        return NBODY_CONSISTENCY_ERROR;
    }

    if (ctx->potentialType == EXTERNAL_POTENTIAL_CUSTOM_LUA)
    {
        mw_printf("Cannot use Lua potential with OpenCL\n");
        return NBODY_UNSUPPORTED;
    }

    devInfo = &st->ci->di;

    if (!nbCheckDevCapabilities(devInfo, ctx, st->nbody))
        return NBODY_CAPABILITY_ERROR;

    if (   nbSetThreadCounts(st->workSizes, devInfo, ctx)
        || nbSetWorkSizes(st->workSizes, devInfo, st->nbody, st->ignoreResponsive))
        return NBODY_ERROR;

    st->effNBody = nbFindEffectiveNBody(st->workSizes, st->usesExact, st->nbody);
    st->maxDepth = nbFindMaxDepthForDevice(devInfo, st->workSizes, ctx->useQuad);

    st->usesConsistentMemory =  (mwIsNvidiaGPUDevice(devInfo) && mwNvidiaInlinePTXAvailable(st->ci->plat))
                              || mwDeviceHasConsistentMemory(devInfo);

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
    {
        mw_printf("Error marshalling initial bodies\n");
        return NBODY_CL_ERROR;
    }

    return NBODY_SUCCESS;
}

#else

NBodyStatus nbInitCL(NBodyState* st, const NBodyCtx* ctx, const CLRequest* clr)
{
    (void) st, (void) ctx, (void) clr;
    return NBODY_CL_ERROR;
}

NBodyStatus nbInitNBodyStateCL(NBodyState* st, const NBodyCtx* ctx)
{
    (void) st, (void) ctx;
    return NBODY_CL_ERROR;
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
        mw_printf("Position difference detected!\n   Difference = [%.15f,%.15f,%.15f]\n", X(Pos(a))-X(Pos(b)),Y(Pos(a))-Y(Pos(b)),Z(Pos(a))-Z(Pos(b)));
        return FALSE;
    }
    if (!equalVector(&Vel(a), &Vel(b)))
    {
        mw_printf("Velocity difference detected!\n   Difference = [%.15f,%.15f,%.15f]\n", X(Vel(a))-X(Vel(b)),Y(Vel(a))-Y(Vel(b)),Z(Vel(a))-Z(Vel(b)));
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
            mw_printf("   Difference = [%.15f,%.15f,%.15f]\n", X(a[i])-X(b[i]),Y(a[i])-Y(b[i]),Z(a[i])-Z(b[i]));
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
        //mw_printf("Body Number: %d\n",i+1);
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

    if (st1->shiftByLMC || st2->shiftByLMC)
    {
        if (!st1->shiftByLMC || !st2->shiftByLMC)
        {
            mw_printf("Comparing non-NULL shiftByLMC to NULL pointer!\n");
            return FALSE;
        }
        if (st1->nShiftLMC != st2->nShiftLMC)
        {
            mw_printf("shiftByLMC Size Difference!\n");
            return FALSE;
        }
        if (!equalVectorArray(st1->shiftByLMC, st2->shiftByLMC, st1->nShiftLMC))
        {
            mw_printf("Different LMC Shifts Detected!\n");
            return FALSE;
        }
    }

    if (!equalVector(&st1->LMCpos, &st2->LMCpos))
    {
        mw_printf("LMC Position difference detected!\n   Difference = [%.15f,%.15f,%.15f]\n", X(st1->LMCpos)-X(st2->LMCpos),Y(st1->LMCpos)-Y(st2->LMCpos),Z(st1->LMCpos)-Z(st2->LMCpos));
        return FALSE;
    }

    if (!equalVector(&st1->LMCvel, &st2->LMCvel))
    {
        mw_printf("LMC Velocity difference detected!\n   Difference = [%.15f,%.15f,%.15f]\n", X(st1->LMCvel)-X(st2->LMCvel),Y(st1->LMCvel)-Y(st2->LMCvel),Z(st1->LMCvel)-Z(st2->LMCvel));
        return FALSE;
    }

    if (st1->orbitTrace || st2->orbitTrace)
    {
        if (!st1->orbitTrace || !st2->orbitTrace)
        {
            mw_printf("Comparing non-NULL orbitTrace to NULL pointer!\n");
            return FALSE;
        }
        if (st1->nOrbitTrace != st2->nOrbitTrace)
        {
            mw_printf("orbitTrace Size Difference!\n");
            return FALSE;
        }
        if (!equalVectorArray(st1->orbitTrace, st2->orbitTrace, st1->nOrbitTrace))
        {
            mw_printf("Different Orbits Detected!\n");
            return FALSE;
        }
    }

    if (!equalBodyArray(st1->bodytab, st2->bodytab, st1->nbody))
    {
        return FALSE;
    }

    if (!equalBodyArray(st1->bestLikelihoodBodyTab, st2->bestLikelihoodBodyTab, st1->nbody))
    {
        return FALSE;
    }

    if (!equalVectorArray(st1->acctab, st2->acctab, st1->nbody))
    {
        mw_printf("Different Accelerations Detected!\n");
        return FALSE;
    }

    return TRUE;
}

/* TODO: Doesn't clone tree or CL stuffs */
void cloneNBodyState(NBodyState* st, const NBodyState* oldSt)
{
    static const NBodyTree emptyTree = EMPTY_TREE;
    unsigned int nbody = oldSt->nbody;
    st->tree = emptyTree;
    st->tree.rsize = oldSt->tree.rsize;

    st->freeCell = NULL;

    st->lastCheckpoint       = oldSt->lastCheckpoint;
    st->step                 = oldSt->step;
    st->nbody                = oldSt->nbody;
    st->effNBody             = oldSt->effNBody;
    st->bestLikelihood       = oldSt->bestLikelihood;
    st->bestLikelihood_count = oldSt->bestLikelihood_count;
    
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
    
    st->bestLikelihoodBodyTab = (Body*) mwMallocA(nbody * sizeof(Body));
    memcpy(st->bestLikelihoodBodyTab, oldSt->bestLikelihoodBodyTab, nbody * sizeof(Body));


    st->acctab = (mwvector*) mwMallocA(nbody * sizeof(mwvector));
    memcpy(st->acctab, oldSt->acctab, nbody * sizeof(mwvector));

    st->previousForwardTime = oldSt->previousForwardTime;

    if (oldSt->orbitTrace)
    {
        st->orbitTrace = (mwvector*) mwMallocA(oldSt->nOrbitTrace * sizeof(mwvector));
        memcpy(st->orbitTrace, oldSt->orbitTrace, oldSt->nOrbitTrace * sizeof(mwvector));
        st->nOrbitTrace = oldSt->nOrbitTrace;
    }

    if (oldSt->shiftByLMC)
    {
        st->shiftByLMC = (mwvector*) mwMallocA(oldSt->nShiftLMC * sizeof(mwvector));
        memcpy(st->shiftByLMC, oldSt->shiftByLMC, oldSt->nShiftLMC * sizeof(mwvector));
        st->nShiftLMC = oldSt->nShiftLMC;
    }

    st->LMCpos = oldSt->LMCpos;
    st->LMCvel = oldSt->LMCvel;

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
        && feqWithNan(h1->mass, h2->mass)
        && feqWithNan(h1->scaleLength, h2->scaleLength)
        && feqWithNan(h1->flattenZ, h2->flattenZ)
        && feqWithNan(h1->flattenY, h2->flattenY)
        && feqWithNan(h1->flattenX, h2->flattenX)
        && feqWithNan(h1->triaxAngle, h2->triaxAngle)
        && feqWithNan(h1->c1, h2->c1)
        && feqWithNan(h1->c2, h2->c2)
        && feqWithNan(h1->c3, h2->c3)
        && feqWithNan(h1->rho0, h2->rho0);
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
        && equalDisk(&p1->disk2, &p2->disk2)
        && equalHalo(&p1->halo, &p2->halo);
}

int equalHistogramParams(const HistogramParams* hp1, const HistogramParams* hp2)
{
    return feqWithNan(hp1->phi, hp2->phi)
        && feqWithNan(hp1->theta, hp2->theta)
        && feqWithNan(hp1->psi, hp2->psi)
        && feqWithNan(hp1->lambdaStart, hp2->lambdaStart)
        && feqWithNan(hp1->lambdaEnd, hp2->lambdaEnd)
        && feqWithNan(hp1->lambdaBins, hp2->lambdaBins)
        && feqWithNan(hp1->betaStart, hp2->betaStart)
        && feqWithNan(hp1->betaEnd, hp2->betaEnd)
        && feqWithNan(hp1->betaBins, hp2->betaBins);
}

int equalNBodyCtx(const NBodyCtx* ctx1, const NBodyCtx* ctx2)
{
    return feqWithNan(ctx1->eps2, ctx2->eps2)
        && feqWithNan(ctx1->theta, ctx2->theta)
        && feqWithNan(ctx1->timestep, ctx2->timestep)
        && feqWithNan(ctx1->timeEvolve, ctx2->timeEvolve)
        && feqWithNan(ctx1->timeBack, ctx2->timeBack)
        && feqWithNan(ctx1->treeRSize, ctx2->treeRSize)
        && feqWithNan(ctx1->sunGCDist, ctx2->sunGCDist)
        && feqWithNan(ctx1->criterion, ctx2->criterion)
        && (ctx1->potentialType == ctx2->potentialType)
        && feqWithNan(ctx1->useQuad, ctx2->useQuad)
        && feqWithNan(ctx1->allowIncest, ctx2->allowIncest)
        && feqWithNan(ctx1->useBestLike, ctx2->useBestLike)
        && feqWithNan(ctx1->useVelDisp, ctx2->useVelDisp)
        && feqWithNan(ctx1->useBetaDisp, ctx2->useBetaDisp)
        && feqWithNan(ctx1->useBetaComp, ctx2->useBetaComp)
        && feqWithNan(ctx1->useVlos, ctx2->useVlos)
        && feqWithNan(ctx1->useDist, ctx2->useDist)
        && feqWithNan(ctx1->BestLikeStart, ctx2->BestLikeStart)
        && feqWithNan(ctx1->Nstep_control, ctx2->Nstep_control)
        && feqWithNan(ctx1->Ntsteps, ctx2->Ntsteps)
        && feqWithNan(ctx1->MultiOutput, ctx2->MultiOutput)
        && feqWithNan(ctx1->OutputFreq, ctx2->OutputFreq)
        && feqWithNan(ctx1->BetaSigma, ctx2->BetaSigma)
        && feqWithNan(ctx1->VelSigma, ctx2->VelSigma)
        && feqWithNan(ctx1->DistSigma, ctx2->DistSigma)
        && feqWithNan(ctx1->IterMax, ctx2->IterMax)
        && feqWithNan(ctx1->BetaCorrect, ctx2->BetaCorrect)
        && feqWithNan(ctx1->VelCorrect, ctx2->VelCorrect)
        && feqWithNan(ctx1->DistCorrect, ctx2->DistCorrect)
        && feqWithNan(ctx1->quietErrors, ctx2->quietErrors)
        && ctx1->checkpointT == ctx2->checkpointT
        && feqWithNan(ctx1->nStep, ctx2->nStep)
        && equalPotential(&ctx1->pot, &ctx2->pot)
        && feqWithNan(ctx1->LMC, ctx2->LMC)
        && feqWithNan(ctx1->LMCmass, ctx2->LMCmass)
        && feqWithNan(ctx1->LMCscale, ctx2->LMCscale)
        && feqWithNan(ctx1->LMCDynaFric, ctx2->LMCDynaFric)
        && feqWithNan(ctx1->calibrationRuns, ctx2->calibrationRuns);
}

