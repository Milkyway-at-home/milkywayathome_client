/*
  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
  Copyright (c) 2010 Matthew Arsenault, Travis Desell, Boleslaw
    Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
    Rensselaer Polytechnic Institute.
  Copyright (c) 2002-2006 John M. Fregeau, Richard Campbell, Jeff Molofee

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

#include "nbody.h"
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "nbody_grav.h"
#include "nbody_show.h"
#include "nbody_lua.h"
#include "milkyway_cpp_util.h"

#if USE_SHMEM
  #include <sys/shm.h>
#endif

#if USE_SHMEM

static void* createSharedMemory(key_t key, size_t size, int* shmIdOut)
{
    int shmId;
    void* p;

    shmId = shmget(key, size, IPC_CREAT | SHM_W | SHM_R);
    if (shmId < 0)
    {
        perror("Error getting shared memory");
        return NULL;
    }

    p = shmat(shmId, NULL, 0);
    if (!p || p == (void*) -1)
    {
        perror("Getting shared memory");
        return NULL;
    }

    memset(p, 0, size);

    if (shmIdOut)
        *shmIdOut = shmId;

    return p;
}

static int createSharedScene(NBodyState* st, const char* inputFile)
{
    key_t key;
    size_t size;
    int shmId = -1;

    size = sizeof(scene_t) + st->nbody * sizeof(FloatPos);
    key = DEFAULT_SHMEM_KEY;
    //key = ftok(inputFile, getpid());
    if (key < 0)
    {
        key = DEFAULT_SHMEM_KEY;
        perror("Error getting key");
        return 1;
    }

    st->scene = (scene_t*) createSharedMemory(key, size, &shmId);
    if (!st->scene)
        return 1;

    st->shmId = shmId;
    st->scene->nbody = st->nbody;

    return 0;
}

static int visualizerIsAttached(const NBodyState* st)
{
    struct shmid_ds buf;

    if (shmctl(st->shmId, IPC_STAT, &buf) < 0)
    {
        perror("Finding attached visualizers");
        return 0;
    }

    return buf.shm_nattch > 1;
}

#elif USE_BOINC_SHMEM

static int createSharedScene(NBodyState* st, const char* inputFile)
{
    size_t size = sizeof(scene_t) + st->nbody * sizeof(FloatPos);

    st->scene = (scene_t*) boinc_graphics_make_shmem("milkyway_nbody", (int) size);
    if (!st->scene)
    {
        warn("Failed to get shmem of size %d\n", size);
        return 1;
    }

    memset(scene, 0, size);
    st->scene->nbody = st->nbody;

    return 0;
}

static int visualizerIsAttached(const NBodyState* st)
{
    /* TODO */
    return 1;
}

#else

static int visualizerIsAttached(const NBodyState* st)
{
    return 0;
}

#endif /* USE_BOINC_SHMEM */


static void updateDisplayedBodies(NBodyState* st)
{
    const Body* b;
    FloatPos* r;
    unsigned int i = 0;
    const unsigned int nbody = st->nbody;
    scene_t* scene = st->scene;

    if (!scene)
        return;

    r = scene->r;
    scene->usleepcount += scene->usleepdt;

    /* read data if not paused */
    if (scene->usleepcount >= scene->dt && (!scene->paused || scene->step == 1))
    {
        scene->usleepcount = 0.0;
        scene->step = 0;

      #ifdef _OPENMP
        #pragma omp parallel for private(i, b) schedule(static)
      #endif
        for (i = 0; i < nbody; ++i)
        {
            b = &st->bodytab[i];
            r[i].x = (float) X(Pos(b));
            r[i].y = (float) Y(Pos(b));
            r[i].z = (float) Z(Pos(b));
        }
    }

    scene->changed = TRUE;
}

static inline int nbodyTimeToCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
  #if BOINC_APPLICATION
    return boinc_time_to_checkpoint();
  #else
    time_t now;

    if (ctx->checkpointT < 0 || ((now = time(NULL)) - st->lastCheckpoint) < ctx->checkpointT)
        return FALSE;
    else
    {
        st->lastCheckpoint = now;
        return TRUE;
    }
  #endif /* BOINC_APPLICATION */
}

static inline void nbodyCheckpoint(const NBodyCtx* ctx, NBodyState* st)
{
    if (nbodyTimeToCheckpoint(ctx, st))
    {
        if (writeCheckpoint(ctx, st))
            fail("Failed to write checkpoint\n");

      #if BOINC_APPLICATION
        boinc_checkpoint_completed();
      #endif /* BOINC_APPLICATION */
    }

  #if BOINC_APPLICATION
    boinc_fraction_done(st->tnow / ctx->timeEvolve);
  #endif /* BOINC_APPLICATION */
}

#ifndef _WIN32
static void launchVisualizer(NBodyState* st)
{
    static const char* const argv[] = { NBODY_GRAPHICS_NAME, NULL };
    pid_t pid;
    char* path = NULL;
    char* newPath = NULL;

    pid = fork();
    if (pid != 0)  /* Parent */
        return;

    /* Child */
    if (shmdt(st->scene) < 0)
    {
        /* Hack to close the shared memory access we inherit so we
         * don't count it when the visualizer actually opens it again */
        perror("Detaching child from shared");
        return;
    }

    /* Put places convenient for testing. Not essential, failure of
     * any of these is OK */
    path = getenv("PATH");
    if (!path)
    {
        perror("Error getting PATH");
    }
    else
    {
        if (asprintf(&newPath, ".:../bin/:%s", path) < 0)
        {
            perror("Appending to path");
        }
        else
        {
            if (setenv("PATH", newPath, TRUE) < 0)
            {
                perror("Error setting PATH");
            }
            free(newPath);
        }
    }

    if (execvp(argv[0], argv) < 0)
    {
        perror("Failed to launch visualizer");
    }
}

#else

static void launchVisualizer(NBodyState* st)
{
    warn("Launching visualizer from main application not unimplemented on Windows\n");
}

#endif /* _WIN32 */

static int runSystem(const NBodyCtx* ctx, NBodyState* st, int visualizer)
{
    const real tstop = ctx->timeEvolve - ctx->timestep / 1024.0;

    if (visualizer)
        launchVisualizer(st);

    while (st->tnow < tstop)
    {
        if (visualizerIsAttached(st))
        {
            updateDisplayedBodies(st);
        }

        if (stepSystem(ctx, st))   /* advance N-body system */
            return 1;

        nbodyCheckpoint(ctx, st);
    }

    if (BOINC_APPLICATION || ctx->checkpointT >= 0)
    {
        mw_report("Making final checkpoint\n");
        if (writeCheckpoint(ctx, st))
            return warn1("Failed to write final checkpoint\n");
    }

    return 0;
}

static void endRun(NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf, const real chisq)
{
    finalOutput(ctx, st, nbf, chisq);
    destroyNBodyState(st);
}

static NBodyStatus setupRun(NBodyCtx* ctx, NBodyState* st, HistogramParams* hp, const NBodyFlags* nbf)
{
    /* If the checkpoint exists, try to use it */
    if (nbf->ignoreCheckpoint || !resolvedCheckpointExists(st))
    {
        if (setupNBody(ctx, st, hp, nbf))
            return warn1("Failed to read input parameters file\n");
    }
    else
    {
        mw_report("Checkpoint exists. Attempting to resume from it.\n");

        if (nbf->inputFile && !BOINC_APPLICATION)
            warn("Warning: input file '%s' unused\n", nbf->inputFile);

        if (readCheckpoint(ctx, st))
        {
            mw_report("Failed to read checkpoint\n");
            destroyNBodyState(st);
            return NBODY_CHECKPOINT_ERROR;
        }
        else
        {
            mw_report("Successfully read checkpoint\n");
        }
    }

    return gravMap(ctx, st); /* Start 1st step */
}

/* Set context fields read from command line flags */
static inline void nbodySetCtxFromFlags(NBodyCtx* ctx, const NBodyFlags* nbf)
{
    ctx->checkpointT = nbf->checkpointPeriod;
}

static int verifyFile(const NBodyFlags* nbf)
{
    int rc;
    NBodyCtx ctx  = EMPTY_NBODYCTX;
    NBodyState st = EMPTY_NBODYSTATE;

    rc = setupNBody(&ctx, &st, &ctx.histogramParams, nbf);
    if (rc)
        warn("File failed\n");
    else
    {
        warn("File is OK\n");
        printNBodyCtx(&ctx);
        printHistogramParams(&ctx.histogramParams);
    }

    destroyNBodyState(&st);

    return rc;
}

/* FIXME */
static NBodyCtx _ctx = EMPTY_NBODYCTX;
static NBodyState _st = EMPTY_NBODYSTATE;

int runNBodySimulation(const NBodyFlags* nbf)
{
    NBodyCtx* ctx = &_ctx;
    NBodyState* st = &_st;

    int rc = 0;
    real chisq;
    double ts = 0.0, te = 0.0;

    if (nbf->verifyOnly)
        return verifyFile(nbf);

    if (resolveCheckpoint(st, nbf->checkpointFileName))
        return warn1("Failed to resolve checkpoint\n");

    if (setupRun(ctx, st, &ctx->histogramParams, nbf))
        return warn1("Failed to setup run\n");

    nbodySetCtxFromFlags(ctx, nbf);
    if (initOutput(st, nbf))
        return warn1("Failed to open output files\n");

    if (nbf->printTiming)     /* Time the body of the calculation */
        ts = mwGetTime();


    if (createSharedScene(st, nbf->inputFile))
    {
        return warn1("Failed to create shared scene\n");
    }

    rc = runSystem(ctx, st, nbf->visualizer);
    if (rc)
        return warn1("Error running system\n");

    if (nbf->printTiming)
    {
        te = mwGetTime();
        printf("<run_time> %g </run_time>\n", te - ts);
    }

    /* Get the likelihood */
    chisq = nbodyChisq(ctx, st, nbf, &ctx->histogramParams);
    if (isnan(chisq))
    {
        warn("Failed to calculate chisq\n");
        rc = 1;
    }

    endRun(ctx, st, nbf, chisq);

    return rc;
}

