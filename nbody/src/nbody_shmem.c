/*
 * Copyright (c) 2011, 2012 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody.h"
#include "nbody_priv.h"
#include "milkyway_util.h"
#include "nbody_grav.h"
#include "nbody_show.h"
#include "nbody_lua.h"
#include "nbody_shmem.h"
#include "nbody_defaults.h"

#if USE_SHMEM
  #if HAVE_SYS_MMAN_H
    #include <sys/mman.h>
  #endif

  #if HAVE_SYS_STAT_H
    #include <sys/stat.h>
  #endif

  #if HAVE_SYS_WAIT_H
    #include <sys/wait.h>
  #endif

  #if HAVE_FCNTL_H
    #include <fcntl.h>
  #endif

  #include <errno.h>
  #include <err.h>
#endif /* USE_SHMEM */


#define MAX_INSTANCES 256


static const char nbodyGraphicsName[] = NBODY_GRAPHICS_NAME;

static void nbPrepareSceneFromState(const NBodyCtx* ctx, const NBodyState* st)
{
    st->scene->nbodyMajorVersion = NBODY_VERSION_MAJOR;
    st->scene->nbodyMinorVersion = NBODY_VERSION_MINOR;
    st->scene->nbody = st->nbody;
    st->scene->nSteps = ctx->nStep;
    st->scene->hasInfo = TRUE;
    st->scene->hasGalaxy = (ctx->potentialType == EXTERNAL_POTENTIAL_DEFAULT);
}

static size_t nbFindShmemSize(int nbody)
{
    size_t snapshotSize = sizeof(NBodyCircularQueue) + nbody * sizeof(FloatPos);
    return sizeof(scene_t) + NBODY_CIRC_QUEUE_SIZE * snapshotSize;
}

#if USE_SHMEM

/* Map the header of a shared memory segment and see if it is owned by
 * an active process already. */
static int nbSegmentIsOwned(int shmId)
{
    void* p;
    int owner;
    struct stat sb;

    if (fstat(shmId, &sb) < 0)
    {
        return TRUE; /* Can't be sure, assume it is owned */
    }

    if (sb.st_size < (off_t) sizeof(scene_t))
    {
        /* This is impossibly small so assume it isn't valid */
        return FALSE;
    }

    p = mmap(NULL, sizeof(scene_t), PROT_READ, MAP_SHARED, shmId, 0);
    if (p == MAP_FAILED)
    {
        mwPerror("mmap: Failed to mmap shared memory for ownership check");
        return FALSE;
    }

    /* We have the tip of the existing segment mapped. Test if this is actually in use */
    owner = OPA_load_int(&((scene_t*) p)->ownerPID);
    munmap(p, sizeof(scene_t));

    return mwProcessIsAlive(owner);
}

/* Resize and map the fully sized shared memory segment */
static scene_t* nbMapSharedSegment(const char* name, int shmId, size_t size)
{
    void* p;

    if (ftruncate(shmId, (off_t) size) < 0) /* Make the segment the correct size */
    {
        mwPerror("Error ftruncate() shared memory '%s' to size "ZU, name, size);
        return NULL;
    }

    p = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, shmId, 0);
    if (p == MAP_FAILED)
    {
        mwPerror("Error mapping shared memory segment '%s'", name);
        return NULL;
    }

    return (scene_t*) p;
}

/* Resize a segment by destroying and recreating it to work around OS X */
static int nbRecreateSegment(const char* name, int shmIdOld, int mode)
{
    int shmId;
    struct stat sb;

    if (!fstat(shmIdOld, &sb))
    {
        /* Assume that if the segment is empty to start it was new */
        if (sb.st_size == 0)
        {
            return shmIdOld;
        }
    }

    if (close(shmIdOld) < 0 || shm_unlink(name) < 0)
    {
        mwPerror("Error closing/unlink shared memory segment '%s' for recreation", name);
        return -1;
    }

    shmId = shm_open(name, O_CREAT | O_RDWR, mode);
    if (shmId < 0)
    {
        mwPerror("Error reopening shared memory '%s'", name);
    }

    return shmId;
}

/* Try to open shared memory segment and resize it if it isn't owned
 * by a living process */
static scene_t* nbOpenMappedSharedSegment(const char* name, size_t size)
{
    int shmId;
    const int mode = S_IRUSR | S_IWUSR | S_IRGRP;
    scene_t* scene = NULL;

    /* Try to create a new segment or test if it exists already */
    shmId = shm_open(name, O_CREAT | O_RDWR, mode);
    if (shmId < 0)
    {
        mwPerror("Error opening shared memory '%s'", name);
        return NULL;
    }

    if (nbSegmentIsOwned(shmId))
    {
        shm_unlink(name);
        return NULL;
    }

  #ifdef __APPLE__
    /* With OS X's noncomformant POSIX shm functions, ftruncate()
     * doesn't work except on segment creation, so remove the old one
     * and recreate to allow resizing it.
     *
     * Also on OS X read() doesn't work with shm_open fd's so that's
     * why we do all this nonsense just to mmap it to test for segment
     * ownership
     */
    shmId = nbRecreateSegment(name, shmId, mode);
    if (shmId < 0)
    {
        return NULL;
    }
  #endif

    /* Not owned, so resize it to what we need and map it */
    scene = nbMapSharedSegment(name, shmId, size);
    if (!scene)
    {
        shm_unlink(name);
    }

    close(shmId);
    return scene;
}

/* Create the next available segment of the form /milkyway_nbody_n n = 0 .. 127 */
int nbCreateSharedScene(NBodyState* st, const NBodyCtx* ctx)
{
    int pid;
    int instanceId;
    char name[128];
    scene_t* scene = NULL;
    size_t size = nbFindShmemSize(st->nbody);

    /* Try looking for the next available segment of the form /milkyway_nbody_<n> */
    for (instanceId = 0; instanceId < MAX_INSTANCES; ++instanceId)
    {
        if (snprintf(name, sizeof(name), "/milkyway_nbody_%d", instanceId) == sizeof(name))
            mw_panic("Buffer too small for shared memory name\n");

        scene = nbOpenMappedSharedSegment(name, size);
        if (scene)
            break;
    }

    if (!scene)
    {
        mw_printf("Could not open new shm segment in %d tries\n", instanceId);
        return 1;
    }


    pid = (int) getpid();
    mw_report("Process %d created scene instance %d\n", pid, instanceId);

    /* Wipe out any possibly remaining queue state, etc. */
    memset(scene, 0, sizeof(scene_t));

    st->scene = (scene_t*) scene;
    st->scene->instanceId = instanceId;
    OPA_store_int(&st->scene->ownerPID, pid);
    strncpy(st->scene->shmemName, name, sizeof(st->scene->shmemName));
    nbPrepareSceneFromState(ctx, st);

    return 0;
}

#elif USE_BOINC_SHMEM

int nbCreateSharedScene(NBodyState* st, const NBodyCtx* ctx)
{
    size_t size = nbFindShmemSize(st->nbody);

    st->scene = (scene_t*) mw_graphics_make_shmem(NBODY_BIN_NAME, (int) size);
    if (!st->scene)
    {
        mw_printf("Failed to get shmem of size %d\n", (int) size);
        return 1;
    }

    memset(st->scene, 0, size);
    nbPrepareSceneFromState(ctx, st);

    return 0;
}

#else

int nbCreateSharedScene(NBodyState* st, const NBodyCtx* ctx)
{
    (void) st, (void) ctx;
    mw_printf("Creating shared scene unimplemented for this system\n");
    return 0;
}

int visualizerIsAttached(const NBodyState* st)
{
    (void) st;
    return 0;
}

#endif /* USE_SHMEM */


#ifndef _WIN32

void nbLaunchVisualizer(NBodyState* st, const char* visArgs)
{
    pid_t pid;
    const char* path = NULL;
    char* newPath = NULL;
    int argc = 0;
    char* buf = NULL;
    char* p = NULL;
    char** argv = NULL;
    size_t argvSize = 0;
    size_t visArgsLen = 0;
    char idArg[128];

    if (!st->scene) /* If there's no scene to share, there's no point */
        return;

    pid = fork();
    if (pid != 0)  /* Parent */
    {
        int status;
        pid_t result;
        int attached;

        /* Wait until the visualizer has exited or successfully
         * attached
         * TODO: maybe this should timeout?
         */
        do
        {
            mwMilliSleep(10);

            attached = !!OPA_load_int(&st->scene->attachedPID);
            result = waitpid(pid, &status, WNOHANG);
            if (result < 0)
            {
                perror("Error waiting for child process");
                return;
            }
        }
        while (!attached && (result == 0));

        return;
    }

    /* Child */

    /* Put places convenient for testing. Not essential, failure of
     * any of these is OK */
    path = getenv("PATH");
    if (!path)
    {
        mwPerror("Error getting PATH");
    }
    else
    {
        if (asprintf(&newPath, ".:../bin/:%s", path) < 0)
        {
            mwPerror("Appending to path");
        }
        else
        {
            if (setenv("PATH", newPath, TRUE) < 0)
            {
                mwPerror("Error setting PATH");
            }
            free(newPath);
        }
    }

    if (snprintf(idArg, sizeof(idArg), "--instance-id=%d ", st->scene->instanceId) == sizeof(idArg))
        mw_panic("Buffer too small for --instance-id visualizer argument\n");

    /* Stick the program name at the head of the arguments passed in */
    visArgsLen = visArgs ? strlen(visArgs) : 0;
    argvSize = visArgsLen + sizeof(idArg) + sizeof(nbodyGraphicsName) + 2; /* arguments + program name + space + null */
    buf = mwCalloc(argvSize, sizeof(char));

    p = stpcpy(buf, nbodyGraphicsName);
    p = stpcpy(p, " ");
    p = stpcpy(p, idArg);
    if (visArgs)
    {
        stpcpy(p, visArgs);
    }

    if (poptParseArgvString(buf, &argc, (const char***) &argv))
    {
        mw_printf("Error parsing arguments for visualizer '%s'\n", visArgs);
        free(buf);
        return;
    }

    if (execvp(argv[0], argv) < 0)
    {
        mwPerror("Failed to launch visualizer '%s'", argv[0]);
    }

    free(buf);
    free(argv);
    mw_finish(EXIT_SUCCESS);  /* Unnecessary */
}

#else

void nbLaunchVisualizer(NBodyState* st, const char* visArgs)
{
    PROCESS_INFORMATION pInfo;
    STARTUPINFO startInfo;
    size_t visArgsLen, argvSize;
    char* buf;
    int attached;
    int ret;

    (void) st;

    memset(&pInfo, 0, sizeof(pInfo));
    memset(&startInfo, 0, sizeof(startInfo));
    startInfo.cb = sizeof(startInfo);

    visArgsLen = visArgs ? strlen(visArgs) : 0;
    argvSize = visArgsLen + sizeof(nbodyGraphicsName) + 2; /* arguments + program name + space + null */
    buf = mwCalloc(argvSize, sizeof(char));

    strcat(buf, nbodyGraphicsName);
    strcat(buf, " ");
    if (visArgs)
    {
        strcat(buf, visArgs);
    }

    if (!CreateProcess(NULL,
                       buf,
                       NULL,
                       NULL,
                       FALSE,
                       NORMAL_PRIORITY_CLASS,
                       NULL,
                       NULL,
                       &startInfo,
                       &pInfo))
    {
        mw_printf("Error creating visualizer process: %ld\n", GetLastError());
    }

    /* Wait for child to attach or exit */
    do
    {
        ret = WaitForSingleObject(pInfo.hProcess, 10);
        attached = !!OPA_load_int(st->scene->attachedPID);
    }
    while (!attached && (ret == WAIT_TIMEOUT));

    free(buf);
}

#endif /* _WIN32 */

static void nbWriteSnapshot(NBodyCircularQueue* queue, int buffer, const NBodyCtx* ctx, const NBodyState* st)
{
    int i;
    const Body* b;
    mwvector cmPos;
    int nbody = st->nbody;
    SceneInfo* info = &queue->info[buffer];
    FloatPos* r = &queue->bodyData[buffer * nbody];

    if (st->tree.root)
    {
        cmPos = Pos(st->tree.root);
    }
    else
    {
        /* If we are using exact nbody or haven't constructed the tree
         * yet we need to calculate the center of mass on our own */
        cmPos = nbCenterOfMass(st);
    }

    info->currentTime = (float) st->step * (float) ctx->timestep;
    info->timeEvolve = (float) ctx->timeEvolve;

    info->rootCenterOfMass[0] = (float) X(cmPos);
    info->rootCenterOfMass[1] = (float) Y(cmPos);
    info->rootCenterOfMass[2] = (float) Z(cmPos);

  #ifdef _OPENMP
    #pragma omp parallel for private(i, b) shared(r) schedule(static)
  #endif
    for (i = 0; i < nbody; ++i)
    {
        b = &st->bodytab[i];
        r[i].x = (float) X(Pos(b));
        r[i].y = (float) Y(Pos(b));
        r[i].z = (float) Z(Pos(b));
        r[i].ignore = ignoreBody(b);
    }
}

static int nbPushCircularQueue(NBodyCircularQueue* queue, const NBodyCtx* ctx, const NBodyState* st)
{
    int head, tail, nextTail;
    tail = OPA_load_int(&queue->tail);
    head = OPA_load_int(&queue->head);

    nextTail = (tail + 1) % NBODY_CIRC_QUEUE_SIZE;
    if (nextTail != head)
    {
        nbWriteSnapshot(queue, tail, ctx, st);
        OPA_store_int(&queue->tail, nextTail);
        return TRUE;
    }
    else /* write failed, queue full */
    {
        return FALSE;
    }
}

static void nbReleaseSceneLocks(scene_t* scene)
{
    OPA_store_int(&scene->paused, 0);
    OPA_store_int(&scene->blockSimulationOnGraphics, 0);
    OPA_store_int(&scene->attachedLock, 0);
    OPA_store_int(&scene->attachedPID, 0);
}

/* TODO: Should we quit if we are supposed to be blocking on graphics
 * and the graphics dies? */
NBodyStatus nbUpdateDisplayedBodies(const NBodyCtx* ctx, NBodyState* st)
{
    int pid;
    scene_t* scene = st->scene;

    if (!scene)
    {
        return NBODY_SUCCESS;
    }

    /* No copying when no screensaver attached */
    pid = OPA_load_int(&scene->attachedPID);
    if (pid == 0)
    {
        return NBODY_SUCCESS;
    }

    if (OPA_load_int(&scene->blockSimulationOnGraphics))
    {
        int updated = FALSE;
        double dt = -1.0;
        double startTime = mwGetTime();

        /* FIXME: This needs to change if we allow changing
         * whether the simulation should block
         */

        do
        {
            int attempt = 0;

            /* Keep trying to push to the queue as long as the
             * process is still attached. */
            while (   !(updated = nbPushCircularQueue(&scene->queue, ctx, st))
                   && ((pid = OPA_load_int(&scene->attachedPID)) != 0)
                   && (attempt < NBODY_QUEUE_WAIT_PERIODS))
            {
                mwMilliSleep(NBODY_QUEUE_SLEEP_INTERVAL);
                ++attempt;
            }

            /* We did update successfully, we are done */
            if (updated)
            {
                return NBODY_SUCCESS;
            }

            /* The PID was reset by the graphics to 0, the process quit normally */
            if (pid == 0)
            {
                mw_report("Graphics process quit while waiting (waited %f seconds)\n",
                          mwGetTime() - startTime);
                return NBODY_GRAPHICS_DEAD;
            }

            /* It's been a while, so the graphics may have
             * crashed. If the process is dead, we give up and
             * remove it's locks */
            if (!mwProcessIsAlive(pid))
            {
                mw_report("Graphics process %d is dead, releasing its locks\n", pid);
                nbReleaseSceneLocks(scene);
                return NBODY_GRAPHICS_DEAD;
            }

            /* FIXME: What if it goes into an infinite loop while
             * paused? */
        }
        while (OPA_load_int(&scene->paused) || (dt = mwGetTime() - startTime) < NBODY_QUEUE_TIMEOUT);

        /* The process is supposedly still alive, but it has
         * been a long time and it isn't paused.
         * TODO: Should we kill it?
         */
        mw_report("Waiting on graphics timed out (%f seconds)\n", mwGetTime() - startTime);

        return NBODY_GRAPHICS_TIMEOUT;
    }
    else
    {
        nbPushCircularQueue(&scene->queue, ctx, st);
    }

    return NBODY_SUCCESS;
}

/* Report the simulation has ended as a hint to the graphics to quit */
void nbReportSimulationComplete(NBodyState* st)
{
    scene_t* scene = st->scene;

    if (!scene)
        return;

    /* Disown the scene */
    OPA_store_int(&scene->ownerPID, 0);
}

