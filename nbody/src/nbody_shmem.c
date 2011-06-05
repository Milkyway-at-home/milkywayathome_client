/*
Copyright (c) 2011 Matthew Arsenault
Copyright (c) 2011 Rensselaer Polytechnic Institute.

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
#include "nbody_shmem.h"

#if USE_SHMEM
  #include <sys/shm.h>
#endif


static const char nbodyGraphicsName[] = NBODY_GRAPHICS_NAME;

static void prepareSceneFromState(const NBodyCtx* ctx, const NBodyState* st)
{
    st->scene->nbody = st->nbody;
    st->scene->info.timeEvolve = (float) ctx->timeEvolve;
    st->scene->drawGalaxy = (ctx->potentialType == EXTERNAL_POTENTIAL_DEFAULT);
}


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

int createSharedScene(NBodyState* st, const NBodyCtx* ctx, const char* inputFile)
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
    prepareSceneFromState(ctx, st);

    return 0;
}

#elif USE_BOINC_SHMEM

int createSharedScene(NBodyState* st, const NBodyCtx* ctx, const char* inputFile)
{
    size_t size = sizeof(scene_t) + st->nbody * sizeof(FloatPos);

    st->scene = (scene_t*) mw_graphics_make_shmem(NBODY_BIN_NAME, (int) size);
    if (!st->scene)
    {
        warn("Failed to get shmem of size %d\n", (int) size);
        return 1;
    }

    memset(st->scene, 0, size);
    prepareSceneFromState(ctx, st);

    return 0;
}

#else

int createSharedScene(NBodyState* st, const char* inputFile)
{
    warn("Creating shared scene unimplemented for this system\n");
    return 0;
}

int visualizerIsAttached(const NBodyState* st)
{
    return 0;
}

#endif /* USE_SHMEM */


#ifndef _WIN32

void launchVisualizer(NBodyState* st, const char* visArgs)
{
    pid_t pid;
    char* path = NULL;
    char* newPath = NULL;
    int argc = 0;
    char* buf = NULL;
    char* p = NULL;
    char** argv = NULL;
    size_t argvSize = 0;
    size_t visArgsLen = 0;

    if (!st->scene) /* If there's no scene to share, there's no point */
        return;

    pid = fork();
    if (pid != 0)  /* Parent */
        return;

    /* Child */

    /* Hack to close the shared memory access we inherit so we
     * don't count it when the visualizer actually opens it again */
    if (detachSharedScene(st))
    {
        warn("Error detaching child from shared");
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

    /* Stick the program name at the head of the arguments passed in */
    visArgsLen = visArgs ? strlen(visArgs) : 0;
    argvSize = visArgsLen + sizeof(nbodyGraphicsName) + 2; /* arguments + program name + space + null */
    buf = mwCalloc(argvSize, sizeof(char));

    p = stpcpy(buf, nbodyGraphicsName);
    p = stpcpy(p, " ");
    if (visArgs)
    {
        stpcpy(p, visArgs);
    }

    if (poptParseArgvString(buf, &argc, (const char***) &argv))
    {
        free(buf);
        free(path);
        warn("Error parsing arguments for visualizer '%s'\n", visArgs);
        return;
    }

    if (execvp(argv[0], argv) < 0)
    {
        perror("Failed to launch visualizer");
    }

    free(buf);
    free(argv);

    mw_finish(EXIT_SUCCESS);  /* Unnecessary */
}

#else

void launchVisualizer(NBodyState* st, const char* visArgs)
{
    PROCESS_INFORMATION pInfo;
    STARTUPINFO startInfo;
    size_t visArgsLen, argvSize;
    char* buf;

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
        warn("Error creating visualizer process: %ld\n", GetLastError());
    }

    free(buf);
}

#endif /* _WIN32 */

void updateDisplayedBodies(NBodyState* st)
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
    scene->info.currentTime = (float) st->tnow;

    /* Read data if not paused. No copying when no screensaver attached */
    if (scene->attached && scene->usleepcount >= scene->dt && (!scene->paused || scene->step))
    {
        scene->usleepcount = 0.0;
        scene->step = FALSE;

      #ifdef _OPENMP
        #pragma omp parallel for private(i, b) schedule(static)
      #endif
        for (i = 0; i < nbody; ++i)
        {
            b = &st->bodytab[i];
            r[i].x = (float) X(Pos(b));
            r[i].y = (float) Y(Pos(b));
            r[i].z = (float) Z(Pos(b));
            r[i].ignore = ignoreBody(b);
        }

        nbodyGraphicsSetOff(&scene->attached);
    }

    scene->changed = TRUE;
}


