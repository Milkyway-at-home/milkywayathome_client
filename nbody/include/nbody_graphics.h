/*
 * Copyright (c) 2011-2012 Matthew Arsenault
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

#ifndef _NBODY_GRAPHICS_H_
#define _NBODY_GRAPHICS_H_

#include "nbody_config.h"

#include <stdint.h>
#include <opa_primitives.h>

#ifndef NAME_MAX
  #define NAME_MAX 255
#endif

#if USE_POSIX_SHMEM
  #define NBODY_SHMEM_NAME_FMT_STR "/milkyway_nbody_%d"
#elif USE_WIN32_SHARED_MAP
  #define NBODY_SHMEM_NAME_FMT_STR "Local\\milkyway_nbody_%d"
#else
  #define NBODY_SHMEM_NAME_FMT_STR "milkyway_nbody_%d"
#endif

#define NBODY_CIRC_QUEUE_SIZE 3

/* Number of milliseconds to sleep if we are blocking the simulation
 * waiting for the queue to clear */
#define NBODY_QUEUE_SLEEP_INTERVAL 16

/* number of NBODY_QUEUE_SLEEP_INTERVALS to wait before checking if the process is dead */
#define NBODY_QUEUE_WAIT_PERIODS 7

/* Number of seconds to wait on a supposedly living attached garphics
   process before giving up on it
 */
#define NBODY_QUEUE_TIMEOUT 10.0

typedef struct
{
    float x, y, z;
    int32_t ignore;
} FloatPos;

/* Mostly for progress information */
typedef struct
{
    unsigned int currentStep;
    float currentTime;
    float timeEvolve;
    float rootCenterOfMass[3];     /* Center of mass of the system  */
} SceneInfo;

typedef struct
{
    OPA_int_t head;
    OPA_int_t tail;
    SceneInfo info[NBODY_CIRC_QUEUE_SIZE];
} NBodyCircularQueue;

/* the scene structure */
typedef struct
{
    /* The Win32 API is retarded and doesn't provide a way to get the
       size of a segment so we store it ourselves first and hope it
       works out OK
    */
    size_t sceneSize;
    char shmemName[NAME_MAX + 1];
    int instanceId;
    OPA_int_t ownerPID;
    int nbodyMajorVersion;
    int nbodyMinorVersion;

    /* We require exclusive access. One visualizer per running simulation
     *
     * We need 2 copies set at different times to signal to forking
     * process that visualizer is ready with setup.
     *
     * We need to acquire the attached lock before we can set
     * settings from the graphics such as the
     * blockSimulationOnGraphics, but then we need to indicate when
     * we have finished setup or otherwise there can be a visible
     * jump at the start
     */
    OPA_int_t attachedPID;
    OPA_int_t attachedLock;
    OPA_int_t paused;

    /* Last time queue was pushed to in seconds */
    OPA_int_t lastUpdateTime;
    OPA_int_t updatePeriod;

    /*
      Optionally block the simulation while the graphics catches up.
      This would let you create smoother visualizations etc. at the
      expense of slowing the simulation.
    */
    OPA_int_t blockSimulationOnGraphics;

    int nbody;
    unsigned int nSteps;
    int hasGalaxy;
    int hasInfo;
    int staticScene;

    NBodyCircularQueue queue;
    FloatPos sceneData[1]; /* Space for orbit trace then space for actual data for the queue */
} scene_t;

/* Get the starting position of the given queue position accounting for the orbit trace offset */
static inline FloatPos* nbSceneGetQueueBuffer(scene_t* scene, int buffer)
{
    return &scene->sceneData[scene->nSteps + buffer * scene->nbody];
}

/* Get the starting position of the orbit trace in the scene data */
static inline FloatPos* nbSceneGetOrbitTrace(scene_t* scene)
{
    return &scene->sceneData[0];
}

static inline size_t nbFindShmemSize(int nbody, int nSteps)
{
    size_t snapshotSize = nbody * sizeof(FloatPos);
    return sizeof(scene_t) + nSteps * sizeof(FloatPos) + NBODY_CIRC_QUEUE_SIZE * snapshotSize;
}

#endif /* _NBODY_GRAPHICS_H_ */


