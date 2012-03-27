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

#include <stdint.h>
#include <opa_primitives.h>

#ifndef NAME_MAX
  #define NAME_MAX 255
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
    float currentTime;
    float timeEvolve;
    float rootCenterOfMass[3];     /* Center of mass of the system  */
} SceneInfo;

/* This should only be used as the last element of the scene_t */
typedef struct
{
    OPA_int_t head;
    OPA_int_t tail;
    SceneInfo info[NBODY_CIRC_QUEUE_SIZE];
    FloatPos bodyData[1];
} NBodyCircularQueue;

/* the scene structure */
typedef struct
{
    char shmemName[NAME_MAX + 1];
    int instanceId;
    int ownerPID;
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
} scene_t;

#endif /* _NBODY_GRAPHICS_H_ */


