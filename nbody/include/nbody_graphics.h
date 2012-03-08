/*
 * Copyright (c) 2002-2006 John M. Fregeau, Richard Campbell, Jeff Molofee
 * Copyright (c) 2011 Matthew Arsenault
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

#define MAX_DRAW_TRACE_POINTS 256

/* "mw_nbody" */
#define DEFAULT_SHMEM_KEY ((key_t) 0x6d775f6e626f6479)

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
} SceneInfo;

typedef enum
{
    MOUSE_MODE_NONE,
    MOUSE_MODE_MOVE,
    MOUSE_MODE_ZOOM
} NBodyMouseMode;

#ifndef NAME_MAX
  #define NAME_MAX 255
#endif

/* the scene structure */
typedef struct
{
    char shmemName[NAME_MAX + 1];
    int instanceId;
    int nbodyMajorVersion;
    int nbodyMinorVersion;

    int nbody;
    int hasGalaxy;

    /* Terrible way of checking that a screensaver is attached.  It
       avoids copying all of the bodies if it isn't running. Also only
       semi-functional way of making sure only one screensaver is
       attached at a time to the simulation. There's a small window
       where you can end up with multiple at a time which could
       potentially break, but it shouldn't matter that much.
     */
    OPA_int_t attachedCount;


    float rootCenterOfMass[3];     /* Center of mass of the system  */
    float startingPositionHint[3];
    float startingAngleHint[3];
    SceneInfo info;

    int currentTracePoint;
    FloatPos orbitTrace[MAX_DRAW_TRACE_POINTS];
    FloatPos rTrace[];
} scene_t;

#if defined(__GNUC__) && (__GNUC__ >= 4 && __GNUC_MINOR__ >= 1)
  // #define nbodyGraphicsAtomicIncrement(x) __sync_fetch_and_add((x), 1)
  // #define nbodyGraphicsAtomicDecrement(x) __sync_fetch_and_sub((x), 1)
    #define nbodyGraphicsSetOn(x) __sync_fetch_and_or((x), 1)
    #define nbodyGraphicsSetOff(x) __sync_fetch_and_and((x), 0)
    #define nbodyGraphicsTestVal(x) __sync_add_and_fetch((x), 0)
#elif defined(_MSC_VER)
  // #define nbodyGraphicsAtomicIncrement(x) InterlockedIncrement((x))
  // #define nbodyGraphicsAtomicDecrement(x) InterlockedDecrement((x))

  #define nbodyGraphicsSetOn(x) InterlockedBitTestAndSet((x), 1) // 1st bit
  #define nbodyGraphicsSetOff(x) InterlockedBitTestAndReset((x), 1)
  #define nbodyGraphicsTestVal(x) InterlockedOr((x), 0);
#else
  #error Need atomics for compiler
#endif /* defined(__GNUC__) */

#endif /* _NBODY_GRAPHICS_H_ */


