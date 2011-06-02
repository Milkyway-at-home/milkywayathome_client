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

/* "mw_nbody" */
#define DEFAULT_SHMEM_KEY ((key_t) 0x6d775f6e626f6479)

typedef struct
{
    float x, y, z;
    int ignore;
} FloatPos;

/* Mostly for progress information */
typedef struct
{
    float currentTime;
    float timeEvolve;
} SceneInfo;

/* the scene structure */
typedef struct
{
    int nbody;
    int drawGalaxy;
    float z;
    float xrot;
    float yrot;
    float starsize;

    /* Terrible way of checking that a screensaver is attached.  It
       avoids copying all of the bodies if it isn't running. Also only
       semi-functional way of making sure only one screensaver is
       attached at a time to the simulation. There's a small window
       where you can end up with multiple at a time which could
       potentially break, but it shouldn't matter that much.
     */
  #ifndef _MSC_VER
    volatile int attached;
  #else
    volatile LONG attached;
  #endif /* _MSC_VER */

    int fullscreen;
    int drawaxes;
    int ntri;
    int paused;
    int step;
    int mousemode;
    int changed;
    double t;
    double dt;
    double usleepcount;
    double usleepdt;

    float startingPositionHint[3];
    float startingAngleHint[3];
    SceneInfo info;
    FloatPos r[];
} scene_t;

#if defined(__GNUC__)
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


