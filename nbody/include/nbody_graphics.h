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

#endif /* _NBODY_GRAPHICS_H_ */


