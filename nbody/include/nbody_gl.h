/*
 * Copyright (c) 2012 Matthew Arsenault
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

#ifndef _NBODY_GL_H_
#define _NBODY_GL_H_

#include "nbody_graphics.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    int fullscreen;
    int plainFullscreen;
    int width;
    int height;
    int blockSimulation;

    int monochrome;
    int untexturedPoints;
    int originCenter;
    int drawAxes;
    int noFloat;

    int pid;
    char* file;
    int instanceId;
} VisArgs;

#define EMPTY_VIS_ARGS { FALSE, FALSE, 0, 0, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 0, NULL, -1 }

int nbglRunGraphics(scene_t* scene, const VisArgs* args);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_GL_H_ */


