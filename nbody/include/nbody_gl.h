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

#define DEFAULT_TEXTURED_POINT_SIZE 250.0f
#define DEFAULT_POINT_POINT_SIZE 40.0f
#define DEFAULT_FLOAT_SPEED 5.0f
#define DEFAULT_FLOAT FALSE

#define DEFAULT_EVENT_POLL_PERIOD ((int) (1000.0 / 30.0))
#define MAX_EVENT_POLL_PERIOD ((int) 5000)

#define DEFAULT_UNTEXTURED_POINTS FALSE
#define DEFAULT_SHOW_AXES FALSE


/* The center of mass points don't persist, so the simulation will
 * only show points from when the graphics started. This should be
 * fixed before enabling it for the screensaver. */
#define DEFAULT_SHOW_ORBIT_TRACE FALSE
#define DEFAULT_SHOW_INFO TRUE
#define DEFAULT_ORIGIN_CENTERED FALSE
#define DEFAULT_MONOCHROMATIC FALSE
#define DEFAULT_BLOCK_SIMULATION FALSE
#define DEFAULT_UPDATE_PERIOD 0
#define DEFAULT_QUIT_ON_COMPLETE FALSE
#define DEFAULT_PRINTFRAMES FALSE

typedef struct
{
    int fullscreen;
    int plainFullscreen;
    int width;
    int height;
    int eventPollPeriod;
    int printFrames;

    int quitOnComplete;
    int blockSimulation;
    int updatePeriod;
    int noFloat;
    float floatSpeed;
    float texturedPointSize;
    float pointPointSize;
    int untexturedPoints;
    int monochromatic;
    int originCentered;
    int noDrawInfo;
    int drawAxes;
    int drawOrbitTrace;

    int pid;
    char* file;
    int instanceId;
} VisArgs;

#define EMPTY_VIS_ARGS { FALSE, FALSE, 0, 0, 0, FALSE, 0, FALSE, FALSE, FALSE, 0.0f, 0.0f, 0.0f, FALSE, FALSE, FALSE, FALSE, FALSE, 0, NULL, -1 }

int nbglRunGraphics(scene_t* scene, const VisArgs* args);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_GL_H_ */


