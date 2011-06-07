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

#ifndef _NBODY_GL_H_
#define _NBODY_GL_H_

#ifndef __APPLE__
  #include <GL/glut.h>
  #include <GL/gl.h>
  #include <GL/glu.h>
#else
  #include <GLUT/glut.h>
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
#endif /* __APPLE__ */

#include "nbody_graphics.h"

typedef struct
{
    int fullscreen;
    int width;
    int height;

    int monochrome;
    int notUseGLPoints;

    /* pid_t pid */
    /* char* nbody bin */
} VisArgs;

#define EMPTY_VIS_ARGS { FALSE, 0, 0, FALSE, FALSE }


int connectSharedScene();
int nbodyGLSetup(const VisArgs* args);
void nbodyGLCleanup();

#ifndef _WIN32
int nbodyInitShmemKey(const char* progName);
#endif

#endif /* _NBODY_GL_H_ */


