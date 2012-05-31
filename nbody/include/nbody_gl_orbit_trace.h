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

#ifndef _NBODY_GL_ORBIT_TRACE_H_
#define _NBODY_GL_ORBIT_TRACE_H_

#include "nbody_gl_includes.h"
#include "nbody_graphics.h"

class OrbitTrace
{
private:
    GLuint vao;
    GLuint cmPosBuffer;

    struct OrbitTraceProgramData
    {
        GLuint program;
        GLint cmPosLoc;
        GLint modelToCameraMatrixLoc;
        GLint cameraToClipMatrixLoc;
    } progData;

    unsigned int nPoints;
    unsigned int maxPoints;

    void createBuffer();
    void prepareVAO();
    void loadShader();

public:
    void drawTrace(const glm::mat4& modelMatrix);
    void updatePoints(const FloatPos* cmList, GLuint step);

    OrbitTrace(const scene_t* scene);
    ~OrbitTrace();
};

#endif /* _NBODY_GL_ORBIT_TRACE_H_ */

