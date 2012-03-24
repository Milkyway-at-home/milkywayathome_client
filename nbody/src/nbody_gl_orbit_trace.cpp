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

#include "nbody_gl_orbit_trace.h"
#include "nbody_gl_shaders.h"
#include "nbody_gl_util.h"
#include "nbody_gl_private.h"

OrbitTrace::OrbitTrace(const scene_t* scene) : nPoints(0), maxPoints(scene->nSteps)
{
    this->loadShader();
    this->createBuffer();
    this->prepareVAO();
}

OrbitTrace::~OrbitTrace()
{
    glDeleteProgram(this->progData.program);
    glDeleteBuffers(1, &this->cmPosBuffer);
    glDeleteVertexArrays(1, &this->vao);
}

void OrbitTrace::loadShader()
{
    this->progData.program = nbglCreateProgram("orbit trace program",
                                               (const char*) orbit_trace_vertex_glsl,
                                               (const char*) orbit_trace_fragment_glsl,
                                               (GLint) orbit_trace_vertex_glsl_len,
                                               (GLint) orbit_trace_fragment_glsl_len);
    GLuint program = this->progData.program;

    this->progData.cmPosLoc = glGetAttribLocation(program, "cmPos");
    this->progData.modelToCameraMatrixLoc = glGetUniformLocation(program, "modelToCameraMatrix");
    this->progData.cameraToClipMatrixLoc = glGetUniformLocation(program, "cameraToClipMatrix");
}

void OrbitTrace::drawTrace(const glm::mat4& modelMatrix)
{
    glUseProgram(this->progData.program);
    glUniformMatrix4fv(this->progData.modelToCameraMatrixLoc, 1, GL_FALSE, glm::value_ptr(modelMatrix));
    glUniformMatrix4fv(this->progData.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(cameraToClipMatrix));

    glDepthMask(GL_FALSE);

    glBindVertexArray(this->vao);
    glDrawArrays(GL_LINE_STRIP, 0, this->nPoints);

    glDepthMask(GL_TRUE);
    glBindVertexArray(0);
    glUseProgram(0);
}

void OrbitTrace::createBuffer()
{
    glGenBuffers(1, &this->cmPosBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, this->cmPosBuffer);
    glBufferData(GL_ARRAY_BUFFER, 4 * this->maxPoints * sizeof(GLfloat), NULL, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void OrbitTrace::prepareVAO()
{
    glGenVertexArrays(1, &this->vao);
    glBindVertexArray(this->vao);
    glEnableVertexAttribArray(this->progData.cmPosLoc);
    glBindBuffer(GL_ARRAY_BUFFER, this->cmPosBuffer);
    glVertexAttribPointer(this->progData.cmPosLoc, 4, GL_FLOAT, GL_FALSE, 0, 0);
    glBindVertexArray(0);
}

void OrbitTrace::addPoint(const GLfloat cm[3])
{
    if (this->nPoints >= this->maxPoints)
        return;

    GLfloat cm4[4] = { cm[0], cm[1], cm[2], 1.0f };

    glBindBuffer(GL_ARRAY_BUFFER, this->cmPosBuffer);
    glBufferSubData(GL_ARRAY_BUFFER,
                    this->nPoints * 4 * sizeof(GLfloat),
                    4 * sizeof(GLfloat),
                    cm4);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    this->nPoints++;
}

