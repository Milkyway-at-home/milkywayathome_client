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

#include "nbody_gl_axes.h"
#include "nbody_gl_resources.h"
#include "nbody_gl_util.h"
#include "nbody_gl_private.h"

NBodyAxes::NBodyAxes()
{
    this->loadShader();
    this->createBuffers();
    this->prepareVAO();
}

NBodyAxes::~NBodyAxes()
{
    glDeleteProgram(this->axesProgramData.program);

    glDeleteBuffers(1, &this->axesBuffer);
    glDeleteBuffers(1, &this->axesColorBuffer);
    glDeleteVertexArrays(1, &this->axesVAO);
}

void NBodyAxes::loadShader()
{
    this->axesProgramData.program = nbglCreateProgram("axes program",
                                                      (const char*) axes_vertex_glsl,
                                                      (const char*) axes_fragment_glsl,
                                                      (GLint) axes_vertex_glsl_len,
                                                      (GLint) axes_fragment_glsl_len);

    GLuint program = this->axesProgramData.program;

    this->axesProgramData.positionLoc = glGetAttribLocation(program, "position");
    this->axesProgramData.colorLoc = glGetAttribLocation(program, "inputColor");
    this->axesProgramData.modelToCameraMatrixLoc = glGetUniformLocation(program, "modelToCameraMatrix");
    this->axesProgramData.cameraToClipMatrixLoc = glGetUniformLocation(program, "cameraToClipMatrix");
}

//#define AXES_LENGTH 10.0f
//#define AXES_LENGTH 15.33f
#define AXES_LENGTH 0.459899f

void NBodyAxes::createBuffers()
{
    static const GLfloat axes[6][4] =
        {
            { 0.0f,        0.0f,        0.0f,        1.0f },
            { AXES_LENGTH, 0.0f,        0.0f,        1.0f },
            { 0.0f,        0.0f,        0.0f,        1.0f },
            { 0.0f,        AXES_LENGTH, 0.0f,        1.0f },
            { 0.0f,        0.0f,        0.0f,        1.0f },
            { 0.0f,        0.0f,        AXES_LENGTH, 1.0f }
        };

    static const GLfloat colors[6][4] =
        {
            { 1.0f, 0.0f, 0.0f, 0.75f },
            { 1.0f, 0.0f, 0.0f, 0.75f },
            { 0.0f, 1.0f, 0.0f, 0.75f },
            { 0.0f, 1.0f, 0.0f, 0.75f },
            { 0.0f, 0.0f, 1.0f, 0.75f },
            { 0.0f, 0.0f, 1.0f, 0.75f }
        };

    glGenBuffers(1, &this->axesBuffer);
    glGenBuffers(1, &this->axesColorBuffer);

    glBindBuffer(GL_ARRAY_BUFFER, this->axesBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(axes), (const GLfloat*) axes, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, this->axesColorBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(colors), (const GLfloat*) colors, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void NBodyAxes::prepareVAO()
{
    glGenVertexArrays(1, &this->axesVAO);
    glBindVertexArray(this->axesVAO);
    glEnableVertexAttribArray(this->axesProgramData.positionLoc);
    glEnableVertexAttribArray(this->axesProgramData.colorLoc);

    glBindBuffer(GL_ARRAY_BUFFER, this->axesBuffer);
    glVertexAttribPointer(this->axesProgramData.positionLoc, 4, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, this->axesColorBuffer);
    glVertexAttribPointer(this->axesProgramData.colorLoc, 4, GL_FLOAT, GL_FALSE, 0, 0);

    glBindVertexArray(0);
}

void NBodyAxes::draw(const glm::mat4& modelMatrix)
{
    glDisable(GL_DEPTH_TEST);

    glUseProgram(this->axesProgramData.program);
    glUniformMatrix4fv(this->axesProgramData.modelToCameraMatrixLoc, 1, GL_FALSE, glm::value_ptr(modelMatrix));
    glUniformMatrix4fv(this->axesProgramData.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(cameraToClipMatrix));

    glBindVertexArray(this->axesVAO);
    glDrawArrays(GL_LINES, 0, 6);
    glBindVertexArray(0);
    glUseProgram(0);
    glEnable(GL_DEPTH_TEST);
}

