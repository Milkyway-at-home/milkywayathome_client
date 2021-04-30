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

#include "nbody_config.h"

#include "nbody_gl_galaxy_model.h"
#include "nbody_gl_shaders.h"
#include "nbody_gl_util.h"
#include "nbody_milkyway_image.h"
#include "nbody_gl_private.h"
#include "milkyway_math.h"

GalaxyModel::GalaxyModel()
{
    this->loadShaders();
    this->loadGalaxyTexture();
    this->setMilkywayModelParameters();
    this->generateModel();
    this->bufferData();
    this->prepareVAO();
}

GalaxyModel::~GalaxyModel()
{
    glDeleteProgram(this->programData.program);
    glDeleteVertexArrays(1, &this->vao);
    glDeleteBuffers(1, &this->buffer);
    glDeleteTextures(1, &this->texture);

    delete[] this->points;
}

void GalaxyModel::loadGalaxyTexture()
{
    size_t imgSize = milkywayImageParams.width * milkywayImageParams.height;
    unsigned char* buf = new unsigned char[imgSize * milkywayImageParams.bytes_per_pixel];

    glGenTextures(1, &this->texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    MILKYWAYIMAGE_RUN_LENGTH_DECODE(buf, milkywayImageRLEPixelData, imgSize, milkywayImageParams.bytes_per_pixel);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, milkywayImageParams.width, milkywayImageParams.height, 0, GL_RGB, GL_UNSIGNED_BYTE, buf);

    glGenerateMipmap(GL_TEXTURE_2D);

    delete[] buf;
}

void GalaxyModel::loadShaders()
{
    this->programData.program = nbglCreateProgram("galaxy program",
                                                  (const char*) galaxy_vertex_glsl,
                                                  (const char*) galaxy_fragment_glsl,
                                                  (GLint) galaxy_vertex_glsl_len,
                                                  (GLint) galaxy_fragment_glsl_len);

    this->programData.positionLoc = glGetAttribLocation(this->programData.program, "position");
    this->programData.modelToCameraMatrixLoc = glGetUniformLocation(this->programData.program, "modelToCameraMatrix");
    this->programData.cameraToClipMatrixLoc = glGetUniformLocation(this->programData.program, "cameraToClipMatrix");

    this->programData.galaxyTextureLoc = glGetUniformLocation(this->programData.program, "galaxyTexture");
    this->programData.invGalaxyDiameterLoc = glGetUniformLocation(this->programData.program, "invGalaxyDiameter");
}

void GalaxyModel::prepareVAO()
{
    glGenVertexArrays(1, &this->vao);

    glBindVertexArray(this->vao);
    glEnableVertexAttribArray(this->programData.positionLoc);
    glBindBuffer(GL_ARRAY_BUFFER, this->buffer);
    glVertexAttribPointer(this->programData.positionLoc, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindVertexArray(0);
}

void GalaxyModel::bufferData()
{
    assert(this->nPoints != 0);

    glGenBuffers(1, &this->buffer);

    glBindBuffer(GL_ARRAY_BUFFER, this->buffer);
    glBufferData(GL_ARRAY_BUFFER, this->size(), (const GLfloat*) this->points, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void GalaxyModel::setMilkywayModelParameters()
{
    this->radialSlices = 15;
    this->axialSlices = 50;
    this->axialSliceSize = (float) M_2PI / (float) this->axialSlices;

    this->bulgeRadius = 0.5f * 8.308f;
    this->bulgeHeight = 0.7f;

    this->diskScale = 7.7f;
    //this->diskCoeff = 2.5f;
    this->diskCoeff = 1.5f * M_E * 0.459899f;

    this->totalDiameter = 2.0f * 15.33f;

    GLuint nDiameterSlices = 2 * this->radialSlices;

    // for texturing
    this->invActualDiameter = 1.0f / (float) this->totalDiameter;
    this->diameterSlice = this->totalDiameter / (float) nDiameterSlices;

    this->diskEdgeZ = shapeFunction(0.5f * this->totalDiameter);
}

float GalaxyModel::shapeFunction(float r)
{
    float rp = std::fabs(r);

    float bulge = this->bulgeHeight * expf(-7.67f * pow(rp / this->bulgeRadius, 4.0f));

    float disk = this->diskCoeff * expf((-1.0f / this->diskScale) * rp);

    return bulge + disk;
}

void GalaxyModel::makePoint(NBodyVertex& point, bool neg, float r, float theta)
{
    float z;

    point.x = std::fabs(r) * cosf(theta);
    point.y = std::fabs(r) * sinf(theta);

    z = shapeFunction(r) - this->diskEdgeZ;

    point.z = neg ? -z : z;
}

// generate top or bottom half of galaxy model
void GalaxyModel::generateSegment(bool neg)
{
    GLint start, end, inc;

    if (neg)
    {
        start = this->radialSlices - 1;
        end = -1;
        inc = -1;
    }
    else
    {
        start = 0;
        end = this->radialSlices;
        inc = 1;
    }

    for (GLint i = start; i != end; i += inc)
    {
        float r, r1;

        if (neg) // next ring is smaller
        {
            r1 = i * this->diameterSlice;
            r = (std::abs(i) + 1) * this->diameterSlice;
        }
        else   // next ring is bigger
        {
            r = i * this->diameterSlice;
            r1 = (std::abs(i) + 1) * this->diameterSlice;
        }

        for (GLint j = 0; j < this->axialSlices + 1; ++j)
        {
            float theta = this->axialSliceSize * (float) j;

            makePoint(this->points[this->count++], neg, r, theta);
            makePoint(this->points[this->count++], neg, r1, theta);
        }
    }
}

void GalaxyModel::generateModel()
{
    GLuint segmentPoints = 2 * this->radialSlices * (this->axialSlices + 1);
    this->nPoints = 2 * segmentPoints;

    this->points = new NBodyVertex[this->nPoints];
    this->count = 0;

    generateSegment(false);
    generateSegment(true);

    assert(this->count == this->nPoints);
}

void GalaxyModel::draw(const glm::mat4& modelMatrix) const
{
    glUseProgram(this->programData.program);
    glUniformMatrix4fv(this->programData.modelToCameraMatrixLoc, 1, GL_FALSE, glm::value_ptr(modelMatrix));
    glUniformMatrix4fv(this->programData.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(cameraToClipMatrix));

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, this->texture);
    glUniform1i(this->programData.galaxyTextureLoc, 2);
    glUniform1f(this->programData.invGalaxyDiameterLoc, this->invActualDiameter);

    glBindVertexArray(this->vao);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, this->nPoints);

    glBindVertexArray(0);
    glUseProgram(0);
}

