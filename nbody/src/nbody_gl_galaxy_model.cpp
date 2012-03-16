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

#include "nbody_gl_galaxy_model.h"
#include "nbody_gl_util.h"
#include "nbody_gl_resources.h"
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
    glGenTextures(1, &this->texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, milkywayImage.width, milkywayImage.height, 0, GL_RGB, GL_UNSIGNED_BYTE, milkywayImage.pixel_data);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, milkywayImage.width, milkywayImage.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, milkywayImage.pixel_data);
    glGenerateMipmap(GL_TEXTURE_2D);
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

    //glBindBuffer(GL_ARRAY_BUFFER, this->colorBuffer);
    //glVertexAttribPointer(this->programData.colorLoc, 3, GL_FLOAT, GL_FALSE, 0, 0);

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
    this->axialSliceSize = M_2PI / (double) this->axialSlices;

    this->bulgeRadius = 0.5 * 8.308;
    this->smallBulgeRadius = 0.7;

    this->diskScale = 7.7;
    this->diskCoeff = 2.5;

    this->totalDiameter = 2.0 * 15.33;

    GLuint nDiameterSlices = 2 * this->radialSlices;

    // for texturing
    this->invActualDiameter = 1.0f / (float) this->totalDiameter;
    this->diameterSlice = this->totalDiameter / (double) nDiameterSlices;

    this->diskEdgeZ = diskShapeFunction(0.5 * this->totalDiameter);
}

double GalaxyModel::diskShapeFunction(double r)
{
    double rd = this->diskScale;
    double a = this->diskCoeff;

    double rp = std::fabs(r);

    if (std::fabs(r) <= this->bulgeRadius)
    {
        rp += this->smallBulgeRadius;
    }

    return a * exp((-1.0 / rd) * rp);
}

void GalaxyModel::makePoint(NBodyVertex& point, bool neg, double r, double theta)
{
    double z;

    point.x = std::fabs(r) * cos(theta);
    point.y = std::fabs(r) * sin(theta);

    double gr = this->bulgeRadius;

    if (std::fabs(r) <= gr && false)
    {
        //z = 0.25 * diskShapeFunction(r) * sqrt(sqr(gr) - sqr(point.x) - sqr(point.y));

        //z = 0.66f * std::fabs(sin(point.y / r));
        //z = gr * cos(point.x / r);

        //z = sqrt(sqr(gr) - sqr(point.x) - sqr(point.y));

        z = sqrt(sqr(this->smallBulgeRadius) - sqr(point.x) - sqr(point.y));
        printf("Bulge z %f, r = %f\n", neg ? -z : z, r);
    }
    else
    {
        z = diskShapeFunction(r);

        //printf("Disk z %f, r = %f\n", neg ? -z : z, r);
    }

    z -= this->diskEdgeZ;

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
        double r, r1;

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
            double theta = this->axialSliceSize * (double) j;

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

    //glBlendFunc(GL_SRC_COLOR, GL_ONE);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_DST_ALPHA);
    //glBlendFunc(GL_CONSTANT_COLOR_EXT, GL_ONE_MINUS_SRC_COLOR);


    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, this->texture);
    glUniform1i(this->programData.galaxyTextureLoc, 2);
    glUniform1f(this->programData.invGalaxyDiameterLoc, this->invActualDiameter);

    glBindVertexArray(this->vao);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, this->nPoints);

    glBindVertexArray(0);
    glUseProgram(0);

    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

