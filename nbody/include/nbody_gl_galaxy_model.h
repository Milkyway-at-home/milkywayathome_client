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

#ifndef _NBODY_GL_GALAXY_MODEL_H_
#define _NBODY_GL_GALAXY_MODEL_H_

#include "nbody_gl_includes.h"

struct NBodyVertex
{
    GLfloat x, y, z;

    NBodyVertex() : x(0.0f), y(0.0f), z(0.0f) { }
    NBodyVertex(GLfloat xx, GLfloat yy, GLfloat zz) : x(xx), y(yy), z(zz) { }
};

class GalaxyModel
{
private:
    NBodyVertex* points;
    GLuint nPoints;
    GLuint count;

    GLuint vao;
    GLuint buffer;

    float bulgeRadius;
    float bulgeHeight;
    float diskScale;
    float diskCoeff;
    float diskEdgeZ;

    // what the rounded out diameter will be
    // after producing shape so we can sample the texture
    float invActualDiameter;

    float totalDiameter;
    GLint radialSlices;
    GLint axialSlices;
    float axialSliceSize;
    float diameterSlice;

    GLuint texture;

    struct GalaxyProgramData
    {
        GLuint program;

        GLint positionLoc;
        GLint modelToCameraMatrixLoc;
        GLint cameraToClipMatrixLoc;
        GLint galaxyTextureLoc;
        GLint invGalaxyDiameterLoc;
    } programData;

    void makePoint(NBodyVertex& point, bool neg, float r, float theta);
    void generateSegment(bool neg);
    void setMilkywayModelParameters();
    float shapeFunction(float r);
    void loadGalaxyTexture();
    void loadShaders();
    void prepareVAO();
    void bufferData();
    void generateModel();

    size_t size() const
    {
        return this->nPoints * sizeof(NBodyVertex);
    }

public:
    void draw(const glm::mat4& modelMatrix) const;

    GalaxyModel();
    ~GalaxyModel();
};

#endif /* _NBODY_GL_GALAXY_MODEL_H_ */


