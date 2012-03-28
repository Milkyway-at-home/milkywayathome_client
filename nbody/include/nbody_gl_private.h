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

#ifndef _NBODY_GL_PRIVATE_H_
#define _NBODY_GL_PRIVATE_H_

#include "nbody_gl_includes.h"

struct SceneData
{
    float currentTime;
    float timeEvolve;
    glm::vec3 centerOfMass;
    bool staticScene;

    SceneData(bool isStatic) : currentTime(0.0f),
                               timeEvolve(0.0f),
                               centerOfMass(glm::vec3(0.0f, 0.0f, 0.0f)),
                               staticScene(isStatic) { }
};


extern glm::mat4 cameraToClipMatrix;
extern glm::mat4 textCameraToClipMatrix;

#endif /* _NBODY_GL_PRIVATE_H_ */

