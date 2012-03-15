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

#include <vector>
#include <OpenGL/gl3.h>

#include "nbody_graphics.h"
#include "Ubuntu-R.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#pragma GCC diagnostic pop

// FIXME: Move me
struct SceneData
{
    float currentTime;
    float timeEvolve;
    glm::vec3 centerOfMassView;

    SceneData() : currentTime(0.0f), timeEvolve(0.0f), centerOfMassView(glm::vec3(0.0f, 0.0f, 0.0f)) { }
};

struct TextPen
{
    float x, y;

    TextPen() : x(0.0f), y(0.0f) { }
};

class NBodyTextItem;

class NBodyText
{
private:
    /* Text which may change every frame */
    NBodyTextItem* varText;
    NBodyTextItem* constText;
    GLuint textTexture;

    TextPen penEndConst;

    struct TextProgramData
    {
        GLuint program;

        GLint positionLoc;
        GLint textTextureLoc;

        GLint modelToCameraMatrixLoc;
        GLint cameraToClipMatrixLoc;
    } textProgram;

    void loadTextTexture();
    void loadShader();

public:
    void prepareConstantText(const scene_t* scene);
    void drawProgressText(const SceneData& sceneData);

    NBodyText();
    ~NBodyText();
};

class NBodyTextItem
{
private:
    GLuint verticesBuffer;
    GLuint indicesBuffer;
    GLuint vao;
    const texture_font_t* font;

    __attribute__((packed))
    struct TextVertex
    {
        float x, y; // position
        float s, t; // texture position
    };

    __attribute__((packed))
    struct TextVertexArray
    {
        TextVertex data[4];
    };

    __attribute__((packed))
    struct TextVertexIndexArray
    {
        GLuint data[6];
    };

    std::vector<TextVertexIndexArray> indices;
    std::vector<TextVertexArray> vertices;


public:
    NBodyTextItem(GLuint textProgramPositionLoc, const texture_font_t* font);
    ~NBodyTextItem();

    void addText(const wchar_t* text, TextPen& pen);
    void uploadText();
    void clearText();
    void drawTextItem();
};

