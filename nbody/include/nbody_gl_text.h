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


#ifndef _NBODY_GL_TEXT_H_
#define _NBODY_GL_TEXT_H_

#include "nbody_config.h"

#include "nbody_gl_includes.h"
#include "nbody_graphics.h"
#include "nbody_gl_private.h"
#include "Roboto_Regular_12.h"

#include <vector>

struct TextPen
{
    float x, y;
    TextPen() : x(0.0f), y(0.0f) { }
};

class NBodyTextItem
{
private:
    GLuint verticesBuffer;
    GLuint indicesBuffer;
    GLuint vao;
    const texture_font_t* font;

    PACKED
    struct TextVertex
    {
        float x, y; // position
        float s, t; // texture position
    };

    PACKED
    struct TextVertexArray
    {
        TextVertex data[4];
    };

    PACKED
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

class NBodyText
{
private:
    /* Text which may change every frame */
    NBodyTextItem* varText;
    NBodyTextItem* constText;
    const texture_font_t* font;
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

    NBodyText(const texture_font_t* textFont);
    ~NBodyText();
};

#endif /* _NBODY_GL_TEXT_H_ */

