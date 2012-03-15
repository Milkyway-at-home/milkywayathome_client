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

#include <assert.h>
#include <cmath>
#include <string>
#include <cstdio>

#include "nbody_gl_text.h"
#include "nbody_shaders.h"
#include "nbody_gl_util.h"

glm::mat4 textCameraToClipMatrix(1.0f);

NBodyTextItem::NBodyTextItem(GLuint textProgramPositionLoc, const texture_font_t* textureFont)
{
    this->font = textureFont;

    glGenVertexArrays(1, &this->vao);

    glGenBuffers(1, &this->verticesBuffer);
    glGenBuffers(1, &this->indicesBuffer);

    glBindVertexArray(this->vao);
    glEnableVertexAttribArray(textProgramPositionLoc);
    glBindBuffer(GL_ARRAY_BUFFER, this->verticesBuffer);
    glVertexAttribPointer(textProgramPositionLoc, 4, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->indicesBuffer);

    glBindVertexArray(0);
}

NBodyTextItem::~NBodyTextItem()
{
    glDeleteVertexArrays(1, &this->vao);
    glDeleteBuffers(1, &this->verticesBuffer);
    glDeleteBuffers(1, &this->indicesBuffer);
}

void NBodyTextItem::drawTextItem()
{
    glBindVertexArray(this->vao);
    glBindBuffer(GL_ARRAY_BUFFER, this->verticesBuffer);
    glDrawElements(GL_TRIANGLES, 6 * this->indices.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}

void NBodyTextItem::clearText()
{
    this->vertices.clear();
    this->indices.clear();
}

static float makefont_texture_glyph_get_kerning(const texture_glyph_t* glyph, wchar_t charcode)
{
    const size_t kerning_count = glyph->kerning_count;

    for (size_t i = 0; i < kerning_count; ++i)
    {
        if (glyph->kerning[i].charcode == charcode)
        {
            return glyph->kerning[i].kerning;
        }
    }

    return 0.0f;
}

void NBodyTextItem::addText(const wchar_t* text, TextPen& pen)
{
    float left = pen.x;
    const size_t n = wcslen(text);

    for (size_t i = 0; i < n; ++i)
    {
        if (text[i] == L'\n')
        {
            pen.x = left;
            pen.y -= (this->font->linegap + this->font->height);
            continue;
        }

        int charIdx = text[i] - 32;
        // We only use ASCII characters in the correct order
        if (charIdx >= 0 && charIdx < 96)
        {
            const texture_glyph_t* glyph = &font->glyphs[charIdx];
            float kerning = 0.0;
            if (i > 0)
            {
                kerning = makefont_texture_glyph_get_kerning(glyph, text[i - 1]);
            }

            pen.x += kerning;

            float x0 = std::floor(pen.x + glyph->offset_x);
            float y0 = std::floor(pen.y + glyph->offset_y);
            float x1 = std::floor(x0 + glyph->width) ;
            float y1 = std::floor(y0 - glyph->height);

            GLuint idx = (GLuint) 4 * this->vertices.size();

            TextVertexIndexArray charIndices =
                {
                    {
                        idx, idx + 1, idx + 2,
                        idx, idx + 2, idx + 3
                    }
                };

            TextVertexArray charVertices =
                {
                    {
                        { x0, y0,  glyph->s0, glyph->t0 },
                        { x0, y1,  glyph->s0, glyph->t1 },
                        { x1, y1,  glyph->s1, glyph->t1 },
                        { x1, y0,  glyph->s1, glyph->t0 }
                    }
                };

            this->indices.push_back(charIndices);
            this->vertices.push_back(charVertices);

            pen.x += glyph->advance_x;
        }
    }
}

void NBodyTextItem::uploadText()
{
    glBindBuffer(GL_ARRAY_BUFFER, this->verticesBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(TextVertexArray) * this->vertices.size(), this->vertices.data(), GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->indicesBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(TextVertexIndexArray) * this->indices.size(), this->indices.data(), GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


void NBodyText::loadShader()
{
    this->textProgram.program = nbglCreateProgram("text program",
                                                  (const char*) text_vertex_glsl,
                                                  (const char*) text_fragment_glsl,
                                                  (GLint) text_vertex_glsl_len,
                                                  (GLint) text_fragment_glsl_len);

    this->textProgram.positionLoc = glGetAttribLocation(this->textProgram.program, "position");
    this->textProgram.textTextureLoc = glGetUniformLocation(this->textProgram.program, "textTexture");
    this->textProgram.modelToCameraMatrixLoc = glGetUniformLocation(this->textProgram.program, "modelToCameraMatrix");
    this->textProgram.cameraToClipMatrixLoc = glGetUniformLocation(this->textProgram.program, "cameraToClipMatrix");
}

NBodyText::NBodyText()
{
    this->loadShader();
    this->loadTextTexture();

    this->varText = new NBodyTextItem(this->textProgram.positionLoc, &ubuntuFontR);
    this->constText = new NBodyTextItem(this->textProgram.positionLoc, &ubuntuFontR);
}

NBodyText::~NBodyText()
{
    if (this->textProgram.program != 0)
        glDeleteProgram(this->textProgram.program);

    delete this->varText;
    delete this->constText;

    glDeleteTextures(1, &this->textTexture);
}

void NBodyText::drawProgressText(const SceneData& sceneData)
{
    wchar_t buf[1024];

    /* Start right after the constant portion */
    TextPen pen = this->penEndConst;

    swprintf(buf, sizeof(buf),
             L"Time: %4.3f / %4.3f Gyr (%4.3f %%)\n",
             sceneData.currentTime,
             sceneData.timeEvolve,
             100.0f * sceneData.currentTime / sceneData.timeEvolve
        );

    this->varText->clearText();
    this->varText->addText(buf, pen);
    this->varText->uploadText();

    // Fix black boxes appearing behind letters
    // also keep text looking same when things are behind it
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(this->textProgram.program);
    glUniformMatrix4fv(this->textProgram.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(textCameraToClipMatrix));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, this->textTexture);
    glUniform1i(this->textProgram.textTextureLoc, 0);
    this->constText->drawTextItem();
    this->varText->drawTextItem();

    glUseProgram(0);
    glBindVertexArray(0);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void NBodyText::prepareConstantText(const scene_t* scene)
{
    TextPen pen;
    wchar_t buf[1024];

    swprintf(buf, sizeof(buf),
             L"N-body simulation (%d particles)\n",
             scene->nbody);

    this->constText->addText(buf, pen);
    this->constText->uploadText();

    this->penEndConst = pen;
}

void NBodyText::loadTextTexture()
{
    const texture_font_t* font = &ubuntuFontR;

    glGenTextures(1, &this->textTexture);
    glBindTexture(GL_TEXTURE_2D, this->textTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, font->tex_width, font->tex_height,
                 0, GL_RED, GL_UNSIGNED_BYTE, font->tex_data);
    glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, this->textTexture);
}


