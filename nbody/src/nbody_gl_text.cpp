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
#include "nbody_gl_shaders.h"
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
    glDrawElements(GL_TRIANGLES, 6 * (GLsizei) this->indices.size(), GL_UNSIGNED_INT, 0);
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

// returns the farthest right and farthest down the added text reached
TextPen NBodyTextItem::addText(const wchar_t* text, TextPen& pen)
{
    float left = pen.x;
    const size_t n = wcslen(text);

    TextPen limits = pen;

    for (size_t i = 0; i < n; ++i)
    {
        if (text[i] == L'\n')
        {
            pen.x = left;
            pen.y -= (this->font->linegap + this->font->height);
            continue;
        }

        int charIdx = text[i] - 32 + 1;
        // We only use ASCII characters in the correct order
        // + 1 for -1 character at start
        if (charIdx > 0 && charIdx < 96 + 1)
        {
            const texture_glyph_t* glyph = &this->font->glyphs[charIdx];
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

            GLuint idx = 4 * (GLuint) this->vertices.size();

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

            // text grows right
            limits.x = std::max(pen.x, limits.x);

            // text grows down
            limits.y = std::min(pen.y, limits.y);
        }
    }

    return limits;
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

NBodyText::NBodyText(const texture_font_t* textFont)
{
    this->font = textFont;

    this->loadShader();
    this->loadTextTexture();

    this->varText = new NBodyTextItem(this->textProgram.positionLoc, this->font);
    this->constText = new NBodyTextItem(this->textProgram.positionLoc, this->font);
    this->helpTextLeftColumn = new NBodyTextItem(this->textProgram.positionLoc, this->font);
    this->helpTextRightColumn = new NBodyTextItem(this->textProgram.positionLoc, this->font);

    this->prepareHelpText();
}

NBodyText::~NBodyText()
{
    glDeleteProgram(this->textProgram.program);

    delete this->varText;
    delete this->constText;
    delete this->helpTextLeftColumn;
    delete this->helpTextRightColumn;

    glDeleteTextures(1, &this->textTexture);
}

void NBodyText::useTextProgram()
{
    // Fix black boxes appearing behind letters
    // also keep text looking same when things are behind it
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(this->textProgram.program);
    glUniformMatrix4fv(this->textProgram.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(textCameraToClipMatrix));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, this->textTexture);
    glUniform1i(this->textProgram.textTextureLoc, 0);
}

void NBodyText::resetTextProgram()
{
    glUseProgram(0);
    glBindVertexArray(0);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void NBodyText::drawProgressText(const SceneData& sceneData)
{
    this->useTextProgram();

    if (sceneData.staticScene)
    {
        this->constText->drawTextItem();
    }
    else
    {
        wchar_t buf[256];

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

        this->constText->drawTextItem();
        this->varText->drawTextItem();
    }

    this->resetTextProgram();
}


void NBodyText::drawHelpText()
{
    this->useTextProgram();
    this->helpTextLeftColumn->drawTextItem();
    this->helpTextRightColumn->drawTextItem();
    this->resetTextProgram();
}

// we split the help text into 2 pieces so that we have key,
// description columns line up without using a non-monospace font

// key and maybe mnemonic column
static const wchar_t* nbglHelpTextLeft =
    L"\n"
    L"Mouse controls:\n"
    L"  Drag\n"
    L"  [Alt]  Drag\n"
    L"  [Ctrl] Drag\n"
    L"    + [Shift]\n"
    L"  Scroll\n"
    L"  [Shift] Scroll\n"
    L"\n"
    L"Key controls:\n"
    L"  ? / h (help)\n"
    L"  q, Escape\n"
    L"  space / p\n"
    L"  a (axes)\n"
    L"  o (origin)\n"
    L"  f (float)\n"
    L"  > / z\n"
    L"  < / x\n"
    L"  t (trace)\n"
    L"  i (info)\n"
    L"  c (color)\n"
    L"  d (draw mode)\n"
    L"  + / b (bigger)\n"
    L"  - / s (smaller)\n"
    L"  =\n"
    L"\n";

// longer description column
static const wchar_t* nbglHelpTextRight =
    L"\n"
    L"\n"
    L": rotate on X and Y (view local) axes\n"
    L": rotate in view local Z direction\n"
    L": rotate X and Y axes based on mouse distance from start\n"
    L": zoom functions with smaller steps \n"
    L": zoom\n"
    L": zoom with smaller steps\n"
    L"\n"
    L"\n"
    L": display this help text\n"
    L": quit\n"
    L": pause/unpause\n"
    L": toggle displaying axes\n"
    L": toggle center between origin or center of mass\n"
    L": float view around randomly\n"
    L": increase random float speed\n"
    L": decrease random float speed\n"
    L": toggle displaying orbit trace of center of mass\n"
    L": toggle information display\n"
    L": toggle particle color scheme (all white vs. colored)\n"
    L": toggle using textured (prettier) points and plain (slightly faster)\n"
    L": make points bigger\n"
    L": (smaller) make points smaller\n"
    L": reset points to initial size\n"
    L"\n";

void NBodyText::prepareHelpText()
{
    TextPen pen(0.0f, 0.0f);

    TextPen limit = this->helpTextLeftColumn->addText(nbglHelpTextLeft, pen);

    this->helpTextLeftColumn->uploadText();

    // move the pen up past the column
    pen.x = limit.x + 0.5f * this->font->size;
    pen.y = 0.0f;

    this->helpTextRightColumn->addText(nbglHelpTextRight, pen);
    this->helpTextRightColumn->uploadText();
}

void NBodyText::prepareConstantText(const scene_t* scene)
{
    TextPen pen;
    wchar_t buf[256];

    if (scene->staticScene)
    {
        swprintf(buf, sizeof(buf),
                 L"Static N-body scene (%d particles)\n",
                 scene->nbody);
    }
    else
    {
        swprintf(buf, sizeof(buf),
                 L"N-body simulation (%d particles)\n",
                 scene->nbody);
    }

    this->constText->addText(buf, pen);
    this->constText->uploadText();

    this->penEndConst = pen;
}

void NBodyText::loadTextTexture()
{
    glGenTextures(1, &this->textTexture);
    glBindTexture(GL_TEXTURE_2D, this->textTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, (GLsizei) font->tex_width, (GLsizei) font->tex_height,
                 0, GL_RED, GL_UNSIGNED_BYTE, font->tex_data);
    glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, this->textTexture);
}


