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

#include "nbody_gl_util.h"
#include "nbody_gl_includes.h"

#include <cstdlib>
#include <cstdio>
#include <stdexcept>

static void nbglGetProgramLog(GLuint program, const char* name)
{
    GLint logLength = 0;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 0)
    {
        GLchar* log = new GLchar[logLength + 1];
        glGetProgramInfoLog(program, logLength, NULL, log);
        fprintf(stderr,
                "Linker output (%s):\n"
                "--------------------------------------------------------------------------------\n"
                "%s\n"
                "--------------------------------------------------------------------------------\n",
                name,
                log);

        delete[] log;
    }

	GLint status = GL_FALSE;
	glGetProgramiv(program, GL_LINK_STATUS, &status);
	if (!status)
	{
        throw std::runtime_error("Linking shader program failed");
	}
}

static const char* showShaderType(GLenum type)
{
    switch (type)
    {
        case GL_VERTEX_SHADER:
            return "Vertex shader";
        case GL_FRAGMENT_SHADER:
            return "Fragment shader";
        default:
            return "<bad shader type>";
    }
}

static GLuint nbglCreateShaderFromSrc(const char* name, const char* src, GLint len, GLenum type)
{
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &src, &len);
    glCompileShader(shader);

    GLint logLength = 0;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 0)
    {
        GLchar* logBuf = new GLchar[logLength + 1];
        glGetShaderInfoLog(shader, logLength, NULL, logBuf);

        fprintf(stderr,
                "Shader compile log '%s' (%s):\n"
                "--------------------------------------------------------------------------------\n"
                "%s\n"
                "--------------------------------------------------------------------------------\n",
                name,
                showShaderType(type),
                logBuf);

        delete[] logBuf;
    }

    GLint status = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (!status)
    {
        throw std::runtime_error("Error compiling shader");
    }

    return shader;
}

GLuint nbglCreateProgram(const char* name,
                         const char* vertSrc,
                         const char* fragSrc,
                         GLint vertSrcLen,
                         GLint fragSrcLen)
{
    GLuint vertexShader = nbglCreateShaderFromSrc(name, vertSrc, (GLint) vertSrcLen, GL_VERTEX_SHADER);
    GLuint pixelShader = nbglCreateShaderFromSrc(name, fragSrc, (GLint) fragSrcLen, GL_FRAGMENT_SHADER);

    GLuint program = glCreateProgram();

    glAttachShader(program, vertexShader);
    glAttachShader(program, pixelShader);

    glLinkProgram(program);

    glDetachShader(program, vertexShader);
    glDetachShader(program, pixelShader);

    glDeleteShader(vertexShader);
    glDeleteShader(pixelShader);

    nbglGetProgramLog(program, name);

    return program;
}

void nbglPrintGLVersionAndExts()
{
    printf("GL version:     %s\n"
           "Shader version: %s\n"
           "Extensions:\n",
           glGetString(GL_VERSION),
           glGetString(GL_SHADING_LANGUAGE_VERSION));

    GLint n;
    glGetIntegerv(GL_NUM_EXTENSIONS, &n);
    for (GLint i = 0; i < n; i++)
    {
        printf("  %s\n", glGetStringi(GL_EXTENSIONS, i));
    }
}

struct RGBColor
{
    GLfloat r, g, b;
    RGBColor() : r(0.0f), g(0.0f), b(0.0f) { }
    RGBColor(GLfloat rr, GLfloat gg, GLfloat bb) : r(rr), g(gg), b(bb) { }
};

struct HSVColor
{
    GLfloat hue;        /* Hue degree between 0 and 255 */
    GLfloat sat;        /* Saturation between 0 (gray) and 255 */
    GLfloat val;        /* Value between 0 (black) and 255 */

    HSVColor() : hue(0.0f), sat(0.0f), val(0.0f) { }
    HSVColor(GLfloat h, GLfloat s, GLfloat v) : hue(h), sat(s), val(v) { }
};

static void convertHSVToRGB(RGBColor& rgb, const HSVColor& hsv)
{
    if (hsv.sat == 0)
    {
        rgb.r = rgb.g = rgb.b = hsv.val;
        return;
    }

    GLfloat h = fmod(hsv.hue, 360.0f) / 60.0f;

    int region = (int) h;
    GLfloat f = h - (GLfloat) region;
    GLfloat p = hsv.val * (1.0f - hsv.sat);
    GLfloat q = hsv.val * (1.0f - (hsv.sat * f));
    GLfloat t = hsv.val * (1.0f - (hsv.sat * (1.0f - f)));

    switch (region)
    {
        case 0:
            rgb.r = hsv.val;
            rgb.g = t;
            rgb.b = p;
            break;

        case 1:
            rgb.r = q;
            rgb.g = hsv.val;
            rgb.b = p;
            break;

        case 2:
            rgb.r = p;
            rgb.g = hsv.val;
            rgb.b = t;
            break;

        case 3:
            rgb.r = p;
            rgb.g = q;
            rgb.b = hsv.val;
            break;

        case 4:
            rgb.r = t;
            rgb.g = p;
            rgb.b = hsv.val;
            break;

        case 5:
        default:
            rgb.r = hsv.val;
            rgb.g = p;
            rgb.b = q;
            break;
    }
}

static void nbglBetweenHSVColors(RGBColor& rgbMid, const HSVColor& a, const HSVColor& b)
{
    HSVColor mid(glm::linearRand(a.hue, b.hue),
                 glm::linearRand(a.sat, b.sat),
                 glm::linearRand(a.val, b.val));
    convertHSVToRGB(rgbMid, mid);
}

void nbglRandomParticleColor(Color& color, bool ignore)
{
//    static const HSVColor dark(213.0f, 0.24f, 0.36f);

//    static const HSVColor light(202.0f, 0.13f, 0.78f);
//    static const HSVColor dark(212.0f, 0.32f, 0.21f);

//    static const HSVColor dark(212.0f, 0.28f, 0.36f);
//    static const HSVColor dark(213.0f, 0.32f, 0.22f);

//    static const HSVColor dark(218.0f, 0.25f, 0.25f);
//    static const HSVColor dark(240.0f, 0.9f, 0.21f);

//    static const HSVColor light(206.0f, 0.25f, 76.0f);
//    static const HSVColor light(198.0f, 8.0f, 85.0f);

    static const HSVColor light(204.0f, 0.13f, 0.84f);
    static const HSVColor dark(209.0f, 0.41f, 0.50f);

    static const HSVColor coreLight(40.0f, 0.01f, 0.99f);
    static const HSVColor coreDark(33.0f, 0.29f, 0.58f);

    RGBColor rgb;
    if (glm::linearRand(0.0f, 1.0f) > 0.2f)
    {
        // bluish
        nbglBetweenHSVColors(rgb, dark, light);
    }
    else
    {
        // reddish
        nbglBetweenHSVColors(rgb, coreDark, coreLight);
    }

    color.r = rgb.r;
    color.g = rgb.g;
    color.b = rgb.b;
    color.ignore = 1.0f;
    //color.ignore = ignore ? ?????
}

