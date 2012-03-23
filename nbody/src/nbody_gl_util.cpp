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

