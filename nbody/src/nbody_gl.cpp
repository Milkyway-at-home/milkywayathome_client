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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#pragma GCC diagnostic pop


#include "nbody_graphics.h"
#include "nbody_gl.h"
#include "nbody_config.h"
#include "milkyway_util.h"


#include <assert.h>

#include <errno.h>
#include <cstring>
#include <cstdlib>

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <ios>
#include <stdexcept>
#include <stack>

#if !BOINC_APPLICATION
  #include <sys/mman.h>
  #include <sys/stat.h>
  #include <fcntl.h>
  #include <errno.h>
#endif /* !BOINC_APPLICATION */

#define GLFW_INCLUDE_GL3 1
//#define GLFW_NO_GLU 1
#include <GL/glfw3.h>

//#include <glload/gl_3_3.h>
#include <glutil/glutil.h>


#include "galaxy_model.h"

#define AXES_LENGTH 10.0f


static const float zNear = 0.01f;
static const float zFar = 1000.0f;

static glm::mat4 cameraToClipMatrix;

static glutil::ViewData initialViewData =
{
	glm::vec3(0.0f, 0.0f, 0.0f),  // center position
    glm::fquat(1.0f, 0.0f, 0.0f, 0.0f), // view direction
	30.0f,  // radius
	0.0f    // spin rotation of up axis
};

static glutil::ViewScale viewScale =
{
	0.05f, 100.0f, // min, max view radius
	0.5f, 0.1f,   // radius deltas
    4.0f, 1.0f,
	90.0f / 250.0f // rotation scale
};

static glutil::ViewPole viewPole = glutil::ViewPole(initialViewData, viewScale, glutil::MB_LEFT_BTN);

static void errorHandler(int errorCode, const char* msg)
{
    std::cerr << "Error: " << msg << std::endl;
}

static void resizeHandler(GLFWwindow window, int w, int h)
{
    std::cout << "Resize " << w << " " << h << std::endl;
    glutil::MatrixStack persMatrix;
    persMatrix.Perspective(90.0f, (float) w / (float) h, zNear, zFar);
    cameraToClipMatrix = persMatrix.Top();

    //cameraToClipMatrix = glm::perspective(90.0f, 1.0f, zNear, zFar);

    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

static int closeHandler(GLFWwindow window)
{
    std::cout << "Close " << std::endl;
    return 0;
}

static void refreshHandler(GLFWwindow window)
{
    std::cout << "Refresh" << std::endl;
}

static void focusHandler(GLFWwindow window, int focus)
{
    std::cout << "Focus: " << focus << std::endl;
}

static int getGLFWModifiers(GLFWwindow window)
{
    int modifiers = 0;

    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
    {
        modifiers |= glutil::MM_KEY_SHIFT;
    }

    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
    {
        modifiers |= glutil::MM_KEY_CTRL;
    }

    return modifiers;
}

static glutil::MouseButtons glfwButtonToGLUtil(int button)
{
    switch (button)
    {
        case GLFW_MOUSE_BUTTON_LEFT:
            return glutil::MB_LEFT_BTN;
        case GLFW_MOUSE_BUTTON_RIGHT:
            return glutil::MB_RIGHT_BTN;
        case GLFW_MOUSE_BUTTON_MIDDLE:
        default:
            return glutil::MB_MIDDLE_BTN;
    }
}

static void mouseButtonHandler(GLFWwindow window, int button, int action)
{
    int x, y;
    int modifiers = getGLFWModifiers(window);
    glfwGetMousePos(window, &x, &y);
    viewPole.MouseClick(glfwButtonToGLUtil(button), action == GLFW_PRESS, modifiers, glm::ivec2(x, y));
}

static void mousePosHandler(GLFWwindow window, int x, int y)
{
    viewPole.MouseMove(glm::ivec2(x, y));
}

static void scrollHandler(GLFWwindow window, int x, int y)
{
    int direction = y > 0;
    int modifiers = getGLFWModifiers(window);

    viewPole.MouseWheel(direction, modifiers, glm::ivec2(x, y));
}

static void keyHandler(GLFWwindow window, int x, int y)
{
    std::cout << "Key press " << x << " " << y << std::endl;
}

static void nbglSetHandlers()
{
    glfwSetErrorCallback(errorHandler);

    glfwSetWindowSizeCallback(resizeHandler);
    glfwSetWindowCloseCallback(closeHandler);
    glfwSetWindowRefreshCallback(refreshHandler);
    glfwSetWindowFocusCallback(focusHandler);

    glfwSetKeyCallback(keyHandler);
    //glfwSetKeyCallback(charHandler);
    glfwSetMouseButtonCallback(mouseButtonHandler);
    glfwSetMousePosCallback(mousePosHandler);
    glfwSetScrollCallback(scrollHandler);
}

struct Color
{
    GLfloat r, g, b;
};

class NBodyGraphics
{
private:
    const scene_t* scene;
    GLFWwindow window;

    GLuint vao;

    /*
    struct ParticleBuffers
    {
        GLuint position;
        GLuint velocity;
        GLuint acceleration;
        GLuint color;
    } particleBuffers;

    struct GalaxyModelBuffers
    {
        GLuint mesh;
        GLuint colors;
    } galaxyModelBuffers;
    */

    struct MainProgramData
    {
        GLuint program;

        GLint positionLoc;
        GLint colorLoc;

        GLint modelToCameraMatrixLoc;
        GLint cameraToClipMatrixLoc;
    } mainProgram;

    struct GalaxyProgramData
    {
        GLuint program;

        GLint positionLoc;
        GLint modelToCameraMatrixLoc;
        GLint cameraToClipMatrixLoc;
    } galaxyProgram;



    GLuint positionBuffer;
    GLuint velocityBuffer;
    GLuint accelerationBuffer;
    GLuint colorBuffer;

    GLuint axesBuffer;
    GLuint axesColorBuffer;

    GLuint galaxyModelBuffer;


    bool running;

    struct DrawOptions
    {
        bool fullscreen;
        bool screensaverMode;
        bool pause;
        bool drawOrbitTrace;
        bool drawInfo;
        bool drawAxes;
        bool drawParticles;
        bool floatMode;
        bool cmCentered;
        bool drawHelp;
        bool monochromatic;
    } drawOptions;

    GLuint createShaderFromSrc(const char* shaderSrc, GLint len, GLenum type);
    void loadShaders();
    void createBuffers();
    void createPositionBuffer();
    void createAxesBuffers();
    void createGalaxyBuffer();

public:
    NBodyGraphics(const scene_t* scene);
    ~NBodyGraphics();

    void prepareWindow(bool fullscreen);
    void prepareContext();
    void populateBuffers();

    void loadColors();
    void drawAxes();
    void drawGalaxy();
    void drawBodies();

    void mainLoop();
};

NBodyGraphics::NBodyGraphics(const scene_t* scene)
{
    this->window = NULL;
    this->scene = scene;
    this->vao = 0;

    this->positionBuffer = 0;
    this->velocityBuffer = 0;
    this->accelerationBuffer = 0;
    this->colorBuffer = 0;
    this->axesBuffer = 0;
    this->axesColorBuffer = 0;
    this->galaxyModelBuffer = 0;

    this->mainProgram.program = 0;
    this->mainProgram.positionLoc = -1;
    this->mainProgram.colorLoc = -1;
    this->mainProgram.modelToCameraMatrixLoc = -1;
    this->mainProgram.cameraToClipMatrixLoc = -1;

    this->galaxyProgram.program = 0;
    this->galaxyProgram.positionLoc = -1;
    this->galaxyProgram.modelToCameraMatrixLoc = -1;
    this->galaxyProgram.cameraToClipMatrixLoc = -1;


    this->running = false;

    this->drawOptions.fullscreen = false;
    this->drawOptions.screensaverMode = false;
    this->drawOptions.pause = false;
    this->drawOptions.drawOrbitTrace = false;
    this->drawOptions.drawInfo = false;
    this->drawOptions.drawAxes = true;
    this->drawOptions.drawParticles = true;
    this->drawOptions.floatMode = false;
    this->drawOptions.cmCentered = false;
    this->drawOptions.drawHelp = false;
    this->drawOptions.monochromatic = false;
}

NBodyGraphics::~NBodyGraphics()
{
    GLuint buffers[7];

    glDeleteProgram(this->mainProgram.program);
    glDeleteProgram(this->galaxyProgram.program);

    buffers[0] = this->positionBuffer;
    buffers[1] = this->velocityBuffer;
    buffers[2] = this->accelerationBuffer;
    buffers[3] = this->colorBuffer;
    buffers[4] = this->axesBuffer;
    buffers[5] = this->axesColorBuffer;
    buffers[6] = this->galaxyModelBuffer;

    glDeleteBuffers(7, buffers);
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

void NBodyGraphics::createGalaxyBuffer()
{
    glBindBuffer(GL_ARRAY_BUFFER, this->galaxyModelBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(galaxy), (const GLfloat*) galaxy, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

static GLuint createShaderFromSrc(const char* src, GLint len, GLenum type)
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
                "Shader compile log '%s':\n"
                "--------------------------------------------------------------------------------\n"
                "%s\n"
                "--------------------------------------------------------------------------------\n",
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

static void nbglGetProgramLog(GLuint program, const char* name)
{
    GLint logLength = 0;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 0)
    {
        GLchar* log = new GLchar[logLength + 1];
        glGetShaderInfoLog(program, logLength, NULL, log);
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

extern "C" unsigned char vertex_glsl[];
extern "C" size_t vertex_glsl_len;

extern "C" unsigned char fragment_glsl[];
extern "C" size_t fragment_glsl_len;


extern "C" unsigned char galaxy_vertex_glsl[];
extern "C" size_t galaxy_vertex_glsl_len;

extern "C" unsigned char galaxy_fragment_glsl[];
extern "C" size_t galaxy_fragment_glsl_len;


static GLuint nbglCreateProgram(const char* name,
                                const char* vertSrc,
                                const char* fragSrc,
                                GLint vertSrcLen,
                                GLint fragSrcLen)
{
    GLuint vertexShader = createShaderFromSrc(vertSrc, (GLint) vertSrcLen, GL_VERTEX_SHADER);
    GLuint pixelShader = createShaderFromSrc(fragSrc, (GLint) fragSrcLen, GL_FRAGMENT_SHADER);

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

void NBodyGraphics::loadShaders()
{
    this->mainProgram.program = nbglCreateProgram("main program",
                                                  (const char*) vertex_glsl,
                                                  (const char*) fragment_glsl,
                                                  (GLint) vertex_glsl_len,
                                                  (GLint) fragment_glsl_len);

    this->galaxyProgram.program = nbglCreateProgram("galaxy program",
                                                    (const char*) galaxy_vertex_glsl,
                                                    (const char*) galaxy_fragment_glsl,
                                                    (GLint) galaxy_vertex_glsl_len,
                                                    (GLint) galaxy_fragment_glsl_len);
}

void NBodyGraphics::createAxesBuffers()
{
    static const GLfloat axes[6][3] =
        {
            { 0.0f,        0.0f,        0.0f        },
            { AXES_LENGTH, 0.0f,        0.0f        },
            { 0.0f,        0.0f,        0.0f        },
            { 0.0f,        AXES_LENGTH, 0.0f        },
            { 0.0f,        0.0f,        0.0f        },
            { 0.0f,        0.0f,        AXES_LENGTH }
        };

    static const GLfloat colors[6][4] =
        {
            { 1.0f, 0.0f, 0.0f, 0.5f },
            { 1.0f, 0.0f, 0.0f, 0.5f },
            { 0.0f, 1.0f, 0.0f, 0.5f },
            { 0.0f, 1.0f, 0.0f, 0.5f },
            { 0.0f, 0.0f, 1.0f, 0.5f },
            { 0.0f, 0.0f, 1.0f, 0.5f }
        };

    glBindBuffer(GL_ARRAY_BUFFER, this->axesBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(axes), (const GLfloat*) axes, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, this->axesColorBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(colors), (const GLfloat*) colors, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void NBodyGraphics::drawAxes()
{
    glBindBuffer(GL_ARRAY_BUFFER, this->axesBuffer);
    glVertexAttribPointer(this->mainProgram.positionLoc, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, this->axesColorBuffer);
    glVertexAttribPointer(this->mainProgram.colorLoc, 4, GL_FLOAT, GL_FALSE, 0, 0);

    glDrawArrays(GL_LINES, 0, 6);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void NBodyGraphics::drawGalaxy()
{
    glBindBuffer(GL_ARRAY_BUFFER, this->galaxyModelBuffer);
    glVertexAttribPointer(this->galaxyProgram.positionLoc, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);
    //glVertexPointer(3, GL_FLOAT, sizeof(vertexData), &galaxy[0].vertex);
    //glNormalPointer(GL_FLOAT, sizeof(vertexData), &galaxy[0].normal);

    glDrawArrays(GL_TRIANGLES, 0, sizeof(galaxy) / (6 * sizeof(GLfloat)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void NBodyGraphics::drawBodies()
{
    const GLfloat* positions = (const GLfloat*) this->scene->rTrace;
    GLint nbody = this->scene->nbody;

	glBindBuffer(GL_ARRAY_BUFFER, this->positionBuffer);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 4 * nbody * sizeof(GLfloat), positions);

    /* 4th component is not included */
    glVertexAttribPointer(this->mainProgram.positionLoc, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, this->colorBuffer);
    glVertexAttribPointer(this->mainProgram.colorLoc, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glDrawArrays(GL_POINTS, 0, nbody);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void NBodyGraphics::createBuffers()
{
    GLuint buffers[7];

    glGenVertexArrays(1, &this->vao);
    glBindVertexArray(this->vao);

    glGenBuffers(7, buffers);

    this->positionBuffer = buffers[0];
    this->velocityBuffer = buffers[1];
    this->accelerationBuffer = buffers[2];
    this->colorBuffer = buffers[3];
    this->axesBuffer = buffers[4];
    this->axesColorBuffer = buffers[5];
    this->galaxyModelBuffer = buffers[6];
}

void NBodyGraphics::createPositionBuffer()
{
    GLint nbody = this->scene->nbody;
    const GLfloat* positions = (const GLfloat*) this->scene->rTrace;

    glBindBuffer(GL_ARRAY_BUFFER, this->positionBuffer);
	glBufferData(GL_ARRAY_BUFFER, 4 * nbody * sizeof(GLfloat), positions, GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void NBodyGraphics::populateBuffers()
{
    this->createPositionBuffer();
    this->createAxesBuffers();
    this->createGalaxyBuffer();
}

void NBodyGraphics::loadColors()
{
    GLint nbody = this->scene->nbody;

    /* assign random particle colors */
    srand((unsigned int) time(NULL));

    Color* color = new Color[nbody];

    for (GLint i = 0; i < nbody; ++i)
    {

        /*
        if (r[i].ignore)
        {
            R = grey.x;  // TODO: Random greyish color?
            G = grey.z;
            B = grey.y;
        }
        else
        {
            R = ((double) rand()) / ((double) RAND_MAX);
            G = ((double) rand()) / ((double) RAND_MAX) * (1.0 - R);
            B = 1.0 - R - G;
        }
        */

        double R = ((double) rand()) / ((double) RAND_MAX);
        double G = ((double) rand()) / ((double) RAND_MAX) * (1.0 - R);
        double B = 1.0 - R - G;

        double scale;
        if (R >= G && R >= B)
        {
            scale = 1.0 + ((double) rand()) / ((double) RAND_MAX) * (std::min(2.0, 1.0 / R) - 1.0);
        }
        else if (G >= R && G >= B)
        {
            scale = 1.0 + ((double) rand()) / ((double) RAND_MAX) * (std::min(2.0, 1.0 / G) - 1.0);
        }
        else
        {
            scale = 1.0 + ((double) rand()) / ((double) RAND_MAX) * (std::min(2.0, 1.0 / B) - 1.0);
        }

        //color[i].ignore = color[i].ignore;
        color[i].r = (GLfloat) R * scale;
        color[i].g = (GLfloat) G * scale;
        color[i].b = (GLfloat) B * scale;
    }

    glBindBuffer(GL_ARRAY_BUFFER, this->colorBuffer);
	glBufferData(GL_ARRAY_BUFFER, 3 * nbody * sizeof(GLfloat), color, GL_STATIC_DRAW);
    glVertexAttribPointer(this->mainProgram.colorLoc, 3, GL_FLOAT, GL_FALSE, 0, 0);

    delete[] color;
}

void NBodyGraphics::prepareContext()
{
    this->loadShaders();
    this->createBuffers();

    this->mainProgram.positionLoc = glGetAttribLocation(this->mainProgram.program, "position");
    this->mainProgram.colorLoc = glGetAttribLocation(this->mainProgram.program, "inputColor");
    this->mainProgram.modelToCameraMatrixLoc = glGetUniformLocation(this->mainProgram.program, "modelToCameraMatrix");
    this->mainProgram.cameraToClipMatrixLoc = glGetUniformLocation(this->mainProgram.program, "cameraToClipMatrix");

    glEnableVertexAttribArray(this->mainProgram.positionLoc);
    glEnableVertexAttribArray(this->mainProgram.colorLoc);



    this->galaxyProgram.positionLoc = glGetAttribLocation(this->galaxyProgram.program, "position");
    this->galaxyProgram.modelToCameraMatrixLoc = glGetUniformLocation(this->galaxyProgram.program, "modelToCameraMatrix");
    this->galaxyProgram.cameraToClipMatrixLoc = glGetUniformLocation(this->galaxyProgram.program, "cameraToClipMatrix");
    glEnableVertexAttribArray(this->galaxyProgram.positionLoc);


    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth(1.0);

    glDepthFunc(GL_LESS);
    //glEnable(GL_DEPTH_TEST | GL_BLEND | GL_ALPHA_TEST | GL_POINT_SMOOTH | GL_VERTEX_PROGRAM_POINT_SIZE);
    //glEnable(GL_DEPTH_TEST | GL_BLEND | GL_ALPHA_TEST | GL_POINT_SMOOTH);
    glEnable(GL_DEPTH_TEST | GL_BLEND | GL_ALPHA_TEST);

    glPointSize(3.0f);

    // allow changing point size from within shader
    // as well as smoothing them to look more spherical
    //glEnable(GL_POINT_SMOOTH | GL_VERTEX_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SMOOTH);
    //glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SPRITE);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glEnable(GL_POINT_SPRITE_ARB);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CW);

    float maxSmoothPointSize[2];
    glGetFloatv(GL_SMOOTH_POINT_SIZE_RANGE, (float *)&maxSmoothPointSize);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glShadeModel(GL_SMOOTH);

    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDepthFunc(GL_LEQUAL);
    glDepthRange(0.0f, 1.0f);
}

static void requestGL32()
{
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3);
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 2);
    glfwOpenWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwOpenWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);

  #ifndef NDEBUG
    glfwOpenWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
  #endif
}

void NBodyGraphics::prepareWindow(bool fullscreen)
{
    GLFWvidmode vidMode;
    glfwGetDesktopMode(&vidMode);
    int winMode = fullscreen ? GLFW_FULLSCREEN : GLFW_WINDOWED;

    int width = vidMode.width / 2;
    int height = vidMode.height / 2;

    requestGL32();
    this->window = glfwOpenWindow(width, height, winMode, "Milkyway@Home N-body", NULL);
    if (!this->window)
    {
        throw std::runtime_error("Failed to open window");
    }
}

void NBodyGraphics::mainLoop()
{
    this->running = true;

    while (this->running)
    {
        glfwPollEvents();

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClearDepth(1.0);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        glUseProgram(this->mainProgram.program);

        glutil::MatrixStack modelMatrix;
        modelMatrix.SetMatrix(viewPole.CalcMatrix());

        //std::cout << viewPole.GetView().orient[0] << std::endl;
		//modelMatrix.Scale(glm::vec3(0.05f, 0.05, 0.05f));
		//modelMatrix.Translate(glm::vec3(0.0f, 0.5f, 0.0f));

        glUniformMatrix4fv(this->mainProgram.modelToCameraMatrixLoc, 1, GL_FALSE, glm::value_ptr(modelMatrix.Top()));
        glUniformMatrix4fv(this->mainProgram.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(cameraToClipMatrix));

        if (this->drawOptions.drawAxes)
        {
            this->drawAxes();
        }

        if (this->drawOptions.drawParticles)
        {
            this->drawBodies();
        }

        if (true)
        {
            glUseProgram(this->galaxyProgram.program);

            glUniformMatrix4fv(this->galaxyProgram.modelToCameraMatrixLoc, 1, GL_FALSE, glm::value_ptr(modelMatrix.Top()));
            glUniformMatrix4fv(this->galaxyProgram.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(cameraToClipMatrix));


            this->drawGalaxy();
        }

        glUseProgram(0);
        glfwSwapBuffers();
    }
}

static void printVersionAndExts()
{
    std::cout << "GL version: "      << glGetString(GL_VERSION)
              << " Shader version: " << glGetString(GL_SHADING_LANGUAGE_VERSION)
              << std::endl;

    GLint n;
    glGetIntegerv(GL_NUM_EXTENSIONS, &n);
    for (GLint i = 0; i < n; i++)
    {
        printf("%s\n", glGetStringi(GL_EXTENSIONS, i));
    }
}

#if !BOINC_APPLICATION

scene_t* nbConnectSharedScene(int instanceId)
{
    int shmId;
    const int mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
    struct stat sb;
    char name[128];
    scene_t* scene = NULL;

    if (snprintf(name, sizeof(name), "/milkyway_nbody_%d", instanceId) == sizeof(name))
    {
        mw_panic("name buffer too small for shared memory name\n");
    }

    shmId = shm_open(name, O_RDWR, mode);
    if (shmId < 0)
    {
        mwPerror("Error getting shared memory");
        return NULL;
    }

    if (fstat(shmId, &sb) < 0)
    {
        mwPerror("shmem fstat");
        shm_unlink(name);
        return NULL;
    }

    scene = (scene_t*) mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, shmId, 0);
    if (scene == MAP_FAILED)
    {
        mwPerror("mmap: Failed to mmap shared memory");
        if (shm_unlink(name) < 0)
        {
            mwPerror("Unlink shared memory");
        }

        return NULL;
    }

    return scene;
}

#else

/* Returns TRUE if connection succeeds */
static int nbAttemptConnectSharedScene(void)
{
    scene_t* scene = (scene_t*) mw_graphics_get_shmem(NBODY_BIN_NAME);
    if (!scene)
    {
        mw_printf("Failed to connect to shared scene\n");
    }

    return scene;
}

#define MAX_TRIES 5
#define RETRY_INTERVAL 250

/* In case the main application isn't ready yet, try and wait for a while */
scene_t* nbConnectSharedScene(int instanceId)
{
    int tries = 0;

    while (tries < MAX_TRIES)
    {
        if ((scene = nbAttemptConnectSharedScene()))
        {
            return alreadyAttached(); /* Error if something already attached */
        }

        mwMilliSleep(RETRY_INTERVAL);
        ++tries;
    }

    mw_printf("Could not attach to simulation after %d attempts\n", MAX_TRIES);
    return NULL;
}

#endif /* !BOINC_APPLICATION */

int nbCheckConnectedVersion(const scene_t* scene)
{
    if (   scene->nbodyMajorVersion != NBODY_VERSION_MAJOR
        || scene->nbodyMinorVersion != NBODY_VERSION_MINOR)
    {
        mw_printf("Graphics version (%d.%d) does not match application version (%d.%d)\n",
                  NBODY_VERSION_MAJOR,
                  NBODY_VERSION_MINOR,
                  scene->nbodyMajorVersion,
                  scene->nbodyMinorVersion);
        return 1;
    }

    return 0;
}

int nbRunGraphics(const scene_t* scene, const VisArgs* args)
{
    NBodyGraphics graphicsContext(scene);
    int rc = 0;

    glfwInit();

    try
    {
        graphicsContext.prepareWindow(false);
        graphicsContext.prepareContext();
        nbglSetHandlers();
        graphicsContext.populateBuffers();
        graphicsContext.loadColors();

        graphicsContext.mainLoop();
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        rc = 1;
    }

    glfwTerminate();

    return rc;
}

