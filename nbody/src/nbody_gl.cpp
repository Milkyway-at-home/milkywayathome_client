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

#include "nbody_gl_includes.h"
#include "nbody_graphics.h"
#include "nbody_gl.h"
#include "nbody_config.h"
#include "milkyway_util.h"

#include "nbody_gl_util.h"

#include "nbody_gl_axes.h"
#include "nbody_gl_text.h"
#include "nbody_gl_private.h"

#include <assert.h>

#include <errno.h>
#include <cstring>
#include <cstdlib>

#include <iostream>
#include <cmath>
#include <stdexcept>

#if !BOINC_APPLICATION
  #include <sys/mman.h>
  #include <sys/stat.h>
  #include <fcntl.h>
  #include <errno.h>
#endif /* !BOINC_APPLICATION */

#include "nbody_gl_includes.h"

#include "nbody_particle_texture.h"
#include "nbody_gl_resources.h"
#include "nbody_gl_galaxy_model.h"

static const float zNear = 0.01f;
static const float zFar = 1000.0f;

glm::mat4 cameraToClipMatrix(1.0f);


class NBodyGraphics;
// For access from callbacks
static NBodyGraphics* globalGraphicsContext = NULL;


static glutil::ViewData initialViewData =
{
	glm::vec3(0.0f, 0.0f, 0.0f),  // center position
    glm::fquat(1.0f, 0.0f, 0.0f, 0.0f), // view direction
	30.0f,  // radius
	0.0f    // spin rotation of up axis
};

static glutil::ViewScale viewScale =
{
	0.05f, 250.0f, // min, max view radius
	0.5f, 0.1f,   // radius deltas
    4.0f, 1.0f,
	90.0f / 250.0f // rotation scale
};

static glutil::ViewPole viewPole = glutil::ViewPole(initialViewData, viewScale, glutil::MB_LEFT_BTN);

struct Color
{
    GLfloat r, g, b;
    GLfloat ignore;
};

class NBodyGraphics
{
private:
    scene_t* scene;
    GLFWwindow window;

    GLuint particleVAO;
    GLuint whiteParticleVAO;

    struct ParticleTextureProgramData
    {
        GLuint program;

        GLint positionLoc;
        GLint colorLoc;

        GLint modelToCameraMatrixLoc;
        GLint cameraToClipMatrixLoc;

        GLint particleTextureLoc;
        GLint pointSizeLoc;
    } particleTextureProgram;

    struct ParticlePointProgramData
    {
        GLuint program;

        GLint positionLoc;
        GLint colorLoc;

        GLint modelToCameraMatrixLoc;
        GLint cameraToClipMatrixLoc;

        GLint pointSizeLoc;
    } particlePointProgram;

    GLuint positionBuffer;
    GLuint velocityBuffer;
    GLuint accelerationBuffer;
    GLuint colorBuffer;
    GLuint whiteBuffer;
    GLuint particleTexture;

    bool running;

    GalaxyModel* galaxyModel;

    NBodyText text;
    NBodyAxes axes;

    SceneData sceneData;

    enum DrawMode
    {
        POINTS,
        TEXTURED_SPRITES
    };

    struct DrawOptions
    {
        bool screensaverMode;
        bool paused;
        bool floatMode;

        bool cmCentered;
        bool monochromatic;

        DrawMode drawMode;
        // Keep a separate point size for each draw mode
        float texturedSpritePointSize;
        float pointPointSize;

        bool drawOrbitTrace;
        bool drawInfo;
        bool drawAxes;
        bool drawParticles;
        bool drawHelp;

        DrawOptions(const VisArgs* args)
        {
            this->screensaverMode = (args->fullscreen && !args->plainFullscreen);
            this->paused = false;
            this->floatMode = (this->screensaverMode && !args->noFloat);

            this->cmCentered = !args->originCenter;
            this->monochromatic = (bool) args->monochrome;

            this->drawMode = args->untexturedPoints ? POINTS : TEXTURED_SPRITES;
            this->texturedSpritePointSize = 250.0f;
            this->pointPointSize = 5.0f;

            this->drawOrbitTrace = false;
            this->drawInfo = true;
            this->drawAxes = (bool) args->drawAxes;
            this->drawParticles = true;
            this->drawHelp = false;
        }
    } drawOptions;


    void loadShaders();
    void createBuffers();
    void prepareColoredVAO(GLuint& vao, GLuint color);
    void prepareVAOs();
    void createPositionBuffer();
    void drawParticlesTextured(const glm::mat4& modelMatrix);
    void drawParticlesPoints(const glm::mat4& modelMatrix);

public:
    NBodyGraphics(scene_t* scene, const VisArgs* args);
    ~NBodyGraphics();

    void prepareContext();
    void populateBuffers();

    void loadModel(GalaxyModel& model);
    void loadColors();
    void drawAxes();
    void drawParticles(const glm::mat4& modelMatrix);
    bool readSceneData();

    void display();
    void mainLoop();

    void stop()
    {
        this->running = false;
    }

    void toggleDrawMode()
    {
        if (this->drawOptions.drawMode == POINTS)
        {
            this->drawOptions.drawMode = TEXTURED_SPRITES;
        }
        else
        {
            this->drawOptions.drawMode = POINTS;
        }
    }

    void increasePointSize()
    {
        if (this->drawOptions.drawMode == TEXTURED_SPRITES)
        {
            float size = this->drawOptions.texturedSpritePointSize;
            size *= 1.05f;
            if (size > 1000.0f)
            {
                size = 1000.0f;
            }

            this->drawOptions.texturedSpritePointSize = size;
        }
        else
        {
            float size = this->drawOptions.pointPointSize;
            size *= 1.05f;
            if (size > 200.0f)
            {
                size = 200.0f;
            }

            this->drawOptions.pointPointSize = size;
        }
    }

    void decreasePointSize()
    {
        if (this->drawOptions.drawMode == TEXTURED_SPRITES)
        {
            float size = this->drawOptions.texturedSpritePointSize;
            size /= 1.05f;
            if (size < 1.0e-3f)
            {
                size = 1.0e-3f;
            }

            this->drawOptions.texturedSpritePointSize = size;
        }
        else
        {
            float size = this->drawOptions.pointPointSize;

            size /= 1.05f;
            if (size < 1.0e-3f)
            {
                size = 1.0e-3f;
            }

            this->drawOptions.pointPointSize = size;
        }
    }

    void toggleDrawAxes()
    {
        this->drawOptions.drawAxes = !this->drawOptions.drawAxes;
    }

    void toggleDrawInfo()
    {
        this->drawOptions.drawInfo = !this->drawOptions.drawInfo;
    }

    void toggleMonochromatic()
    {
        this->drawOptions.monochromatic = !this->drawOptions.monochromatic;
    }

    void togglePaused()
    {
        this->drawOptions.paused = !this->drawOptions.paused;
    }

    void toggleFloatMode()
    {
        this->drawOptions.floatMode = !this->drawOptions.floatMode;
    }

    void toggleOrigin()
    {
        this->drawOptions.cmCentered = !this->drawOptions.cmCentered;
    }
};

static void errorHandler(int errorCode, const char* msg)
{
    fprintf(stderr, "GLFW error (%d): %s\n", errorCode, msg);
}

static void resizeHandler(GLFWwindow window, int w, int h)
{
    float wf = (float) w;
    float hf = (float) h;
    float aspectRatio = wf / hf;
    cameraToClipMatrix = glm::perspective(90.0f, aspectRatio, zNear, zFar);

    const float fontHeight = 12.0f; // FIXME: hardcoded font height
    textCameraToClipMatrix = glm::ortho(0.0f, wf, -hf, 0.0f);
    textCameraToClipMatrix = glm::translate(textCameraToClipMatrix, glm::vec3(0.0f, -fontHeight, 0.0f));

    glViewport(0, 0, (GLsizei) w, (GLsizei) h);

    if (globalGraphicsContext)
    {
        globalGraphicsContext->display();
        glfwSwapBuffers();
    }
}

static int closeHandler(GLFWwindow window)
{
    if (globalGraphicsContext)
    {
        globalGraphicsContext->stop();
    }

    return 0;
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

static void keyHandler(GLFWwindow window, int key, int pressed)
{
    if (!pressed) // release
        return;

    NBodyGraphics* ctx = globalGraphicsContext;
    if (!ctx)
        return;

    switch (key)
    {
        case GLFW_KEY_A:
            ctx->toggleDrawAxes();
            break;

        case GLFW_KEY_I:
            ctx->toggleDrawInfo();
            break;

        case GLFW_KEY_Q:
            ctx->stop();
            break;

        case GLFW_KEY_B:
            ctx->increasePointSize();
            break;

        case GLFW_KEY_S:
            ctx->decreasePointSize();
            break;

        case GLFW_KEY_L:
            ctx->toggleDrawMode();
            break;

        case GLFW_KEY_N:
            break;

        case GLFW_KEY_H:
      //case '?':
            //scene->drawHelp = !scene->drawHelp;
            break;

        case GLFW_KEY_O: /* Toggle camera following CM or on milkyway center */
            ctx->toggleOrigin();
            break;

        case GLFW_KEY_R: /* Toggle floating */
            ctx->toggleFloatMode();
            break;

        case GLFW_KEY_P:
            ctx->togglePaused();
            break;

        case GLFW_KEY_C:
            ctx->toggleMonochromatic();
            break;

        default:
            return;
    }
}

static void nbglSetHandlers(NBodyGraphics* graphicsContext)
{
    glfwSetErrorCallback(errorHandler);

    glfwSetWindowSizeCallback(resizeHandler);
    glfwSetWindowCloseCallback(closeHandler);

    glfwSetKeyCallback(keyHandler);
    //glfwSetKeyCallback(charHandler);
    glfwSetMouseButtonCallback(mouseButtonHandler);
    glfwSetMousePosCallback(mousePosHandler);
    glfwSetScrollCallback(scrollHandler);

    globalGraphicsContext = graphicsContext;
}

NBodyGraphics::NBodyGraphics(scene_t* scene, const VisArgs* args) : drawOptions(args)
{
    this->window = NULL;
    this->scene = scene;
    this->particleVAO = 0;
    this->whiteParticleVAO = 0;

    this->positionBuffer = 0;
    this->velocityBuffer = 0;
    this->accelerationBuffer = 0;
    this->colorBuffer = 0;

    this->particleTexture = 0;

    this->particleTextureProgram.program = 0;
    this->particleTextureProgram.positionLoc = -1;
    this->particleTextureProgram.colorLoc = -1;
    this->particleTextureProgram.modelToCameraMatrixLoc = -1;
    this->particleTextureProgram.cameraToClipMatrixLoc = -1;
    this->particleTextureProgram.particleTextureLoc = -1;
    this->particleTextureProgram.pointSizeLoc = -1;

    this->particlePointProgram.program = 0;
    this->particlePointProgram.positionLoc = -1;
    this->particlePointProgram.colorLoc = -1;
    this->particlePointProgram.modelToCameraMatrixLoc = -1;
    this->particlePointProgram.pointSizeLoc = -1;

    this->running = false;
    this->galaxyModel = NULL;
}

NBodyGraphics::~NBodyGraphics()
{
    GLuint buffers[5];

    glDeleteProgram(this->particleTextureProgram.program);
    glDeleteProgram(this->particlePointProgram.program);

    buffers[0] = this->positionBuffer;
    buffers[1] = this->velocityBuffer;
    buffers[2] = this->accelerationBuffer;
    buffers[3] = this->colorBuffer;
    buffers[4] = this->whiteBuffer;

    glDeleteBuffers(5, buffers);

    glDeleteVertexArrays(1, &this->particleVAO);
    glDeleteVertexArrays(1, &this->whiteParticleVAO);
    glDeleteTextures(1, &this->particleTexture);
}

void NBodyGraphics::loadShaders()
{
    this->particleTextureProgram.program = nbglCreateProgram("particle texture program",
                                                             (const char*) particle_vertex_glsl,
                                                             (const char*) particle_texture_fragment_glsl,
                                                             (GLint) particle_vertex_glsl_len,
                                                             (GLint) particle_texture_fragment_glsl_len);

    GLuint tprogram = this->particleTextureProgram.program;

    this->particleTextureProgram.positionLoc = glGetAttribLocation(tprogram, "position");
    this->particleTextureProgram.colorLoc = glGetAttribLocation(tprogram, "inputColor");
    this->particleTextureProgram.modelToCameraMatrixLoc = glGetUniformLocation(tprogram, "modelToCameraMatrix");
    this->particleTextureProgram.cameraToClipMatrixLoc = glGetUniformLocation(tprogram, "cameraToClipMatrix");
    this->particleTextureProgram.particleTextureLoc = glGetUniformLocation(tprogram, "particleTexture");
    this->particleTextureProgram.pointSizeLoc = glGetUniformLocation(tprogram, "pointSize");


    this->particlePointProgram.program = nbglCreateProgram("particle point program",
                                                           (const char*) particle_vertex_glsl,
                                                           (const char*) particle_point_fragment_glsl,
                                                           (GLint) particle_vertex_glsl_len,
                                                           (GLint) particle_point_fragment_glsl_len);

    GLuint pprogram = this->particlePointProgram.program;

    this->particlePointProgram.positionLoc = glGetAttribLocation(pprogram, "position");
    this->particlePointProgram.colorLoc = glGetAttribLocation(pprogram, "inputColor");
    this->particlePointProgram.modelToCameraMatrixLoc = glGetUniformLocation(pprogram, "modelToCameraMatrix");
    this->particlePointProgram.cameraToClipMatrixLoc = glGetUniformLocation(pprogram, "cameraToClipMatrix");
    this->particlePointProgram.pointSizeLoc = glGetUniformLocation(pprogram, "pointSize");
}

// create the VAO for the monochrome vs. not scene
void NBodyGraphics::prepareColoredVAO(GLuint& vao, GLuint color)
{
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glEnableVertexAttribArray(this->particleTextureProgram.positionLoc);
    glEnableVertexAttribArray(this->particleTextureProgram.colorLoc);

    glBindBuffer(GL_ARRAY_BUFFER, this->positionBuffer);
    /* 4th component is not included */
    glVertexAttribPointer(this->particleTextureProgram.positionLoc, 4, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, color);
    glVertexAttribPointer(this->particleTextureProgram.colorLoc, 4, GL_FLOAT, GL_FALSE, 0, 0);

    glBindVertexArray(0);
}

void NBodyGraphics::prepareVAOs()
{
    this->prepareColoredVAO(this->particleVAO, this->colorBuffer);
    this->prepareColoredVAO(this->whiteParticleVAO, this->whiteBuffer);
}

static int nbPopCircularQueue(NBodyCircularQueue* queue, int nbody, GLuint positionBuffer, SceneData* sceneData)
{
    int head = OPA_load_int(&queue->head);
    int tail = OPA_load_int(&queue->tail);

    if (head == tail)
    {
        return FALSE;  /* queue is empty */
    }

    const SceneInfo* info = &queue->info[head];
    const GLfloat* bodyData = (const GLfloat*) &queue->bodyData[head * nbody];

    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 4 * nbody * sizeof(GLfloat), bodyData);

    sceneData->currentTime = info->currentTime;
    sceneData->timeEvolve = info->timeEvolve;
    sceneData->centerOfMassView = glm::vec3(-info->rootCenterOfMass[0],
                                            -info->rootCenterOfMass[1],
                                            -info->rootCenterOfMass[2]);

    head = (head + 1) % NBODY_CIRC_QUEUE_SIZE;
    OPA_store_int(&queue->head, head);
    return TRUE;
}

bool NBodyGraphics::readSceneData()
{
    scene_t* scene = this->scene;
    return (bool) nbPopCircularQueue(&scene->queue, scene->nbody, this->positionBuffer,  &this->sceneData);
}

void NBodyGraphics::drawParticlesTextured(const glm::mat4& modelMatrix)
{
    assert(glGetError() == GL_NO_ERROR);
    glPointSize(this->drawOptions.texturedSpritePointSize);
    glUseProgram(this->particleTextureProgram.program);
    glUniformMatrix4fv(this->particleTextureProgram.modelToCameraMatrixLoc, 1, GL_FALSE, glm::value_ptr(modelMatrix));
    glUniformMatrix4fv(this->particleTextureProgram.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(cameraToClipMatrix));


    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, this->particleTexture);
    glUniform1i(this->particleTextureProgram.particleTextureLoc, 1);
    glUniform1f(this->particleTextureProgram.pointSizeLoc, this->drawOptions.texturedSpritePointSize);

    if (this->drawOptions.monochromatic)
    {
        glBindVertexArray(this->whiteParticleVAO);
    }
    else
    {
        glBindVertexArray(this->particleVAO);
    }

    glDepthMask(GL_FALSE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    //glBlendFunc(GL_ONE, GL_ONE);
    //float alpha = 0.5f;
    //glBlendFunc(GL_SRC_COLOR, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_DST_COLOR, GL_ONE_MINUS_SRC_ALPHA);

    //glBlendColor(1.0f - alpha, 1.0f - alpha, 1.0f - alpha, 1.0f);

    //glDepthMask(GL_FALSE);

    //glEnable(GL_POINT_SPRITE);
    //printf("enable point sprite %d\n", glGetError());
    //glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
    //printf("tex envi %d\n", glGetError());

    glDrawArrays(GL_POINTS, 0, this->scene->nbody);

    //glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    //glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    glBindVertexArray(0);
    glUseProgram(0);
}

void NBodyGraphics::drawParticlesPoints(const glm::mat4& modelMatrix)
{
    glPointSize(this->drawOptions.pointPointSize);

    glUseProgram(this->particlePointProgram.program);
    glUniformMatrix4fv(this->particlePointProgram.modelToCameraMatrixLoc, 1, GL_FALSE, glm::value_ptr(modelMatrix));
    glUniformMatrix4fv(this->particlePointProgram.cameraToClipMatrixLoc, 1, GL_FALSE, glm::value_ptr(cameraToClipMatrix));
    glUniform1f(this->particlePointProgram.pointSizeLoc, this->drawOptions.pointPointSize);

    if (this->drawOptions.monochromatic)
    {
        glBindVertexArray(this->whiteParticleVAO);
    }
    else
    {
        glBindVertexArray(this->particleVAO);
    }

    glDrawArrays(GL_POINTS, 0, this->scene->nbody);

    glBindVertexArray(0);
    glUseProgram(0);
}


void NBodyGraphics::drawParticles(const glm::mat4& modelMatrix)
{
    switch (this->drawOptions.drawMode)
    {
        case POINTS:
            this->drawParticlesPoints(modelMatrix);
            break;

        case TEXTURED_SPRITES:
            this->drawParticlesTextured(modelMatrix);
            break;

        default:
            mw_panic("Invalid draw mode\n");
    }
}

void NBodyGraphics::createBuffers()
{
    GLuint buffers[5];

    glGenBuffers(5, buffers);

    this->positionBuffer = buffers[0];
    this->velocityBuffer = buffers[1];
    this->accelerationBuffer = buffers[2];
    this->colorBuffer = buffers[3];
    this->whiteBuffer = buffers[4];
}

void NBodyGraphics::createPositionBuffer()
{
    GLint nbody = this->scene->nbody;

    glBindBuffer(GL_ARRAY_BUFFER, this->positionBuffer);
	glBufferData(GL_ARRAY_BUFFER, 4 * nbody * sizeof(GLfloat), NULL, GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void NBodyGraphics::populateBuffers()
{
    this->createPositionBuffer();
    this->particleTexture = createParticleTexture(32);
}

void NBodyGraphics::loadModel(GalaxyModel& model)
{
    model.generateModel();
    model.bufferData();
    model.loadShaders();
    model.prepareVAO();
    model.loadGalaxyTexture();
    this->galaxyModel = &model;
}

void NBodyGraphics::loadColors()
{
    GLint nbody = this->scene->nbody;

    /* assign random particle colors */
    srand((unsigned int) time(NULL));

    Color* color = new Color[nbody];

    for (GLint i = 0; i < nbody; ++i)
    {
        color[i].r = color[i].g = color[i].b = 1.0f;
        color[i].ignore = 1.0f;
    }

    // create a white buffer now
    // TODO: There is probably a better way to set a buffer to all 1
    glBindBuffer(GL_ARRAY_BUFFER, this->whiteBuffer);
    glBufferData(GL_ARRAY_BUFFER, 4 * nbody * sizeof(GLfloat), color, GL_STATIC_DRAW);

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

        if (false)
        {
            color[i].r = (GLfloat) R * scale;
            color[i].g = (GLfloat) G * scale;
            color[i].b = (GLfloat) B * scale;
        }
        else
        {
            color[i].r = (double) rand() / (double) RAND_MAX;
            color[i].g = (double) rand() / (double) RAND_MAX;
            color[i].b = (double) rand() / (double) RAND_MAX;
        }

        // TODO: Doesn't do anything yet?
        color[i].ignore = 1.0f;
    }

    glBindBuffer(GL_ARRAY_BUFFER, this->colorBuffer);
    glBufferData(GL_ARRAY_BUFFER, 4 * nbody * sizeof(GLfloat), color, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    delete[] color;
}

void NBodyGraphics::prepareContext()
{
    this->createBuffers();
    this->loadShaders();
    this->prepareVAOs();

    this->text.prepareConstantText(this->scene);

    this->axes.loadShader();
    this->axes.createBuffers();
    this->axes.prepareVAO();

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth(1.0);


    glEnable(GL_MULTISAMPLE);

    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE); // not sure what difference is between these
    glEnable(GL_PROGRAM_POINT_SIZE);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);

    float maxSmoothPointSize[2];
    glGetFloatv(GL_SMOOTH_POINT_SIZE_RANGE, (GLfloat*) &maxSmoothPointSize);
    printf("point size range %f %f\n", maxSmoothPointSize[0], maxSmoothPointSize[1]);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDepthRange(0.0f, 1.0f);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_CONSTANT_COLOR_EXT, GL_ONE_MINUS_SRC_COLOR);
    //glBlendFunc(GL_CONSTANT_COLOR_EXT, GL_ONE_MINUS_SRC_COLOR);

    //glEnable(GL_ALPHA_TEST);
    ///glEnable(GL_COLOR_MATERIAL);
}

static void requestGL32()
{
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3);
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 2);
    glfwOpenWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwOpenWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwOpenWindowHint(GLFW_DEPTH_BITS, 24);

    glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);

  #ifndef NDEBUG
    glfwOpenWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
  #endif
}

static GLFWwindow nbglPrepareWindow(const VisArgs* args)
{
    GLFWvidmode vidMode;
    glfwGetDesktopMode(&vidMode);
    int winMode = (args->fullscreen || args->plainFullscreen) ? GLFW_FULLSCREEN : GLFW_WINDOWED;

    int width = args->width == 0 ? 3 * vidMode.width / 4 : args->width;
    int height = args->height == 0 ? 3 * vidMode.height / 4 : args->height;

    requestGL32();
    return glfwOpenWindow(width, height, winMode, "Milkyway@Home N-body", NULL);
}

void NBodyGraphics::display()
{
    glm::mat4 modelMatrix = viewPole.CalcMatrix();

    if (this->drawOptions.cmCentered)
    {
        modelMatrix = glm::translate(modelMatrix, this->sceneData.centerOfMassView);
    }

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClearDepth(1.0);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    if (this->galaxyModel)
    {
        this->galaxyModel->draw(modelMatrix);
    }

    if (this->drawOptions.drawAxes)
    {
        this->axes.draw(modelMatrix);
    }

    if (this->drawOptions.drawParticles)
    {
        this->drawParticles(modelMatrix);
    }

    if (this->drawOptions.drawInfo)
    {
        // TODO: Check we have info / aren't a static scene
        this->text.drawProgressText(this->sceneData);
    }
}

void NBodyGraphics::mainLoop()
{
    this->running = true;

    while (this->running)
    {
        double t1 = glfwGetTime();

        glfwPollEvents();
        if (!this->drawOptions.paused)
        {
            this->readSceneData();
        }
        this->display();
        glfwSwapBuffers();

        double t2 = glfwGetTime();

        double dt = t2 - t1;
        dt = dt - (1.0 / 60.0);
        if ((int) dt > 0)
        {
            mwMilliSleep((int) dt);
        }
    }
}

int nbglGetExclusiveSceneAccess(scene_t* scene)
{
    int oldPID = OPA_cas_int(&scene->attachedPID, 0, (int) getpid());
    if (oldPID != 0)
    {
        mw_printf("Could not get exclusive access to simulation shared segment "
                  "(Owned by process %d)\n",
                  oldPID);

        /* TODO: We could check if this process is actually alive in
           case something went wrong and steal it if it is dead
         */
        return 1;
    }

    return 0;
}

#if !BOINC_APPLICATION

/* FIXME: Duplicated in nbody_shmem.c */
static size_t nbFindShmemSize(int nbody)
{
    size_t snapshotSize = sizeof(NBodyCircularQueue) + nbody * sizeof(FloatPos);
    return sizeof(scene_t) + NBODY_CIRC_QUEUE_SIZE * snapshotSize;
}

scene_t* nbglConnectSharedScene(int instanceId)
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

    if (sb.st_size < sizeof(scene_t) || sb.st_size < nbFindShmemSize(scene->nbody))
    {
        mw_printf("Shared memory segment is impossibly small ("ZU")\n", (size_t) sb.st_size);
        if (shm_unlink(name) < 0)
        {
            mwPerror("Unlink shared memory");
        }

        return NULL;
    }

    return scene;
}

#else

static scene_t* nbAttemptConnectSharedScene(void)
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
scene_t* nbglConnectSharedScene(int instanceId)
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

int nbglCheckConnectedVersion(const scene_t* scene)
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

static void nbglSetSceneSettings(scene_t* scene, const VisArgs* args)
{
    OPA_store_int(&scene->blockSimulationOnGraphics, args->blockSimulation);
}

int nbglRunGraphics(scene_t* scene, const VisArgs* args)
{
    if (!scene)
    {
        mw_printf("No scene to display\n");
        return 1;
    }

    if (!glfwInit())
    {
        mw_printf("Failed to initialize GLFW\n");
        return 1;
    }

    GLFWwindow window = nbglPrepareWindow(args);
    if (!window)
    {
        mw_printf("Failed to open window: %s\n", glfwErrorString(glfwGetError()));
        return 1;
    }

    nbglSetSceneSettings(scene, args);

    {
        // GL context needs to be open or else destructors will crash
        NBodyGraphics graphicsContext(scene, args);

        try
        {
            graphicsContext.prepareContext();
            nbglSetHandlers(&graphicsContext);
            graphicsContext.populateBuffers();
            graphicsContext.loadColors();

            GalaxyModel model;
            graphicsContext.loadModel(model);
            graphicsContext.mainLoop();
        }
        catch (const std::exception& e)
        {
            std::cerr << "Exception: " << e.what() << std::endl;
            glfwTerminate();
            return 1;
        }
    }

    globalGraphicsContext = NULL;
    glfwTerminate();

    return 0;
}

