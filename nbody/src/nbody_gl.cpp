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

#include "nbody_config.h"

#include "nbody_gl_includes.h"
#include "nbody_gl.h"
#include "nbody_graphics.h"
#include "nbody_gl_util.h"
#include "nbody_gl_axes.h"
#include "nbody_gl_text.h"
#include "nbody_gl_galaxy_model.h"
#include "nbody_gl_private.h"
#include "nbody_gl_includes.h"
#include "nbody_particle_texture.h"
#include "nbody_gl_shaders.h"
#include "milkyway_util.h"

#include <assert.h>
#include <iostream>
#include <stdexcept>

static const float zNear = 0.01f;
static const float zFar = 1000.0f;

static const glm::vec3 xAxis(1.0f, 0.0f, 0.0f);
static const glm::vec3 yAxis(0.0f, 1.0f, 0.0f);
static const glm::vec3 zAxis(0.0f, 0.0f, 1.0f);
static const glm::mat4 identityMatrix(1.0f);

static const glm::vec2 zeroVec2(0.0f, 0.0f);

glm::mat4 cameraToClipMatrix(1.0f);

static const float pointSizeChangeFactor = 1.05f;

// For access from callbacks
class NBodyGraphics;
static NBodyGraphics* globalGraphicsContext = NULL;


static const glm::vec3 origin(0.0f, 0.0f, 0.0f);
static const glm::fquat startOrientation = glm::fquat(1.0f, 0.0f, 0.0f, 0.0f);

static const float startRadius = 2.0f * 30.0f;

// angular component only
static const float startFloatSpeed = 0.2f;
static const float minFloatSpeed = 1.0e-3f;
static const float maxFloatSpeed = 10.0f;
static const float floatSpeedChangeFactor = 1.1f;


static const glutil::ViewData initialViewData =
{
    origin,  // center position
    startOrientation, // view direction
    startRadius,  // radius
    0.0f    // spin rotation of up axis
};

// limits of how far you can zoom to
static const float minViewRadius = 0.05f;
static const float maxViewRadius = 350.0f;


static const double minFloatTime = 3.0;
static const double maxFloatTime = 10.0;

static const float maxRadialFloatSpeed = 0.2f;

static glutil::ViewScale viewScale =
{
    minViewRadius, maxViewRadius,
    0.5f, 0.1f,    // radius deltas
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

    GalaxyModel* galaxyModel;

    NBodyText text;
    NBodyAxes axes;

    SceneData sceneData;

    // state of auto rotate store
    struct FloatState
    {
        glm::fquat orient;
        glm::vec2 angleDiff;
        double lastTime;
        double floatTime; // time this float should last
        float speed;
        float radius;
        float rSpeed;

        FloatState() : orient(startOrientation),
                       angleDiff(glm::vec2(0.0f, 0.0f)),
                       lastTime(-FLT_MAX),
                       floatTime(glm::linearRand(minFloatTime, maxFloatTime)),
                       speed(startFloatSpeed),
                       radius(0.5f * startRadius),
                       rSpeed(0.0f) { }
    } floatState;

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

    bool running;

    void loadShaders();
    void createBuffers();
    void prepareColoredVAO(GLuint& vao, GLuint color);
    void prepareVAOs();
    void createPositionBuffer();
    void drawParticlesTextured(const glm::mat4& modelMatrix);
    void drawParticlesPoints(const glm::mat4& modelMatrix);
    void loadColors();
    void calculateModelToCameraMatrix(glm::mat4& matrix);
    void newFloatDirection();

public:
    NBodyGraphics(scene_t* scene, const VisArgs* args);
    ~NBodyGraphics();

    void prepareContext();

    void drawAxes();
    void drawParticles(const glm::mat4& modelMatrix);
    bool readSceneData();
    void floatMotion();

    void display();
    void mainLoop();

    void stop()
    {
        this->running = false;
    }

    bool isScreensaver()
    {
        return this->drawOptions.screensaverMode;
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
            size *= pointSizeChangeFactor;
            if (size > 1000.0f)
            {
                size = 1000.0f;
            }

            this->drawOptions.texturedSpritePointSize = size;
        }
        else
        {
            float size = this->drawOptions.pointPointSize;
            size *= pointSizeChangeFactor;
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
            size /= pointSizeChangeFactor;
            if (size < 1.0e-3f)
            {
                size = 1.0e-3f;
            }

            this->drawOptions.texturedSpritePointSize = size;
        }
        else
        {
            float size = this->drawOptions.pointPointSize;

            size /= pointSizeChangeFactor;
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

    void toggleHelp()
    {
        this->drawOptions.drawHelp = !this->drawOptions.drawHelp;
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

    void increaseFloatSpeed()
    {
        this->floatState.speed *= floatSpeedChangeFactor;
        this->floatState.speed = glm::clamp(this->floatState.speed, minFloatSpeed, maxFloatSpeed);
    }

    void decreaseFloatSpeed()
    {
        this->floatState.speed /= floatSpeedChangeFactor;
        this->floatState.speed = glm::clamp(this->floatState.speed, minFloatSpeed, maxFloatSpeed);
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

    if (   glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS
        || glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS)
    {
        modifiers |= glutil::MM_KEY_SHIFT;
    }

    if (   glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS
        || glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS)
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
    if (globalGraphicsContext->isScreensaver())
    {
        globalGraphicsContext->stop();
    }

    int x, y;
    int modifiers = getGLFWModifiers(window);
    glfwGetMousePos(window, &x, &y);
    viewPole.MouseClick(glfwButtonToGLUtil(button), action == GLFW_PRESS, modifiers, glm::ivec2(x, y));
}

static void mousePosHandler(GLFWwindow window, int x, int y)
{
    if (globalGraphicsContext->isScreensaver())
    {
        int w, h;

        /* In fullscreen mode the cursor seems to start at the corner
         * and get constantly reset to the midpoint */
        glfwGetWindowSize(window, &w, &h);

        if (   !(x == w / 2 && y == h / 2)  /* Not midpoint */
            && !(x == 0 && y == 0))         /* Not corner */
        {
            globalGraphicsContext->stop();
            return;
        }
    }

    viewPole.MouseMove(glm::ivec2(x, y));
}

static void scrollHandler(GLFWwindow window, int x, int y)
{
    int direction = y > 0;
    viewPole.MouseWheel(direction, getGLFWModifiers(window), glm::ivec2(x, y));
}

static void keyHandler(GLFWwindow window, int key, int pressed)
{
    if (!pressed) // release
        return;

    NBodyGraphics* ctx = globalGraphicsContext;

    if (ctx->isScreensaver())
    {
        ctx->stop();
    }

    switch (key)
    {
        case GLFW_KEY_A:
            ctx->toggleDrawAxes();
            break;

        case GLFW_KEY_I:
            ctx->toggleDrawInfo();
            break;

        case GLFW_KEY_ESC:
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
        case GLFW_KEY_SLASH: /* Same as question mark key */
            ctx->toggleHelp();
            break;

        case GLFW_KEY_O: /* Toggle camera following CM or on milkyway center */
            ctx->toggleOrigin();
            break;

        case GLFW_KEY_F: /* Toggle floating */
            ctx->toggleFloatMode();
            break;

        case GLFW_KEY_Z:
            ctx->increaseFloatSpeed();
            break;

        case GLFW_KEY_X:
            ctx->decreaseFloatSpeed();
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

static void nbglSetHandlers()
{
    glfwSetWindowSizeCallback(resizeHandler);
    glfwSetWindowCloseCallback(closeHandler);

    glfwSetKeyCallback(keyHandler);
    //glfwSetKeyCallback(charHandler);
    glfwSetMouseButtonCallback(mouseButtonHandler);
    glfwSetMousePosCallback(mousePosHandler);
    glfwSetScrollCallback(scrollHandler);
}

NBodyGraphics::NBodyGraphics(scene_t* scene, const VisArgs* args)
    : sceneData(SceneData((bool) scene->staticScene)),
      drawOptions(args)
{
    this->scene = scene;
    this->loadShaders();
    this->createBuffers();
    this->prepareVAOs();

    this->createPositionBuffer();
    this->loadColors();
    this->particleTexture = createParticleTexture(32);

    this->galaxyModel = scene->hasGalaxy ? new GalaxyModel() : NULL;

    this->running = true;
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

    delete this->galaxyModel;
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
    this->text.prepareConstantText(this->scene);

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
}

void NBodyGraphics::calculateModelToCameraMatrix(glm::mat4& matrix)
{
    // we skip using viewPole.CalcMatrix() because
    // we want to combine the orientation directly from the view pole
    // with the orientation we randomly generate while floating


    // We should only get really close or really far away if manually zooming.
    // The float radius should be kept < the current view radius
    // It should stay a reasonable viewing distance if left on its own
    float effRadius = 0.5f * (viewPole.GetView().radius + this->floatState.radius);

    //In this space, we are facing in the correct direction. Which means that the camera point
    //is directly behind us by the radius number of units.
    matrix = glm::translate(identityMatrix, glm::vec3(0.0f, 0.0f, -effRadius));

    //Rotate the world to look in the right direction..
    glm::fquat totalOrientation = this->floatState.orient * viewPole.GetView().orient;
    glm::fquat fullRotation = glm::angleAxis(-viewPole.GetView().degSpinRotation, zAxis) * totalOrientation;

    matrix *= glm::mat4_cast(fullRotation);

    //Translate the world by the negation of the lookat point, placing the origin at the
    //lookat point.

    if (this->drawOptions.cmCentered)
    {
        matrix = glm::translate(matrix, this->sceneData.centerOfMassView);
    }
}

void NBodyGraphics::display()
{
    glm::mat4 modelMatrix;

    this->calculateModelToCameraMatrix(modelMatrix);

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
        this->text.drawProgressText(this->sceneData);
    }
}

void NBodyGraphics::mainLoop()
{
    while (this->running)
    {
        double t1 = glfwGetTime();

        glfwPollEvents();
        if (!this->drawOptions.paused)
        {
            this->readSceneData();
        }

        if (this->drawOptions.floatMode)
        {
            this->floatMotion();
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

void NBodyGraphics::newFloatDirection()
{
    glm::vec2& angleDiff = this->floatState.angleDiff;

    angleDiff.x = glm::linearRand(0.0f, this->floatState.speed);
    angleDiff.y = glm::linearRand(0.0f, this->floatState.speed);
    angleDiff = this->floatState.speed * glm::normalize(angleDiff);

    this->floatState.floatTime = glm::linearRand(minFloatTime, maxFloatTime);

    // make sure we are trying to move out if close
    if (this->floatState.radius <= 2.0f * minViewRadius)
    {
        this->floatState.rSpeed = glm::gaussRand(0.0f, maxRadialFloatSpeed);
    }
    else
    {
        // spend some periods at the same radius
        if (rand() > RAND_MAX / 2)
        {
            this->floatState.rSpeed = glm::gaussRand(-maxRadialFloatSpeed, maxRadialFloatSpeed);
        }
        else
        {
            this->floatState.rSpeed = 0.0f;
        }
    }
}

void NBodyGraphics::floatMotion()
{
    const glm::vec2& diff = this->floatState.angleDiff;

    this->floatState.orient = this->floatState.orient * glm::angleAxis(diff.x, yAxis);
    this->floatState.orient = glm::angleAxis(diff.y, xAxis) * this->floatState.orient;

    this->floatState.radius += this->floatState.rSpeed;
    this->floatState.radius = glm::clamp(this->floatState.radius,
                                         minViewRadius,
                                         viewPole.GetView().radius);

    double t = glfwGetTime();

    if (t - this->floatState.lastTime < this->floatState.floatTime)
    {
        return;
    }

    this->floatState.lastTime = t;

    this->newFloatDirection();
}

static void nbglSetSceneSettings(scene_t* scene, const VisArgs* args)
{
    OPA_store_int(&scene->blockSimulationOnGraphics, args->blockSimulation);
}

static void nbglRequestGL32()
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
    const char* title = "Milkyway@Home N-body";
    GLFWvidmode vidMode;

    glfwGetDesktopMode(&vidMode);
    nbglRequestGL32();

    if (args->fullscreen || args->plainFullscreen)
    {
        int width = args->width == 0 ? vidMode.width : args->width;
        int height = args->height == 0 ? vidMode.height : args->height;

        return glfwOpenWindow(width, height, GLFW_FULLSCREEN, title, NULL);
    }
    else
    {
        int width = args->width == 0 ? 3 * vidMode.width / 4 : args->width;
        int height = args->height == 0 ? 3 * vidMode.height / 4 : args->height;

        return glfwOpenWindow(width, height, GLFW_WINDOWED, title, NULL);
    }
}

int nbglRunGraphics(scene_t* scene, const VisArgs* args)
{
    if (!scene)
    {
        mw_printf("No scene to display\n");
        return 1;
    }

    glfwSetErrorCallback(errorHandler);

    if (!glfwInit())
    {
        mw_printf("Failed to initialize GLFW\n");
        return 1;
    }

    GLFWwindow window = nbglPrepareWindow(args);
    if (!window)
    {
        return 1;
    }

    nbglSetSceneSettings(scene, args);

    try
    {
        // GL context needs to be open or else destructors will crash
        NBodyGraphics graphicsContext(scene, args);

        graphicsContext.prepareContext();
        graphicsContext.readSceneData();

        globalGraphicsContext = &graphicsContext;
        nbglSetHandlers();
        graphicsContext.mainLoop();
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        glfwTerminate();
        return 1;
    }

    globalGraphicsContext = NULL;
    glfwTerminate();

    return 0;
}

