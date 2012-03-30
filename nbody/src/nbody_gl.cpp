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
#include "nbody_gl_orbit_trace.h"
#include "nbody_gl_text.h"
#include "nbody_gl_galaxy_model.h"
#include "nbody_gl_private.h"
#include "nbody_gl_includes.h"
#include "nbody_gl_particle_texture.h"
#include "nbody_gl_shaders.h"
#include "milkyway_util.h"

#include <assert.h>
#include <iostream>
#include <stdexcept>

static const float zNearView = 0.01f;
static const float zFarView = 1000.0f;

static const glm::vec3 xAxis(1.0f, 0.0f, 0.0f);
static const glm::vec3 yAxis(0.0f, 1.0f, 0.0f);
static const glm::vec3 zAxis(0.0f, 0.0f, 1.0f);
static const glm::mat4 identityMatrix(1.0f);

static const glm::vec2 zeroVec2(0.0f, 0.0f);

glm::mat4 cameraToClipMatrix(1.0f);

static const float minPointPointSize = 1.0e-3f;
static const float maxPointPointSize = 500.0f;

static const float minTexturedPointSize = 1.0e-3f;
static const float maxTexturedPointSize = 2000.0f;

static const float pointSizeChangeFactor = 1.05f;

static const glm::vec3 origin(0.0f, 0.0f, 0.0f);
static const glm::fquat identityOrientation = glm::fquat(1.0f, 0.0f, 0.0f, 0.0f);

static const float startRadius = 30.0f;

// angular component only, degrees per second
static const float minFloatSpeed = 1.0e-3f;
static const float maxFloatSpeed = 30.0f;
static const float floatSpeedChangeFactor = 1.1f;


// limits of how far you can zoom to
static const float minViewRadius = 0.05f;
static const float maxViewRadius = 350.0f;

// radius where floating is considered too close and should move away
static const float closeViewRadiusThreshold = 5.0f;


static const double minFloatTime = 2.0;
static const double maxFloatTime = 10.0;

static const float maxRadialFloatSpeed = 0.2f;

static const glutil::ViewScale viewScale =
{
    minViewRadius, maxViewRadius,
    2.5f, 0.5f,    // radius deltas
    4.0f, 1.0f,
    90.0f / 250.0f // rotation scale
};

static const glutil::ViewData initialViewData =
{
    origin,              // center position
    identityOrientation, // view direction
    startRadius,         // radius
    0.0f                 // spin rotation of up axis
};

class NBodyGraphics;
static NBodyGraphics* globalGraphicsContext;


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
    OrbitTrace orbitTrace;

    SceneData sceneData;
    glutil::ViewPole viewPole;

    // state of auto rotate store
    struct FloatState
    {
        glm::vec2 angleVec;
        double lastShiftTime;    // time since direction changed
        double lastUpdateTime;   // time since angle updated
        double floatTime; // time this float should last
        float speed;
        float rSpeed;

        FloatState(const VisArgs* args)
        : angleVec(glm::vec2(0.0f, 0.0f)),
          lastShiftTime(0.0f),
          lastUpdateTime(0.0f),
          floatTime(glm::linearRand(minFloatTime, maxFloatTime)),
          speed(args->floatSpeed),
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
        bool floatMode;

        bool cmCentered;
        bool monochromatic;

        DrawMode drawMode;
        // Keep a separate point size for each draw mode
        float texturedSpritePointSize;
        float pointPointSize;

        // for resetting to the same size as at the start
        float texturedSpritePointSizeOrig;
        float pointPointSizeOrig;

        bool drawInfo;
        bool drawAxes;
        bool drawOrbitTrace;
        bool drawParticles;
        bool drawHelp;
        bool quitOnSimulationComplete;

        DrawOptions(const VisArgs* args)
        : screensaverMode(args->fullscreen && !args->plainFullscreen),
          floatMode((bool) !args->noFloat),
          cmCentered((bool) !args->originCentered),
          monochromatic((bool) args->monochromatic),
          drawMode(args->untexturedPoints ? POINTS : TEXTURED_SPRITES),
          texturedSpritePointSize(args->texturedPointSize),
          pointPointSize(args->pointPointSize),
          texturedSpritePointSizeOrig(args->texturedPointSize),
          pointPointSizeOrig(args->pointPointSize),
          drawInfo((bool) !args->noDrawInfo),
          drawAxes((bool) args->drawAxes),
          drawOrbitTrace((bool) args->drawOrbitTrace),
          drawParticles(true),
          drawHelp(false),
          quitOnSimulationComplete((bool) args->quitOnComplete)
        { }
    } drawOptions;

    bool running;
    bool needsUpdate;
    bool paused;
    int eventPollPeriod;

    void loadShaders();
    void createBuffers();
    void prepareColoredVAO(GLuint& vao, GLuint color);
    void prepareVAOs();
    void createPositionBuffer();
    void drawParticlesTextured(const glm::mat4& modelMatrix);
    void drawParticlesPoints(const glm::mat4& modelMatrix);
    void loadColors();
    void newFloatDirection();
    void resetFloatState();
    void findInitialOrientation();
    void floatMotion();
    void prepareContext();
    bool readSceneData();

public:
    NBodyGraphics(scene_t* scene, const VisArgs* args);
    ~NBodyGraphics();

    inline void mouseClick(glutil::MouseButtons button, bool pressed, int modifiers, int x, int y)
    {
        this->viewPole.MouseClick(button, pressed, modifiers, glm::ivec2(x, y));
    }

    inline void mouseMove(int x, int y)
    {
        this->viewPole.MouseMove(glm::ivec2(x, y));
    }

    inline void mouseWheel(int modifiers, int x, int y)
    {
        int direction = y > 0;
        this->viewPole.MouseWheel(direction, modifiers, glm::ivec2(x, y));
    }

    void drawAxes();
    void drawParticles(const glm::mat4& modelMatrix);

    void display();
    void mainLoop();

    void stop()
    {
        this->running = false;
    }

    inline void markDirty()
    {
        this->needsUpdate = true;
    }

    inline void markClean()
    {
        this->needsUpdate = false;
    }

    bool isScreensaver()
    {
        return this->drawOptions.screensaverMode;
    }

    void toggleDrawMode()
    {
        this->markDirty();

        if (this->drawOptions.drawMode == POINTS)
        {
            this->drawOptions.drawMode = TEXTURED_SPRITES;
        }
        else
        {
            this->drawOptions.drawMode = POINTS;
        }
    }

    void resetPointSize()
    {
        this->markDirty();

        if (this->drawOptions.drawMode == TEXTURED_SPRITES)
        {
            this->drawOptions.texturedSpritePointSize = this->drawOptions.texturedSpritePointSizeOrig;
        }
        else
        {
            this->drawOptions.pointPointSize = this->drawOptions.pointPointSizeOrig;
        }
    }

    void increasePointSize()
    {
        this->markDirty();

        if (this->drawOptions.drawMode == TEXTURED_SPRITES)
        {
            float size = this->drawOptions.texturedSpritePointSize;
            size = glm::clamp(size * pointSizeChangeFactor, minTexturedPointSize, maxTexturedPointSize);
            this->drawOptions.texturedSpritePointSize = size;
        }
        else
        {
            float size = this->drawOptions.pointPointSize;
            size = glm::clamp(size * pointSizeChangeFactor, minPointPointSize, maxPointPointSize);
            this->drawOptions.pointPointSize = size;
        }
    }

    void decreasePointSize()
    {
        this->markDirty();

        if (this->drawOptions.drawMode == TEXTURED_SPRITES)
        {
            float size = this->drawOptions.texturedSpritePointSize;
            size = glm::clamp(size / pointSizeChangeFactor, minTexturedPointSize, maxTexturedPointSize);
            this->drawOptions.texturedSpritePointSize = size;
        }
        else
        {
            float size = this->drawOptions.pointPointSize;
            size = glm::clamp(size / pointSizeChangeFactor, minPointPointSize, maxPointPointSize);
            this->drawOptions.pointPointSize = size;
        }
    }

    void toggleDrawAxes()
    {
        this->markDirty();
        this->drawOptions.drawAxes = !this->drawOptions.drawAxes;
    }

    void toggleDrawOrbitTrace()
    {
        this->markDirty();
        this->drawOptions.drawOrbitTrace = !this->drawOptions.drawOrbitTrace;
    }

    void toggleDrawInfo()
    {
        this->markDirty();
        this->drawOptions.drawInfo = !this->drawOptions.drawInfo;
    }

    void toggleHelp()
    {
        this->markDirty();
        this->drawOptions.drawHelp = !this->drawOptions.drawHelp;
    }

    void toggleMonochromatic()
    {
        this->markDirty();
        this->drawOptions.monochromatic = !this->drawOptions.monochromatic;
    }

    void togglePaused()
    {
        this->markDirty();
        this->paused = !this->paused;
        OPA_store_int(&this->scene->paused, (int) this->paused);

        if (this->drawOptions.floatMode)
        {
            this->resetFloatState();
        }
    }

    void toggleFloatMode()
    {
        this->drawOptions.floatMode = !this->drawOptions.floatMode;

        if (this->drawOptions.floatMode)
        {
            this->resetFloatState();
        }
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
        this->markDirty();
        this->drawOptions.cmCentered = !this->drawOptions.cmCentered;

        if (this->drawOptions.cmCentered)
        {
            this->viewPole.SetOrigin(this->sceneData.centerOfMass);
        }
        else
        {
            this->viewPole.SetOrigin(origin);
        }
    }
};

static void errorHandler(int errorCode, const char* msg)
{
    fprintf(stderr, "GLFW error (%d): %s\n", errorCode, msg);
}

static void resizeHandler(GLFWwindow window, int w, int h)
{
    (void) window;

    float wf = (float) w;
    float hf = (float) h;
    float aspectRatio = wf / hf;
    cameraToClipMatrix = glm::perspective(90.0f, aspectRatio, zNearView, zFarView);

    const float fontHeight = 12.0f; // FIXME: hardcoded font height
    textCameraToClipMatrix = glm::ortho(0.0f, wf, -hf, 0.0f);
    textCameraToClipMatrix = glm::translate(textCameraToClipMatrix, glm::vec3(0.0f, -fontHeight, 0.0f));

    glViewport(0, 0, (GLsizei) w, (GLsizei) h);


    globalGraphicsContext->display();
    glfwSwapBuffers();
}

static int closeHandler(GLFWwindow window)
{
    (void) window;

    globalGraphicsContext->stop();
    return 0;
}

static int getGLFWModifiers(GLFWwindow window)
{
    (void) window;

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
    NBodyGraphics* ctx = globalGraphicsContext;

    if (ctx->isScreensaver())
    {
        ctx->stop();
    }

    int x, y;
    int modifiers = getGLFWModifiers(window);
    glfwGetMousePos(window, &x, &y);
    ctx->mouseClick(glfwButtonToGLUtil(button), action == GLFW_PRESS, modifiers, x, y);
}

static void mousePosHandler(GLFWwindow window, int x, int y)
{
    (void) window;

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

    globalGraphicsContext->mouseMove(x, y);
    globalGraphicsContext->markDirty();
}

static void scrollHandler(GLFWwindow window, int x, int y)
{
    globalGraphicsContext->mouseWheel(getGLFWModifiers(window), x, y);
    globalGraphicsContext->markDirty();
}

static void keyHandler(GLFWwindow window, int key, int pressed)
{
    (void) window;

    if (!pressed) // release
        return;

    NBodyGraphics* ctx = globalGraphicsContext;

    switch (key)
    {
        case GLFW_KEY_ESC:
            ctx->stop();
            break;

        case GLFW_KEY_PAUSE:
            ctx->togglePaused();
            break;

        case GLFW_KEY_F1:
            ctx->toggleHelp();
            break;

            // TODO: We should fix the mouse poles to use these
            // special keys for controls
        case GLFW_KEY_PAGE_UP:
        case GLFW_KEY_PAGE_DOWN:
        case GLFW_KEY_UP:
        case GLFW_KEY_DOWN:
        case GLFW_KEY_LEFT:
        case GLFW_KEY_RIGHT:
            break;

       /* Any non-character keys that aren't used should be here. If a
          key doesn't have a simple function and we are a screensaver,
          we quit
        */
        case GLFW_KEY_ENTER:
        case GLFW_KEY_BACKSLASH:
        case GLFW_KEY_RIGHT_BRACKET:
        case GLFW_KEY_GRAVE_ACCENT:
        case GLFW_KEY_WORLD_1:
        case GLFW_KEY_WORLD_2:
        case GLFW_KEY_TAB:
        case GLFW_KEY_BACKSPACE:
        case GLFW_KEY_INSERT:
        case GLFW_KEY_DELETE:
        case GLFW_KEY_HOME:
        case GLFW_KEY_END:
        case GLFW_KEY_CAPS_LOCK:
        case GLFW_KEY_SCROLL_LOCK:
        case GLFW_KEY_NUM_LOCK:
        case GLFW_KEY_PRINT_SCREEN:

        case GLFW_KEY_F2:
        case GLFW_KEY_F3:
        case GLFW_KEY_F4:
        case GLFW_KEY_F5:
        case GLFW_KEY_F6:
        case GLFW_KEY_F7:
        case GLFW_KEY_F8:
        case GLFW_KEY_F9:
        case GLFW_KEY_F10:
        case GLFW_KEY_F11:
        case GLFW_KEY_F12:
        case GLFW_KEY_F13:
        case GLFW_KEY_F14:
        case GLFW_KEY_F15:
        case GLFW_KEY_F16:
        case GLFW_KEY_F17:
        case GLFW_KEY_F18:
        case GLFW_KEY_F19:
        case GLFW_KEY_F20:
        case GLFW_KEY_F21:
        case GLFW_KEY_F22:
        case GLFW_KEY_F23:
        case GLFW_KEY_F24:
        case GLFW_KEY_F25:

        case GLFW_KEY_KP_0:
        case GLFW_KEY_KP_1:
        case GLFW_KEY_KP_2:
        case GLFW_KEY_KP_3:
        case GLFW_KEY_KP_4:
        case GLFW_KEY_KP_5:
        case GLFW_KEY_KP_6:
        case GLFW_KEY_KP_7:
        case GLFW_KEY_KP_8:
        case GLFW_KEY_KP_9:
        case GLFW_KEY_KP_DECIMAL:
        case GLFW_KEY_KP_DIVIDE:
        case GLFW_KEY_KP_MULTIPLY:
        case GLFW_KEY_KP_SUBTRACT:
        case GLFW_KEY_KP_ADD:
        case GLFW_KEY_KP_ENTER:

        case GLFW_KEY_MENU:
            if (ctx->isScreensaver())
            {
                ctx->stop();
            }
            break;

        default:
            return;
    }
}

static void charHandler(GLFWwindow window, int charCode)
{
    (void) window;

    NBodyGraphics* ctx = globalGraphicsContext;

    // the mousepoles have some controls for manually moving the focus
    // point but we aren't using that

    switch (charCode)
    {
        case 'A':
        case 'a':
            ctx->toggleDrawAxes();
            break;

        case 'T':
        case 't':
            ctx->toggleDrawOrbitTrace();
            break;

        case 'I':
        case 'i':
            ctx->toggleDrawInfo();
            break;

        case 'Q':
        case 'q':
            ctx->stop();
            break;

        case 'B':
        case 'b':
        case '+':
            ctx->increasePointSize();
            break;

        case 'S':
        case 's':
        case '-':
            ctx->decreasePointSize();
            break;

        case '=':
            ctx->resetPointSize();
            break;

        case 'D':
        case 'd':
            ctx->toggleDrawMode();
            break;

        case 'H':
        case 'h':
        case '?':
            ctx->toggleHelp();
            break;

        case 'O':
        case 'o': /* Toggle camera following CM or on milkyway center */
            ctx->toggleOrigin();
            break;

        case 'F':
        case 'f': /* Toggle floating */
            ctx->toggleFloatMode();
            break;

        case 'X':
        case 'x':
        case '>':
            ctx->increaseFloatSpeed();
            break;

        case 'Z':
        case 'z':
        case '<':
            ctx->decreaseFloatSpeed();
            break;

        case 'P':
        case 'p':
        case ' ':
            ctx->togglePaused();
            break;

        case 'C':
        case 'c':
            ctx->toggleMonochromatic();
            break;

        default:
            if (ctx->isScreensaver())
            {
                ctx->stop();
            }

            return;
    }
}

static void nbglSetHandlers()
{
    glfwSetWindowSizeCallback(resizeHandler);
    glfwSetWindowCloseCallback(closeHandler);

    glfwSetKeyCallback(keyHandler);
    glfwSetCharCallback(charHandler);
    glfwSetMouseButtonCallback(mouseButtonHandler);
    glfwSetMousePosCallback(mousePosHandler);
    glfwSetScrollCallback(scrollHandler);
}

NBodyGraphics::NBodyGraphics(scene_t* scene_, const VisArgs* args)
    : scene(scene_),
      text(NBodyText(&robotoRegular12)),
      orbitTrace(OrbitTrace(scene)),
      sceneData(SceneData((bool) scene->staticScene)),
      viewPole(glutil::ViewPole(initialViewData, viewScale, glutil::MB_LEFT_BTN)),
      floatState(FloatState(args)),
      drawOptions(args),
      running(true),
      needsUpdate(true),
      paused(false),
      eventPollPeriod(glm::clamp(args->eventPollPeriod, 0, MAX_EVENT_POLL_PERIOD))
{
    this->loadShaders();
    this->createBuffers();
    this->prepareVAOs();

    this->createPositionBuffer();
    this->loadColors();
    this->particleTexture = nbglCreateParticleTexture(32);

    this->galaxyModel = scene->hasGalaxy ? new GalaxyModel() : NULL;

    this->prepareContext();

    // read initial scene data before starting
    this->readSceneData();
    this->findInitialOrientation();
    this->resetFloatState();
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

// return TRUE if something was popped from the queue, FALSE if it was empty
static int nbPopCircularQueue(NBodyCircularQueue* queue,
                              int nbody,
                              GLuint positionBuffer,
                              SceneData* sceneData,
                              OrbitTrace* trace)
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
    sceneData->centerOfMass = glm::vec3(info->rootCenterOfMass[0],
                                        info->rootCenterOfMass[1],
                                        info->rootCenterOfMass[2]);

    trace->addPoint(info->rootCenterOfMass);

    head = (head + 1) % NBODY_CIRC_QUEUE_SIZE;
    OPA_store_int(&queue->head, head);
    return TRUE;
}

bool NBodyGraphics::readSceneData()
{
    bool success;
    success = (bool) nbPopCircularQueue(&this->scene->queue, this->scene->nbody,
                                        this->positionBuffer, &this->sceneData, &this->orbitTrace);

    if (success)
    {
        if (this->drawOptions.cmCentered)
        {
            this->viewPole.SetOrigin(this->sceneData.centerOfMass);
        }

        this->markDirty();
    }

    return success;
}

void NBodyGraphics::drawParticlesTextured(const glm::mat4& modelMatrix)
{
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

    glDrawArrays(GL_POINTS, 0, this->scene->nbody);

    glDepthMask(GL_TRUE);
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


static inline glm::vec3 nbglBetweenHSVColors(const glm::vec3& a, const glm::vec3& b)
{
    return glm::rgbColor(glm::vec3(glm::linearRand(a.x, b.x),
                                   glm::linearRand(a.y, b.y),
                                   glm::linearRand(a.z, b.z)));
}

static void nbglRandomParticleColor(Color& color, bool ignore)
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

    // HSV colors
    static const glm::vec3 light(204.0f, 0.13f, 0.84f);
    static const glm::vec3 dark(209.0f, 0.41f, 0.50f);

    static const glm::vec3 coreLight(40.0f, 0.01f, 0.99f);
    static const glm::vec3 coreDark(33.0f, 0.29f, 0.58f);

    glm::vec3 rgb;
    if (glm::linearRand(0.0f, 1.0f) > 0.2f)
    {
        // bluish
        rgb = nbglBetweenHSVColors(dark, light);
    }
    else
    {
        // reddish
        rgb = nbglBetweenHSVColors(coreDark, coreLight);
    }

    color.r = rgb.x;
    color.g = rgb.y;
    color.b = rgb.z;
    color.ignore = 1.0f;
    //color.ignore = ignore ? ?????
}

void NBodyGraphics::loadColors()
{
    GLint nbody = this->scene->nbody;
    const bool approxRealStarColors = true;

    /* assign random particle colors */
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

        if (approxRealStarColors)
        {
            // TODO: ignore doesn't do anything yet
            nbglRandomParticleColor(color[i], false);
        }
        else
        {
            if (false)
            {
                color[i].r = (GLfloat) (R * scale);
                color[i].g = (GLfloat) (G * scale);
                color[i].b = (GLfloat) (B * scale);
            }
            else
            {
                color[i].r = (float) rand() / (float) RAND_MAX;
                color[i].g = (float) rand() / (float) RAND_MAX;
                color[i].b = (float) rand() / (float) RAND_MAX;
            }
        }
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

    glEnable(GL_LINE_SMOOTH);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);

    //float maxSmoothPointSize[2];
    //glGetFloatv(GL_SMOOTH_POINT_SIZE_RANGE, (GLfloat*) &maxSmoothPointSize);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDepthRange(0.0f, 1.0f);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void NBodyGraphics::display()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClearDepth(1.0);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    glm::mat4 modelMatrix = this->viewPole.CalcMatrix();

    if (this->galaxyModel)
    {
        this->galaxyModel->draw(modelMatrix);
    }

    if (this->drawOptions.drawAxes)
    {
        this->axes.draw(modelMatrix);
    }

    if (this->drawOptions.drawOrbitTrace)
    {
        this->orbitTrace.drawTrace(modelMatrix);
    }

    if (this->drawOptions.drawParticles)
    {
        this->drawParticles(modelMatrix);
    }

    // FIXME? Maybe? Doesn't show info text when displaying help
    if (this->drawOptions.drawHelp)
    {
        this->text.drawHelpText();
    }
    else
    {
        if (this->drawOptions.drawInfo)
        {
            this->text.drawProgressText(this->sceneData);
        }
    }

    this->markClean();
}

// We need to throw out old float state whenever an event happens
// which could result in the float state not being updated for a
// significant period of time or else there will be a visible jump
void NBodyGraphics::resetFloatState()
{
    double now = glfwGetTime();
    // update the update times avoid a small jump on the first frame
    this->floatState.lastUpdateTime = now;
    this->floatState.lastShiftTime = now;
    this->newFloatDirection();
}

void NBodyGraphics::mainLoop()
{
    while (true)
    {
        while (!this->needsUpdate)
        {
            glfwPollEvents();

            if (!this->running)
            {
                /* Make sure to test right after polling to avoid
                 * waiting for the scene to update before quitting.
                 * This can leave us stuck with quitting not working
                 * after the simulation ends
                 */
                return;
            }

            if (!this->paused)
            {
                this->readSceneData();

                if (this->drawOptions.floatMode)
                {
                    this->floatMotion();
                }

                /* We will not quit if on completion if paused */
                if (this->drawOptions.quitOnSimulationComplete)
                {
                    if (OPA_load_int(&this->scene->ownerPID) == 0)
                    {
                        return;
                    }
                }
            }

            mwMilliSleep(this->eventPollPeriod);
        }

        this->display();
        glfwSwapBuffers();
    }
}

void NBodyGraphics::newFloatDirection()
{
    glm::vec2& angleVec = this->floatState.angleVec;

    angleVec.x = glm::linearRand(0.0f, 1.0f);
    angleVec.y = glm::linearRand(-1.0f, 0.0f);
    angleVec = glm::normalize(angleVec);

    this->floatState.floatTime = glm::linearRand(minFloatTime, maxFloatTime);

    // make sure we are trying to move out if close
    if (this->viewPole.GetView().radius <= closeViewRadiusThreshold)
    {
        this->floatState.rSpeed = glm::linearRand(0.3f * maxRadialFloatSpeed, maxRadialFloatSpeed);
    }
    else
    {
        // spend some periods at the same radius
        if (rand() > RAND_MAX / 2)
        {
            this->floatState.rSpeed = glm::linearRand(-maxRadialFloatSpeed, maxRadialFloatSpeed);
        }
        else
        {
            this->floatState.rSpeed = 0.0f;
        }
    }
}

void NBodyGraphics::floatMotion()
{
    double now = glfwGetTime();
    double dt = now - this->floatState.lastUpdateTime;
    this->floatState.lastUpdateTime = now;

    glm::vec2 angle = this->floatState.speed * this->floatState.angleVec * (float) dt;
    glm::fquat localOrientation = glm::angleAxis(angle.x, yAxis);
    localOrientation = glm::angleAxis(angle.y, xAxis) * localOrientation;

    this->viewPole.ApplyExternalOrientation(localOrientation);

    float dr = this->floatState.rSpeed * (float) dt;

    this->viewPole.ApplyExternalRadiusDelta(dr);


    // time since direction changed
    if (now - this->floatState.lastShiftTime >= this->floatState.floatTime)
    {
        this->floatState.lastShiftTime = now;
        this->newFloatDirection();
    }

    this->markDirty();
}

// we'll use the float orientation to start out oriented facing the
// origin from somewhere behind the center of mass
void NBodyGraphics::findInitialOrientation()
{
    glm::vec3 centerOfMass = this->sceneData.centerOfMass;
    glm::vec3 centerOfMassDir = glm::normalize(centerOfMass);

    // find angle between the center of mass and the plane of the milkyway
    float angle = atanf(fabsf(centerOfMass.z / centerOfMass.x));
    angle = glm::degrees(angle);

    // find an angle somewhat close to that angle
    float rangeAngle = glm::gaussRand(angle, 4.0f);

    // we want to eliminate the z component or else the starting
    // orientation makes rotation feel funny at the start
    glm::vec3 cmCrossZ = glm::normalize(glm::cross(centerOfMass, zAxis));

    float x = glm::linearRand(-1.0f, 1.0f);
    float y = glm::linearRand(-1.0f, 1.0f);
    glm::vec3 rVector = glm::normalize(glm::vec3(x, y, 0.0f));

    // Find a view that looks sort of from behind the center of mass
    // towards the origin but somewhat random.
    //
    // This doesn't do quite what I want but it seems to look good
    // enough most of the time

    float rotateAngle = glm::linearRand(0.0f, 360.0f);
    rVector = glm::rotate(centerOfMassDir, rotateAngle, rVector);
    rVector = glm::perp(rVector, zAxis);

    glm::vec3 shifted = glm::normalize(rVector + cmCrossZ);
    glm::fquat startOrient = glm::angleAxis(rangeAngle, shifted);

    this->viewPole.SetOrientation(startOrient);
}

static void nbglSetSceneSettings(scene_t* scene, const VisArgs* args)
{
    OPA_store_int(&scene->blockSimulationOnGraphics, args->blockSimulation);

    /* Must be last thing set */
    OPA_store_int(&scene->attachedPID, (int) getpid());
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

static void nbglPrintGLVersion()
{
    mw_printf("OpenGL %s, GLSL %s\n",
              glGetString(GL_VERSION),
              glGetString(GL_SHADING_LANGUAGE_VERSION));
}

#if USE_GL3W

// return true on error
static bool nbglInitGL3W()
{
    if (gl3wInit())
    {
        mw_printf("Failed to load OpenGL 3\n");
        return true;
    }

    if (!gl3wIsSupported(3, 2))
    {
        mw_printf("OpenGL 3.2 not supported\n");
        return true;
    }

    return false;
}

#else

static bool nbglInitGL3W()
{
    return false;
}

#endif /* USE_GL3W */

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

    if (nbglInitGL3W())
    {
        return 1;
    }

    nbglPrintGLVersion();
    nbglSetSceneSettings(scene, args);
    srand((unsigned int) time(NULL));

    try
    {
        // GL context needs to be open or else destructors will crash
        NBodyGraphics graphicsContext(scene, args);
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

