/*
 * Copyright (c) 2002-2006 John M. Fregeau, Richard Campbell, Jeff Molofee
 * Copyright (c) 2011 Matthew Arsenault
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

/*
  Derived from glstarview-0.6

  This code was created by Jeff Molofee '99 (ported to Linux/GLUT by Richard Campbell '99)

  If you've found this code useful, please let me know.

  Visit me at www.demonews.com/hosted/nehe
  (email Richard Campbell at ulmont@bellsouth.net)
*/

#include "milkyway_util.h"
#include "milkyway_boinc_util.h"
#include "nbody_gl.h"
#include "nbody_gl_util.h"
#include "nbody_graphics.h"
#include "nbody_types.h"


#if !BOINC_APPLICATION
  #include <sys/mman.h>
  #include <sys/stat.h>
  #include <fcntl.h>
  #include <errno.h>
#endif /* !BOINC_APPLICATION */



/* ASCII keycodes */
#define ESCAPE 27
#define PAGE_UP 73
#define PAGE_DOWN 81
#define UP_ARROW 72
#define DOWN_ARROW 80
#define LEFT_ARROW 75
#define RIGHT_ARROW 77

/* other stuff */
#define STARSIZE 0.01f
#define NTRI 8
#define ORBIT_TRACE_POINT_SIZE 0.03f

#define DELTAXROT 5.0f
#define DELTAYROT 5.0f
#define DELTAR 3.0f

#define SCALE 15.0

#define ROTSCALE 1.0
#define ZOOMSCALE 0.2

/* microseconds */
#define USLEEPDT 10000.0

/* Milliseconds */
#define THRASH_SLEEP_INTERVAL 15


static FloatPos* color = NULL;

static const FloatPos white = { 1.0f, 1.0f, 1.0f, 0 };
static const FloatPos grey = { 0.3294f, 0.3294f, 0.3294f, 1 };

static int window = 0;
static int xlast = 0, ylast = 0;
static int width = 0, height = 0;

static GLUquadricObj* quadratic = NULL;
static scene_t* scene = NULL;
static GLboolean ownScene = GL_FALSE;  /* Is this scene owned by the graphics or shared memory? */
static dsfmt_t rndState;

/* Max/min time in seconds between changing directions when randomly moving */
#define MIN_CHANGE_INTERVAL 10
#define MAX_CHANGE_INTERVAL 60

#define MAX_FLOAT_RATE 0.3f
#define MIN_FLOAT_RATE 0.3f

static const char* helpBindings =
    "KEYBINDINGS:\n"
    "  PAGE (UP/DOWN) : zoom (out/in)\n"
    "  ARROW KEYS     : rotate\n"
    "  h or ?         : display this help text\n"
    "  q, ESCAPE      : quit\n"
    "  p              : pause/unpause\n"
    "  SPACE          : step through frames (while paused)\n"
    "  a              : (axes) toggle axes\n"
    "  o              : (origin) toggle origin of visualization\n"
    "  r              : (rotate) when idle, float view around randomly\n"
    "  t              : (trace) toggle orbit trace\n"
    "  i              : (info) toggle info printout\n"
    "  n              : toggle showing individual particles (can be used to only view CM orbit)\n"
    "  c              : (color) toggle particle color scheme\n"
    "  l              : toggle using GL points (default, faster) or spheres\n"
    "  b              : (bigger) make stars bigger\n"
    "  s              : (smaller) make stars smaller\n"
    "  m              : (more) use more polygons to render stars\n"
    "  f              : (fewer) use fewer polygons to render stars\n"
    "  <              : slow down animation\n"
    "  >              : speed up animation\n"
    "\n"
    "MOUSEBINDINGS:\n"
    "  DRAG                        : rotate\n"
    "  [SHIFT] DRAG (up/down)      : zoom (in/out)\n"
    "  [RIGHTCLICK] DRAG (up/down) : zoom (in/out)\n";

int nbglLoadStaticSceneFromFile(const char* filename)
{
    FILE* f;
    size_t lnCount; /* ~= nbody */
    size_t line = 0;
    int nbody = 0;
    char lnBuf[4096];
    int rc = 0;
    int ignore;
    double x, y, z;
    double vx, vy, vz;
    double lambda;
    FloatPos* r;

    f = fopen(filename, "r");
    if (!f)
    {
        mwPerror("Failed to open file '%s'", filename);
        return 1;
    }

    lnCount = mwCountLinesInFile(f);
    if (lnCount == 0)
    {
        mw_printf("Error counting lines from file '%s'\n", filename);
        fclose(f);
        return 1;
    }

    scene = mwCalloc(sizeof(scene_t) + lnCount * sizeof(FloatPos), sizeof(char));
    r = scene->rTrace;

    /* Skip the 1st line with the # comment */
    fgets(lnBuf, sizeof(lnBuf), f);

    while (rc != EOF)
    {
        ++line;
        rc = fscanf(f,
                    "%d , %lf , %lf , %lf , %lf , %lf , %lf ",
                    &ignore,
                    &x, &y, &z,
                    &vx, &vy, &vz);
        if (rc == 7)
        {
            ++nbody;
            /* May or may not be there */
            rc = fscanf(f, " , %lf \n", &lambda);
            if (rc != 1)
            {
                fscanf(f, " \n");
            }
            assert(line < lnCount);

            r[line].x = (float) x;
            r[line].y = (float) y;
            r[line].z = (float) z;
        }
        else if (rc != EOF)
        {
            mw_printf("Error reading '%s' at line "ZU"\n", filename, line);
        }
    }

    if (rc != EOF)
    {
        fclose(f);
        free(scene);
        scene = NULL;
        return 1;
    }

    scene->nbody = nbody;

    if (fclose(f))
    {
        mwPerror("Failed to close file '%s'", filename);
    }

    return 0;
}



/* A general OpenGL initialization function.  Sets all of the initial parameters. */
static void initGL(int w, int h)
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    glClearDepth(1.0);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST | GL_POINT_SMOOTH | GL_BLEND | GL_ALPHA_TEST);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glShadeModel(GL_SMOOTH);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    quadratic = gluNewQuadric();
    gluQuadricNormals(quadratic, GLU_SMOOTH);
    gluQuadricTexture(quadratic, GL_TRUE);

    gluPerspective(45.0f, (GLfloat)w / (GLfloat)h, 0.1f, 1000.0f);

    glMatrixMode(GL_MODELVIEW);
}

/* The function called when our window is resized */
static void resizeGLScene(int w, int h)
{
    /* Prevent A Divide By Zero If The Window Is Too Small */
    if (h == 0)
    {
        h = 1;
    }

    width = w;
    height = h;

    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(45.0f, (GLfloat)w / (GLfloat)h, 0.1f, 1000.0f);
    glMatrixMode(GL_MODELVIEW);

    scene->changed = TRUE;
}

#define MILKYWAY_RADIUS 15.33f
#define MILKYWAY_BULGE_RADIUS 1.5f
#define MILKYWAY_DISK_THICKNESS 0.66f


#define AXES_LENGTH (1.25f * MILKYWAY_RADIUS / SCALE)

static void drawSimpleGalaxy(float galacticRadius, float galacticBulgeRadius, float galacticDiskThickness)
{
    glColor4f(0.643f, 0.706f, 0.867f, 0.85f);
    gluCylinder(quadratic,
                galacticRadius / SCALE,
                galacticRadius / SCALE,
                galacticDiskThickness / SCALE,
                100,
                20);

    /* Put caps on the disk */
    glTranslatef(0.0f, 0.5f * galacticDiskThickness / SCALE, 0.0f);
    gluDisk(quadratic, 0.0, (GLdouble) galacticRadius / SCALE, 100, 50);

    glTranslatef(0.0f, -0.5f * galacticDiskThickness / SCALE, 0.0f);
    gluDisk(quadratic, 0.0, (GLdouble) galacticRadius / SCALE, 100, 50);

    glColor4f(0.98f, 0.835f, 0.714f, 0.85f);
    gluSphere(quadratic, galacticBulgeRadius / SCALE, 50, 50);
}

static void drawAxes(void)
{
    glColor3f(1.0, 0.0, 0.0);  /* x axis */
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(AXES_LENGTH, 0.0, 0.0);
    glEnd();

    glColor3f(0.0, 1.0, 0.0);  /* y axis */
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, AXES_LENGTH, 0.0);
    glEnd();

    glColor3f(0.0, 0.0, 1.0); /* z axis */
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, AXES_LENGTH);
    glEnd();
}

static void drawOrbitTrace(void)
{
    int i;
    const int n = scene->currentTracePoint;
    const FloatPos* trace = &scene->orbitTrace[0];

    glColor3f(0.0f, 1.0f, 0.0f); /* Draw starting point as different color */
    for (i = 0; i < n; ++i)
    {
        glPushMatrix();

        glTranslatef(trace[i].x / SCALE, trace[i].y / SCALE, trace[i].z / SCALE);

        gluSphere(quadratic, ORBIT_TRACE_POINT_SIZE, scene->ntri, scene->ntri);
        glColor3f(white.x, white.y, white.z);

        glPopMatrix();
    }
}

static inline void setIgnoredColor(int ignore)
{
    if (ignore)
    {
        glColor3f(grey.x, grey.y, grey.z);
    }
    else
    {
        glColor3f(white.x, white.y, white.z);
    }
}

/* draw stars */
static void drawParticles(void)
{
    int i;
    int nbody = scene->nbody;
    const FloatPos* r = &scene->rTrace[0];

    if (scene->useGLPoints)
    {
        glPointSize(scene->starsize);  /* Is this actually working? */
        glBegin(GL_POINTS);

        if (scene->monochromatic)
        {
            for (i = 0; i < nbody; ++i)
            {
                setIgnoredColor(r[i].ignore);
                glVertex3f(r[i].x / SCALE, r[i].y / SCALE, r[i].z / SCALE);
            }
        }
        else
        {
            for (i = 0; i < nbody; ++i)
            {
                glColor3f(color[i].x, color[i].y, color[i].z);
                glVertex3f(r[i].x / SCALE, r[i].y / SCALE, r[i].z / SCALE);
            }
        }

        glEnd();
    }
    else
    {
        /* FIXME: Use modern OpenGL stuff */
        for (i = 0; i < nbody; ++i)
        {
            glPushMatrix();
            glTranslatef(r[i].x / SCALE, r[i].y / SCALE, r[i].z / SCALE);

            if (scene->monochromatic)
            {
                setIgnoredColor(r[i].ignore);
            }
            else
            {
                glColor3f(color[i].x, color[i].y, color[i].z);
            }
            /* glutSolidSphere(scene->starsize, scene->ntri, scene->ntri); */
            gluSphere(quadratic, scene->starsize, scene->ntri, scene->ntri);
            /* glTranslatef(-r[i].x/SCALE, -r[i].y/SCALE, -r[i].z/SCALE); */

            glPopMatrix();
        }
    }
}


static void setOrthographicProjection(void)
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, width, height, 0);
    glMatrixMode(GL_MODELVIEW);
}

static void restorePerspectiveProjection(void)
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

/* http://www.lighthouse3d.com/tutorials/glut-tutorial/bitmap-fonts-and-orthogonal-projections/ */
static void drawWords(const char* str, int x, int y)
{
    setOrthographicProjection();
    glPushMatrix();
    glLoadIdentity();

    glColor3f(1.0f, 1.0f, 1.0f);
    glRasterPos2i(x, y);
    nbody_glutBitmapStringHelvetica((unsigned char*) str);

    glPopMatrix();
    restorePerspectiveProjection();
}

static void drawHelp()
{
    /* FIXME: Things don't line up since not using monospace font */
    drawWords(helpBindings, width - 500, 20);
}

static void drawInfo(void)
{
    char buf[1024];

    snprintf(buf, sizeof(buf),
             "Time: %4.3f / %4.3f Gyr (%4.3f %%)\n",
             scene->info.currentTime,
             scene->info.timeEvolve,
             100.0f * scene->info.currentTime / scene->info.timeEvolve
        );

    drawWords(buf, 20, 20);
}

/* Center of mass to center of galaxy */
static float centerOfMassDistance(void)
{
    return sqrtf(sqr(scene->rootCenterOfMass[0]) + sqr(scene->rootCenterOfMass[1]) + sqr(scene->rootCenterOfMass[2]));
}

/* The main drawing function.  Technically there's a race condition
   between drawing and writing new values from the simulation for the
   body positions. I'm lazy and It doesn't need to be displayed 100%
   accurately so this should be good enough.
 */
static void drawGLScene(void)
{
    if (!scene)
        return;

    mwMicroSleep(scene->usleepdt);

    /* draw scene if necessary */
    if (scene->changed)
    {
        nbodyGraphicsSetOn(&scene->attached);
        scene->changed = FALSE;

        /* erase scene */
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

        /* get ready to render 3D objects */
        glLoadIdentity();

        /* Spin around mostly randomly when fullscreen.
           This should be smarter than it is, and always try to look at something interesting
        */
        if (scene->floatMode)
        {
            /* Select a random angle that is a multiple of 45 degrees */
            static float floatRate = 0.0f;
            static float thetaFloatRate = 0.1f; /* Movement rates when moving on its own */
            static float phiFloatRate = 0.1f;
            static float zFloatRate = 0.1f; /* Zoom */

            static float destR = 0.0f;
            static time_t changeInterval = 0;
            static time_t lastChange = 0;
            time_t curT;

            curT = time(NULL);
            if (curT - lastChange >= changeInterval)
            {
                float rCM = centerOfMassDistance();
                lastChange = curT;
                changeInterval = rand() % MAX_CHANGE_INTERVAL; /* Move in one direction between 10 and 60 seconds */
                if (changeInterval < MIN_CHANGE_INTERVAL)
                    changeInterval = MIN_CHANGE_INTERVAL;

                floatRate = mwXrandom(&rndState, MIN_FLOAT_RATE, MAX_FLOAT_RATE);

                thetaFloatRate = mwXrandom(&rndState, 0.0, floatRate);
                phiFloatRate = sqrtf(sqr(floatRate) - sqr(thetaFloatRate));

                destR = mwXrandom(&rndState, 1.1f * rCM, 1.5f * rCM);
            }

            if (fabsf(scene->r) - destR >= 0.1f) /* Not at correct altitude */
            {
                if (fabsf(scene->r) > destR) /* Farther than destination altitude */
                {
                    if (scene->r < 0.0)
                        scene->r += zFloatRate;
                    else
                        scene->r -= zFloatRate;
                }
                else if (fabsf(scene->r) < destR) /* Closer than destination altitude */
                {
                    if (scene->r > 0.0)
                        scene->r += zFloatRate;
                    else
                        scene->r -= zFloatRate;
                }
            }

            scene->xrot += thetaFloatRate;
            scene->yrot += phiFloatRate;
        }

        glTranslatef(0.0f, 0.0f, scene->r + 0.1f);

        /* rotate view--I know these two rotation matrices don't commute */
        glRotatef(scene->yrot, 0.0f, 1.0f, 0.0f);
        glRotatef(scene->xrot, 1.0f, 0.0f, 0.0f);

        /* Focus on the center of mass */
        if (scene->cmCentered)
        {
            glTranslatef(-scene->rootCenterOfMass[0] / SCALE,
                         -scene->rootCenterOfMass[1] / SCALE,
                         -scene->rootCenterOfMass[2] / SCALE);
        }

        if (scene->drawGalaxy)
        {
            drawSimpleGalaxy(MILKYWAY_RADIUS, MILKYWAY_BULGE_RADIUS, MILKYWAY_DISK_THICKNESS);
        }

        if (scene->drawAxes)
        {
            drawAxes();
        }

        if (scene->drawParticles)
        {
            drawParticles();
        }

        if (scene->drawOrbitTrace)
        {
            drawOrbitTrace();
        }

        if (scene->drawHelp)
        {
            drawHelp();
        }

        if (scene->drawInfo)
        {
            drawInfo();
        }

        glutSwapBuffers();
    }
}

/* The function called whenever a key is pressed. */
static void keyPressed(unsigned char key, int x, int y)
{
    (void) x, (void) y;
    /* avoid thrashing this call */
    mwMilliSleep(THRASH_SLEEP_INTERVAL);

    if (scene->screensaverMode)
    {
        mw_finish(EXIT_SUCCESS);
    }

    switch (key)
    {
        case ESCAPE:
        case 'q':
            glutDestroyWindow(window);
            mw_finish(EXIT_SUCCESS);

        case 'h':
        case '?':
            scene->drawHelp = !scene->drawHelp;
            scene->changed = TRUE;
            break;

        case 'o': /* Toggle camera following CM or on milkyway center */
            scene->cmCentered = !scene->cmCentered;
            break;

        case 'r': /* Toggle floating */
            scene->floatMode = !scene->floatMode;
            break;

        case 'p':
            scene->paused = !scene->paused;
            break;

        case ' ':
            scene->step = TRUE;
            break;

        case 'b':
            scene->starsize *= 1.1;
            if (scene->starsize > 100.0)
            {
                scene->starsize = 100.0;
            }
            scene->changed = TRUE;
            break;

        case 's':
            scene->starsize *= 0.9;
            if (scene->starsize < 1.0e-3)
            {
                scene->starsize = 1.0e-3;
            }
            scene->changed = TRUE;
            break;

        case 'm':
            scene->ntri *= 2;
            if (scene->ntri > 24)
            {
                scene->ntri = 24;
            }
            scene->changed = TRUE;
            break;

        case 'n':
            scene->drawParticles = !scene->drawParticles;
            scene->changed = TRUE;
            break;

        case 'l':
            scene->useGLPoints = !scene->useGLPoints;
            scene->changed = TRUE;
            break;

        case 'f':
            scene->ntri /= 2;
            if (scene->ntri < 4)
            {
                scene->ntri = 4;
            }
            scene->changed = TRUE;
            break;

        case 'a':
            scene->drawAxes = !scene->drawAxes;
            scene->changed = TRUE;
            break;

        case 't':
            scene->drawOrbitTrace = !scene->drawOrbitTrace;
            scene->changed = TRUE;
            break;

        case 'i':
            scene->drawInfo = !scene->drawInfo;
            scene->changed = TRUE;
            break;

        case 'c':
            scene->monochromatic = !scene->monochromatic;
            scene->changed = TRUE;
            break;

        case '>':
            scene->dt /= 2.0;
            if (scene->dt < 10.0)
            {
                scene->dt = 10.0;
            }
            break;

        case '<':
            scene->dt *= 2.0;
            if (scene->dt > 1.0e6)
            {
                scene->dt = 1.0e6;
            }
            break;

        default:
            break;
    }
}

/* The function called whenever a normal key is pressed. */
static void specialKeyPressed(int key, int x, int y)
{
    (void) x, (void) y;
    mwMilliSleep(THRASH_SLEEP_INTERVAL);

    switch (key)
    {
        case GLUT_KEY_PAGE_UP:
            scene->r -= DELTAR;
            scene->changed = TRUE;
            break;
        case GLUT_KEY_PAGE_DOWN:
            scene->r += DELTAR;
            scene->changed = TRUE;
            break;
        case GLUT_KEY_UP:
            scene->xrot -= DELTAXROT;
            scene->changed = TRUE;
            break;
        case GLUT_KEY_DOWN:
            scene->xrot += DELTAXROT;
            scene->changed = TRUE;
            break;
        case GLUT_KEY_LEFT:
            scene->yrot -= DELTAYROT;
            scene->changed = TRUE;
            break;
        case GLUT_KEY_RIGHT:
            scene->yrot += DELTAYROT;
            scene->changed = TRUE;
            break;
        default:
            break;
    }
}

static void mouseFunc(int button, int state, int x, int y)
{
    int mods;

    if (scene->screensaverMode)
    {
        mw_finish(EXIT_SUCCESS);
    }

    if (state == GLUT_DOWN)
    {
        xlast = x;
        ylast = y;
        mods = glutGetModifiers();

        if ((button == GLUT_RIGHT_BUTTON) || (mods == GLUT_ACTIVE_SHIFT))
        {
            scene->mouseMode = MOUSE_MODE_ZOOM;
        }
        else
        {
            scene->mouseMode = MOUSE_MODE_MOVE;
        }
    }
}

static void motionFunc(int x, int y)
{
    double dx, dy;

    dx = (double) x - (double) xlast;
    dy = (double) y - (double) ylast;

    xlast = x;
    ylast = y;

    switch (scene->mouseMode)
    {
        case MOUSE_MODE_MOVE:
            scene->xrot += ROTSCALE * dy;
            scene->yrot += ROTSCALE * dx;
            scene->changed = TRUE;
            break;

        case MOUSE_MODE_ZOOM:
            scene->r -= ZOOMSCALE * dy;
            scene->changed = TRUE;
            break;

        case MOUSE_MODE_NONE:
            break;

        default:
            mw_unreachable();
    }
}

static void passiveMotionFunc(int x, int y)
{
    (void) x, (void) y;

    if (scene->screensaverMode)
    {
        mw_finish(EXIT_SUCCESS);
    }
}

static void assignParticleColors(int nbody)
{
    int i;
    double R, G, B, scale;
    const FloatPos* r = scene->rTrace;

    /* assign random particle colors */
    srand((unsigned int) time(NULL));

    for (i = 0; i < nbody; ++i)
    {
        if (r[i].ignore)
        {
            R = grey.x;  /* TODO: Random greyish color? */
            G = grey.z;
            B = grey.y;
        }
        else
        {
            R = ((double) rand()) / ((double) RAND_MAX);
            G = ((double) rand()) / ((double) RAND_MAX) * (1.0 - R);
            B = 1.0 - R - G;
        }

        if (R >= G && R >= B)
        {
            scale = 1.0 + ((double) rand()) / ((double) RAND_MAX) * (MIN(2.0, 1.0 / R) - 1.0);
        }
        else if (G >= R && G >= B)
        {
            scale = 1.0 + ((double) rand()) / ((double) RAND_MAX) * (MIN(2.0, 1.0 / G) - 1.0);
        }
        else
        {
            scale = 1.0 + ((double) rand()) / ((double) RAND_MAX) * (MIN(2.0, 1.0 / B) - 1.0);
        }

        color[i].ignore = r[i].ignore;
        color[i].x = R * scale;
        color[i].y = G * scale;
        color[i].z = B * scale;
    }
}

static int nbodyInitDrawState(void)
{
    color = (FloatPos*) mwCallocA(scene->nbody, sizeof(FloatPos));
    assignParticleColors(scene->nbody);

    return 0;
}

/* Only one screensaver can be allowed since we do touch some things
 * from the simulation although we probably shouldn't */
static int alreadyAttached(void)
{
    int oldVal;

    oldVal = nbodyGraphicsTestVal(&scene->attached);
    if (oldVal != 0)
    {
        mw_printf("Screensaver already attached\n");
        return TRUE;
    }

    return FALSE;
}


#if !BOINC_APPLICATION

int connectSharedScene(int instanceId)
{
    int shmId;
    const int mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
    struct stat sb;
    char name[128];

    if (snprintf(name, sizeof(name), "/milkyway_nbody_%d", instanceId) == sizeof(name))
        mw_panic("name buffer too small for shared memory name\n");

    shmId = shm_open(name, O_RDWR, mode);
    if (shmId < 0)
    {
        mwPerror("Error getting shared memory");
        return errno;
    }

    if (fstat(shmId, &sb) < 0)
    {
        mwPerror("shmem fstat");
        return errno;
    }

    scene = (scene_t*) mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, shmId, 0);
    if (scene == MAP_FAILED)
    {
        mwPerror("mmap: Failed to mmap shared memory");
        if (shm_unlink(name) < 0)
        {
            mwPerror("Unlink shared memory");
        }

        return 1;
    }

    return 0;
}

#else

/* Returns TRUE if connection succeeds */
static int attemptConnectSharedScene(void)
{
    scene = (scene_t*) mw_graphics_get_shmem(NBODY_BIN_NAME);
    if (!scene)
    {
        mw_printf("Failed to connect to shared scene\n");
        return FALSE;
    }

    return TRUE;
}

#define MAX_TRIES 5
#define RETRY_INTERVAL 250

/* In case the main application isn't ready yet, try and wait for a while */
int connectSharedScene(int instanceId)
{
    int tries = 0;

    while (tries < MAX_TRIES)
    {
        if (attemptConnectSharedScene())
        {
            return alreadyAttached(); /* Error if something already attached */
        }

        mwMilliSleep(RETRY_INTERVAL);
        ++tries;
    }

    mw_printf("Could not attach to simulation after %d attempts\n", MAX_TRIES);
    return 1;
}

#endif /* !BOINC_APPLICATION */

static void sceneInit(const VisArgs* args)
{
    width = args->width == 0 ? glutGet(GLUT_SCREEN_WIDTH) / 2 : args->width;
    height = args->height == 0 ? glutGet(GLUT_SCREEN_HEIGHT) / 2 : args->height;
    scene->monochromatic = args->monochrome;
    scene->useGLPoints = !args->notUseGLPoints;

    scene->cmCentered = !args->originCenter;
    scene->floatMode = !args->noFloat;

#if 0
    phi = atanf(scene->rootCenterOfMass[1] / scene->rootCenterOfMass[0]);
    theta = atanf(sqrtf(sqr(scene->rootCenterOfMass[0]) + sqr(scene->rootCenterOfMass[1])) / scene->rootCenterOfMass[2]);
    rho = centerOfMassDistance();

    /* This doesn't actually quite work as intended but that's OK */
    scene->xrot = r2d(theta);
    scene->yrot = -r2d(phi);
    scene->r = -2.0f * rho / SCALE;

    if (scene->r - scene->rootCenterOfMass[2] < 0.0f)
    {
        scene->r = -5.0f * rho / SCALE;
    }
#endif

    scene->xrot = -60.0f;
    scene->yrot = -15.0f;
    scene->r = -8.0f;

    scene->starsize = STARSIZE;
    scene->fullscreen = args->fullscreen || args->plainFullscreen;
    scene->screensaverMode = scene->fullscreen && !args->plainFullscreen;
    scene->drawInfo = TRUE;
    scene->drawHelp = FALSE;
    scene->drawAxes = TRUE;
    scene->drawParticles = TRUE;
    scene->drawOrbitTrace = TRUE;
    scene->ntri = NTRI;
    scene->paused = FALSE;
    scene->step = FALSE;
    scene->mouseMode = MOUSE_MODE_MOVE;
    scene->changed = FALSE;
    scene->dt = 300;
    scene->usleepdt = USLEEPDT;

    dsfmt_init_gen_rand(&rndState, (uint32_t) time(NULL));
}

int nbodyGLSetup(const VisArgs* args)
{
    sceneInit(args);
    nbodyInitDrawState();

    /* prepare rendering */
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutInitWindowSize(width, height);
    glutInitWindowPosition(0, 0);

    window = glutCreateWindow("Milkyway@Home N-body");
    glutDisplayFunc(drawGLScene);

    if (scene->fullscreen)
    {
        glutFullScreen();
        glutPassiveMotionFunc(passiveMotionFunc);
    }

    glutIdleFunc(drawGLScene);
    glutReshapeFunc(resizeGLScene);
    glutKeyboardFunc(keyPressed);
    glutSpecialFunc(specialKeyPressed);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);

    initGL(width, height);

    return 0;
}

/* glut main loop never quits, so actually freeing these is a bad idea */
void nbodyGLCleanup(void)
{
    mwFreeA(color);
    color = NULL;

    if (ownScene)
    {
        free(scene);
        scene = NULL;
        ownScene = GL_FALSE;
    }
}

int checkConnectedVersion(void)
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

