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
#include "nbody_gl.h"
#include "nbody_gl_util.h"
#include "nbody_graphics.h"
#include "nbody_types.h"
#include "milkyway_cpp_util.h"

#if !BOINC_APPLICATION
  #include <sys/shm.h>
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

#define DELTAXROT 15.0f
#define DELTAYROT 15.0f
#define DELTAZ 3.0f

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
static int monochromatic = FALSE;
static int useGLPoints = FALSE;

static GLUquadricObj* quadratic = NULL;
static scene_t* scene = NULL;


/* print the bindings */
static void print_bindings(FILE* f)
{
    fprintf(f,
            "KEYBINDINGS:\n"
            "  PAGE (UP/DOWN) : zoom (out/in)\n"
            "  ARROW KEYS     : rotate\n"
            "  h              : print this help text\n"
            "  q, ESCAPE      : quit\n"
            "  p              : pause/unpause\n"
            "  SPACE          : step through frames (while paused)\n"
            "  b              : make stars bigger\n"
            "  s              : make stars smaller\n"
            "  m              : use more polygons to render stars\n"
            "  f              : use fewer polygons to render stars\n"
            "  a              : toggle axes\n"
            "  t              : toggle time printout\n"
            "  l              : toggle log printout\n"
            "  <              : slow down animation\n"
            "  >              : speed up animation\n"
            "\n"
            "MOUSEBINDINGS:\n"
            "  DRAG                   : rotate\n"
            "  [SHIFT] DRAG (up/down) : zoom (in/out)\n");
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

#define GALACTIC_RADIUS 15.33f
#define GALACTIC_BULGE_RADIUS 1.5f
#define GALACTIC_DISK_THICKNESS 0.66f

#define AXES_LENGTH (1.25f * GALACTIC_RADIUS / SCALE)

static void drawSimpleGalaxy()
{
    glColor4f(0.643f, 0.706f, 0.867f, 0.85f);
    gluCylinder(quadratic,
                GALACTIC_RADIUS / SCALE,
                GALACTIC_RADIUS / SCALE,
                GALACTIC_DISK_THICKNESS / SCALE,
                100,
                20);

    /* Put caps on the disk */
    glTranslatef(0.0f, 0.5f * GALACTIC_DISK_THICKNESS / SCALE, 0.0f);
    gluDisk(quadratic, 0.0, (GLdouble) GALACTIC_RADIUS / SCALE, 100, 50);
    glTranslatef(0.0f, -0.5f * GALACTIC_DISK_THICKNESS / SCALE, 0.0f);
    gluDisk(quadratic, 0.0, (GLdouble) GALACTIC_RADIUS / SCALE, 100, 50);


    glColor4f(0.98f, 0.835f, 0.714f, 0.85f);
    gluSphere(quadratic, GALACTIC_BULGE_RADIUS / SCALE, 50, 50);
}

static void drawAxes()
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
static void drawPoints()
{
    int i;
    int nbody = scene->nbody;
    FloatPos* r = &scene->r[0];

    if (useGLPoints)
    {
        glPointSize(scene->starsize);  /* Is this actually working? */
        glBegin(GL_POINTS);

        if (monochromatic)
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
            glLoadIdentity();
            glTranslatef(0.0f, 0.0f, scene->z);
            glRotatef(scene->xrot, 1.0f, 0.0f, 0.0f);
            glRotatef(scene->yrot, 0.0f, 1.0f, 0.0f);

            glTranslatef(r[i].x / SCALE, r[i].y / SCALE, r[i].z / SCALE);

            if (monochromatic)
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
        }
    }
}


static void setOrthographicProjection()
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, width, height, 0);
    glMatrixMode(GL_MODELVIEW);
}

static void restorePerspectiveProjection()
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

/* http://www.lighthouse3d.com/tutorials/glut-tutorial/bitmap-fonts-and-orthogonal-projections/ */
static void drawInfo()
{
    char buf[1024];

    setOrthographicProjection();
    glPushMatrix();
    glLoadIdentity();

    glColor3f(1.0f, 1.0f, 1.0f);
    glRasterPos2i(20, 20);

    snprintf(buf, sizeof(buf),
             "Time: %4.3f / %4.3f Gyr (%4.3f %%)\n",
             scene->info.currentTime,
             scene->info.timeEvolve,
             100.0f * scene->info.currentTime / scene->info.timeEvolve
        );
    nbody_glutBitmapStringHelvetica(buf);

    glPopMatrix();
    restorePerspectiveProjection();
}

/* The main drawing function.  Technically there's a race condition
   between drawing and writing new values from the simulation for the
   body positions. I'm lazy and It doesn't need to be displayed 100%
   accurately so this should be good enough.
 */
static void drawGLScene()
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
        glTranslatef(0.0f, 0.0f, scene->z + 0.1f);

        /* rotate view--I know these two rotation matrices don't commute */
        glRotatef(scene->xrot, 1.0f, 0.0f, 0.0f);
        glRotatef(scene->yrot, 0.0f, 1.0f, 0.0f);

        if (scene->drawGalaxy)
        {
            drawSimpleGalaxy();
        }

        if (scene->drawaxes)
        {
            drawAxes();
        }

        drawPoints();
        drawInfo();

        glutSwapBuffers();
    }
}

/* The function called whenever a key is pressed. */
static void keyPressed(unsigned char key, int x, int y)
{
    /* avoid thrashing this call */
    mwMilliSleep(THRASH_SLEEP_INTERVAL);

    if (scene->fullscreen)
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
            print_bindings(stderr);
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

        case 'f':
            scene->ntri /= 2;
            if (scene->ntri < 4)
            {
                scene->ntri = 4;
            }
            scene->changed = TRUE;
            break;

        case 'a':
            scene->drawaxes = !scene->drawaxes;
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
    mwMilliSleep(THRASH_SLEEP_INTERVAL);

    switch (key)
    {
        case GLUT_KEY_PAGE_UP:
            scene->z -= DELTAZ;
            scene->changed = TRUE;
            break;
        case GLUT_KEY_PAGE_DOWN:
            scene->z += DELTAZ;
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
    if (scene->fullscreen)
    {
        mw_finish(EXIT_SUCCESS);
    }

    if (state == GLUT_DOWN)
    {
        xlast = x;
        ylast = y;

        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
        {
            scene->mousemode = 2;
        }
        else
        {
            scene->mousemode = 1;
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

    if (scene->mousemode == 1)
    {
        scene->xrot += ROTSCALE * dy;
        scene->yrot += ROTSCALE * dx;
        scene->changed = TRUE;
    }
    else if (scene->mousemode == 2)
    {
        scene->z -= ZOOMSCALE * dy;
        scene->changed = TRUE;
    }
}

static void passiveMotionFunc(int x, int y)
{
    if (scene->fullscreen)
    {
        mw_finish(EXIT_SUCCESS);
    }
}

static void assignParticleColors(unsigned int nbody)
{
    int i;
    double R, G, B, scale;
    const FloatPos* r = scene->r;

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

static int nbodyInitDrawState()
{
    if (!monochromatic)
    {
        color = (FloatPos*) mwCallocA(scene->nbody, sizeof(FloatPos));
        assignParticleColors(scene->nbody);
    }

    return 0;
}

/* Only one screensaver can be allowed since we do touch some things
 * from the simulation although we probably shouldn't */
static int alreadyAttached()
{
    int oldVal;

    oldVal = nbodyGraphicsTestVal(&scene->attached);
    if (oldVal != 0)
    {
        warn("Screensaver already attached\n");
        return TRUE;
    }

    return FALSE;
}


#if !BOINC_APPLICATION

static key_t key = -1;

static int setShmemKey()
{
    /* TODO: pid of simulation, etc. */
    key = DEFAULT_SHMEM_KEY;

    return 0;
}

int connectSharedScene()
{
    int shmId;
    struct shmid_ds buf;

    setShmemKey();

    shmId = shmget(key, 0, 0);
    if (shmId < 0)
    {
        perror("Error getting shared memory");
        return 1;
    }

    if (shmctl(shmId, IPC_STAT, &buf) < 0)
    {
        perror("Finding shared scene");
        return 0;
    }

    if (buf.shm_nattch > 1)
    {
        perror("Screensaver already attached to segment");
        return 1;
    }
    else if (buf.shm_nattch <= 0)
    {
        perror("Simulation not running");
        return 1;
    }

    scene = (scene_t*) shmat(shmId, NULL, 0);
    if (!scene || scene == (scene_t*) -1)
    {
        perror("Attaching to shared memory");
        return 1;
    }

    return 0;
}

#else

/* Returns TRUE if connection succeeds */
static int attemptConnectSharedScene()
{
    scene = (scene_t*) mw_graphics_get_shmem(NBODY_BIN_NAME);
    if (!scene)
    {
        warn("Failed to connect to shared scene\n");
        return FALSE;
    }

    return TRUE;
}

#define MAX_TRIES 5
#define RETRY_INTERVAL 250

/* In case the main application isn't ready yet, try and wait for a while */
int connectSharedScene()
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

    warn("Could not attach to simulation after %d attempts\n", MAX_TRIES);
    return 1;
}

#endif /* !BOINC_APPLICATION */

static void sceneInit(const VisArgs* args)
{
    width = args->width == 0 ? glutGet(GLUT_SCREEN_WIDTH) / 2 : args->width;
    height = args->height == 0 ? glutGet(GLUT_SCREEN_HEIGHT) / 2 : args->height;
    monochromatic = args->monochrome;
    useGLPoints = !args->notUseGLPoints;

    scene->xrot = -60.0f;
    scene->yrot = -15.0f;
    scene->z = -8.0f;
    scene->starsize = STARSIZE;
    scene->fullscreen = args->fullscreen;
    scene->drawaxes = TRUE;
    scene->ntri = NTRI;
    scene->paused = FALSE;
    scene->step = FALSE;
    scene->mousemode = 1;
    scene->changed = FALSE;
    scene->dt = 300;
    scene->usleepdt = USLEEPDT;
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
void nbodyGLCleanup()
{
    mwFreeA(color);
    color = NULL;
}

