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
#include "nbody_types.h"
#include "nbody_priv.h"

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

#define USE_GL_POINTS 0

/* the scene structure */
typedef struct
{
    int nstar;
    GLfloat z;
    GLfloat xrot;
    GLfloat yrot;
    GLfloat starsize;
    int fullscreen;
    int drawaxes;
    int ntri;
    int paused;
    int step;
    int mousemode;
    int changed;
    double t;
    char tstring[128];
    char log[3500];
    double dt;
} scene_t;

static const NBodyState* _drawState = NULL;

/* Set this when colors and everything ready. I'm pretty sure we don't
 * need to lock around this. We just need to wait for it to be nonzero
 * in the drawing thread. */
static volatile int drawStateReady = FALSE;

typedef struct
{
    float x, y, z;
} FloatPos;


static FloatPos* r = NULL;
static FloatPos* color = NULL;


/* So many global variables... */
static int window, xlast = 0, ylast = 0, debug = 0, width = 1024, height = 768;
static const double usleepdt = 100.0;
static GLUquadricObj* quadratic;
static scene_t scene;


/* print the bindings */
static void print_bindings(FILE* stream)
{
    fprintf(stream, "KEYBINDINGS:\n");
    fprintf(stream, "  PAGE (UP/DOWN) : zoom (out/in)\n");
    fprintf(stream, "  ARROW KEYS     : rotate\n");
    fprintf(stream, "  h              : print this help text\n");
    fprintf(stream, "  q, ESCAPE      : quit\n");
    fprintf(stream, "  p              : pause/unpause\n");
    fprintf(stream, "  SPACE          : step through frames (while paused)\n");
    fprintf(stream, "  b              : make stars bigger\n");
    fprintf(stream, "  s              : make stars smaller\n");
    fprintf(stream, "  m              : use more polygons to render stars\n");
    fprintf(stream, "  f              : use fewer polygons to render stars\n");
    fprintf(stream, "  a              : toggle axes\n");
    fprintf(stream, "  t              : toggle time printout\n");
    fprintf(stream, "  l              : toggle log printout\n");
    fprintf(stream, "  <              : slow down animation\n");
    fprintf(stream, "  >              : speed up animation\n");
    fprintf(stream, "\n");
    fprintf(stream, "MOUSEBINDINGS:\n");
    fprintf(stream, "  DRAG                   : rotate\n");
    fprintf(stream, "  [SHIFT] DRAG (up/down) : zoom (in/out)\n");
}

/* A general OpenGL initialization function.  Sets all of the initial parameters. */
static void initGL(int w, int h)
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    glClearDepth(1.0);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST | GL_POINT_SMOOTH);

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

    scene.changed = 1;
}

void updateDisplayedBodies()
{
    const Body* b;
    const NBodyState* st = _drawState;
    static double usleepcount = 0.0;
    unsigned int i = 0;
    const unsigned int nbody = st->nbody;

    usleepcount += usleepdt;

    /* read data if not paused */
    if (usleepcount >= scene.dt && (scene.paused == 0 || scene.step == 1))
    {
        usleepcount = 0.0;
        scene.step = 0;

      #ifdef _OPENMP
        #pragma omp parallel for private(i, b) schedule(static)
      #endif
        for (i = 0; i < nbody; ++i)
        {
            b = &st->bodytab[i];
            r[i].x = (float) X(Pos(b));
            r[i].y = (float) Y(Pos(b));
            r[i].z = (float) Z(Pos(b));
        }
    }

    scene.changed = 1;
}

/* The main drawing function */
static void drawGLScene()
{
    unsigned int i = 0;
    const NBodyState* st = _drawState;

    usleep(usleepdt);

    /* draw scene if necessary */
    if (scene.changed)
    {
        scene.changed = 0;

        /* erase scene */
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        glLoadIdentity();

        /* get ready to render 3D objects */
        glTranslatef(0.0f, 0.0f, scene.z + 0.1);

        /* rotate view--I know these two rotation matrices don't commute */
        glRotatef(scene.xrot, 1.0f, 0.0f, 0.0f);
        glRotatef(scene.yrot, 0.0f, 1.0f, 0.0f);

        /* draw axes */
        if (scene.drawaxes)
        {
            glColor3f(1.0, 0.0, 0.0);
            glBegin(GL_LINES);
            glVertex3f(0.0, 0.0, 0.0);
            glVertex3f(1.0, 0.0, 0.0);
            glEnd();
            glColor3f(0.0, 1.0, 0.0);
            glBegin(GL_LINES);
            glVertex3f(0.0, 0.0, 0.0);
            glVertex3f(0.0, 1.0, 0.0);
            glEnd();
            glColor3f(0.0, 0.0, 1.0);
            glBegin(GL_LINES);
            glVertex3f(0.0, 0.0, 0.0);
            glVertex3f(0.0, 0.0, 1.0);
            glEnd();
        }

        /* draw stars */
    #if !USE_GL_POINTS
        /* FIXME: Use modern OpenGL stuff */
        for (i = 0; i < st->nbody; ++i)
        {

            glLoadIdentity();
            glTranslatef(0.0f, 0.0f, scene.z);
            glRotatef(scene.xrot, 1.0f, 0.0f, 0.0f);
            glRotatef(scene.yrot, 0.0f, 1.0f, 0.0f);

            glTranslatef(r[i].x / SCALE, r[i].y / SCALE, r[i].z / SCALE);
            glColor3f(color[i].x, color[i].y, color[i].z);
            /* glutSolidSphere(scene.starsize, scene.ntri, scene.ntri); */
            gluSphere(quadratic, scene.starsize, scene.ntri, scene.ntri);
            /* glTranslatef(-r[i].x/SCALE, -r[i].y/SCALE, -r[i].z/SCALE); */
        }
    #else
        glPointSize(scene.starsize);  /* Is this actually working? */
        glBegin(GL_POINTS);
        glColor3f(1.0f, 1.0f, 1.0f);
        for (i = 0; i < st->nbody; ++i)
        {
            glVertex3f(r[i].x / SCALE, r[i].y / SCALE, r[i].z / SCALE);
        }
        glEnd();
    #endif /* !USE_GL_POINTS */

        glutSwapBuffers();
    }
}

/* The function called whenever a key is pressed. */
static void keyPressed(unsigned char key, int x, int y)
{
    /* avoid thrashing this call */
    usleep(100);

    if (key == ESCAPE || key == 'q')
    {
        glutDestroyWindow(window);
        mw_finish(EXIT_SUCCESS);
    }
    else if (key == 'h')
    {
        print_bindings(stderr);
    }
    else if (key == 'p')
    {
        scene.paused = !scene.paused;
    }
    else if (key == ' ')
    {
        scene.step = 1;
    }
    else if (key == 'b')
    {
        scene.starsize *= 1.5;
        if (scene.starsize > 100.0)
        {
            scene.starsize = 100.0;
        }
        scene.changed = 1;
    }
    else if (key == 's')
    {
        scene.starsize /= 1.5;
        if (scene.starsize < 1.0e-3)
        {
            scene.starsize = 1.0e-3;
        }
        scene.changed = 1;
    }
    else if (key == 'm')
    {
        scene.ntri *= 2;
        if (scene.ntri > 24)
        {
            scene.ntri = 24;
        }
        scene.changed = 1;
    }
    else if (key == 'f')
    {
        scene.ntri /= 2;
        if (scene.ntri < 4)
        {
            scene.ntri = 4;
        }
        scene.changed = 1;
    }
    else if (key == 'a')
    {
        scene.drawaxes = !scene.drawaxes;
        scene.changed = 1;
    }
    else if (key == '>')
    {
        scene.dt /= 2.0;
        if (scene.dt < 10.0)
        {
            scene.dt = 10.0;
        }
    }
    else if (key == '<')
    {
        scene.dt *= 2.0;
        if (scene.dt > 1.0e6)
        {
            scene.dt = 1.0e6;
        }
    }
}

/* The function called whenever a normal key is pressed. */
static void specialKeyPressed(int key, int x, int y)
{
    /* avoid thrashing this procedure */
    usleep(100);

    switch (key)
    {
        case GLUT_KEY_PAGE_UP:
            scene.z -= DELTAZ;
            scene.changed = 1;
            break;
        case GLUT_KEY_PAGE_DOWN:
            scene.z += DELTAZ;
            scene.changed = 1;
            break;
        case GLUT_KEY_UP:
            scene.xrot -= DELTAXROT;
            scene.changed = 1;
            break;
        case GLUT_KEY_DOWN:
            scene.xrot += DELTAXROT;
            scene.changed = 1;
            break;
        case GLUT_KEY_LEFT:
            scene.yrot -= DELTAYROT;
            scene.changed = 1;
            break;
        case GLUT_KEY_RIGHT:
            scene.yrot += DELTAYROT;
            scene.changed = 1;
            break;
        default:
            break;
    }
}

static void mouseFunc(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN)
    {
        xlast = x;
        ylast = y;

        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
        {
            scene.mousemode = 2;
        }
        else
        {
            scene.mousemode = 1;
        }
    }
}

static void motionFunc(int x, int y)
{
    double dx, dy;

    dx = ((double) x) - ((double) xlast);
    dy = ((double) y) - ((double) ylast);

    xlast = x;
    ylast = y;

    if (scene.mousemode == 1)
    {
        scene.xrot += ROTSCALE * dy;
        scene.yrot += ROTSCALE * dx;
        scene.changed = 1;
    }
    else if (scene.mousemode == 2)
    {
        scene.z -= ZOOMSCALE * dy;
        scene.changed = 1;
    }
}

static void assignParticleColors(unsigned int nbody)
{
    int i;
    double R, G, B, scale;

    /* assign particle colors */
    srand((unsigned int) time(NULL));

    color[0].x = 1.0;
    color[0].y = 0.0;
    color[0].z = 0.0;

    color[1].x = 0.0;
    color[1].y = 1.0;
    color[1].z = 0.0;

    color[2].x = 0.0;
    color[2].y = 0.0;
    color[2].z = 1.0;

    color[3].x = 1.0;
    color[3].y = 1.0;
    color[3].z = 0.0;

    color[4].x = 1.0;
    color[4].y = 0.0;
    color[4].z = 1.0;

    color[5].x = 0.0;
    color[5].y = 1.0;
    color[5].z = 1.0;

    for (i = 6; i < nbody; ++i)
    {
        R = ((double) rand()) / ((double) RAND_MAX);
        G = ((double) rand()) / ((double) RAND_MAX) * (1.0 - R);
        B = 1.0 - R - G;
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

        color[i].x = R * scale;
        color[i].y = G * scale;
        color[i].z = B * scale;
    }

}

int nbodyInitDrawState(const NBodyState* st)
{
    _drawState = st;

    r = (FloatPos*) mwCallocA(st->nbody, sizeof(FloatPos));
    color = (FloatPos*) mwCallocA(st->nbody, sizeof(FloatPos));

    assignParticleColors(st->nbody);

    drawStateReady = TRUE;

    return 0;
}

int nbodyGLSetup(int* argc, char** argv)
{
#if 0

    const char* short_opts = "fpatldVh";
    const struct option long_opts[] =
    {
        {"fullscreen", no_argument, NULL, 'f'},
        {"paused", no_argument, NULL, 'p'},
        {"noaxes", no_argument, NULL, 'a'},
        {"time", no_argument, NULL, 't'},
        {"debug", no_argument, NULL, 'd'},
        {"version", no_argument, NULL, 'V'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };
#endif

    /* set a few initial variables */
    scene.nstar = 0;
    scene.xrot = 0.0f;
    scene.yrot = 0.0f;
    scene.z = -6.0f;
    scene.starsize = STARSIZE;
    scene.fullscreen = 0;
    scene.drawaxes = 1;
    scene.ntri = NTRI;
    scene.paused = 0;
    scene.step = 0;
    scene.mousemode = 1;
    scene.changed = 0;
    scene.dt = 300;

    if (0)
    {
        switch (0)
        {
            case 'f':
                scene.fullscreen = 1;
                break;
            case 'p':
                scene.paused = 1;
                break;
            case 'a':
                scene.drawaxes = 0;
                break;
            case 'd':
                debug = 1;
                break;
            default:
                break;
        }
    }

    glutInit(argc, argv);

    /* prepare rendering */
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutInitWindowPosition(0, 0);
    window = glutCreateWindow("Milkyway@Home N-body");
    glutDisplayFunc(drawGLScene);

    if (scene.fullscreen)
    {
        glutFullScreen();
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
    mwFreeA(r);
    r = NULL;

    mwFreeA(color);
    color = NULL;

    _drawState = NULL;
    drawStateReady = FALSE;  /* If we were to actually use this again, probably should use a lock */
}

/* Wait for it to be safe to draw and then run graphics loop */
void nbodyRunDisplayWhenReady()
{
    while (!drawStateReady);
    glutMainLoop();
}


#ifndef _WIN32

/* Launch main simulation as separate worker thread */

static void* runNBodySimulationInThreadWorker(void* nbf)
{
    /* Just return something non-null if error */
    return runNBodySimulation((const NBodyFlags*) nbf) ? nbf : NULL;
}

int runNBodySimulationInThread(NBodyFlags* nbf, NBodyThreadID* tid)
{
    if (pthread_create(tid, NULL, runNBodySimulationInThreadWorker, nbf))
    {
        perror("Creating worker thread");
        return 1;
    }

    return 0;
}

int nbodyWaitForSimThread(NBodyThreadID tid)
{
    void* ret;

    if (pthread_join(tid, &ret) != 0)
    {
        perror("Joining worker thread");
        return 1;
    }

    return ret != NULL;
}

#else

static DWORD __stdcall runNBodySimulationInThreadWorker(LPVOID nbf)
{
    return (DWORD) runNBodySimulation((const NBodyFlags*) nbf);
}

int runNBodySimulationInThread(NBodyFlags* nbf, NBodyThreadID* tid)
{
    *tid = CreateThread(NULL, 0, runNBodySimulationInThreadWorker, NULL, 0, nbf))
    if (!(*tid))
    {
        warn("Error creating display thread: %ld\n", GetLastError());
        return 1;
    }

    return 0;
}

int nbodyWaitForSimThread(NBodyThreadID tid)
{
    DWORD rc;

    if (WaitForSingleObject(tid, INFINITE) != WAIT_OBJECT_0)
    {
        warn("Error waiting for worker thread: %ld\n", GetLastError());
        return 1;
    }

    if (!GetExitCodeThread(tid, &rc))
    {
        warn("Failed to get worker thread exit status: %ld\n", GetLastError());
        return 1;
    }

    return (int) rc;
}

#endif /* _WIN32 */

