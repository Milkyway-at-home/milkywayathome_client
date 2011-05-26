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

#define USLEEPDT 100.0

#define USE_GL_POINTS 0



static FloatPos* color = NULL;

static int window, xlast = 0, ylast = 0, debug = 0, width = 1024, height = 768;

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

    scene->changed = TRUE;
}

/* The main drawing function */
static void drawGLScene()
{
    unsigned int i = 0;
    int nbody = scene->nbody;
    FloatPos* r;

    if (!scene)
        return;

    usleep(scene->usleepdt);

    /* draw scene if necessary */
    if (scene->changed)
    {
        scene->changed = FALSE;

        /* erase scene */
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        glLoadIdentity();

        /* get ready to render 3D objects */
        glTranslatef(0.0f, 0.0f, scene->z + 0.1);

        /* rotate view--I know these two rotation matrices don't commute */
        glRotatef(scene->xrot, 1.0f, 0.0f, 0.0f);
        glRotatef(scene->yrot, 0.0f, 1.0f, 0.0f);

        /* draw axes */
        if (scene->drawaxes)
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

        r = &scene->r[0];
        /* draw stars */
    #if !USE_GL_POINTS
        /* FIXME: Use modern OpenGL stuff */
        for (i = 0; i < nbody; ++i)
        {
            glLoadIdentity();
            glTranslatef(0.0f, 0.0f, scene->z);
            glRotatef(scene->xrot, 1.0f, 0.0f, 0.0f);
            glRotatef(scene->yrot, 0.0f, 1.0f, 0.0f);

            glTranslatef(r[i].x / SCALE, r[i].y / SCALE, r[i].z / SCALE);
            glColor3f(color[i].x, color[i].y, color[i].z);
            /* glutSolidSphere(scene->starsize, scene->ntri, scene->ntri); */
            gluSphere(quadratic, scene->starsize, scene->ntri, scene->ntri);
            /* glTranslatef(-r[i].x/SCALE, -r[i].y/SCALE, -r[i].z/SCALE); */
        }
    #else
        glPointSize(scene->starsize);  /* Is this actually working? */
        glBegin(GL_POINTS);
        glColor3f(1.0f, 1.0f, 1.0f);
        for (i = 0; i < nbody; ++i)
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
            scene->step = 1;
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
    /* avoid thrashing this procedure */
    usleep(100);

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

    dx = ((double) x) - ((double) xlast);
    dy = ((double) y) - ((double) ylast);

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

static void assignParticleColors(unsigned int nbody)
{
    int i;
    double R, G, B, scale;

    /* assign random particle colors */
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

int nbodyInitDrawState()
{
    color = (FloatPos*) mwCallocA(scene->nbody, sizeof(FloatPos));
    assignParticleColors(scene->nbody);

    return 0;
}

#if !BOINC_APPLICATION

static int key = -1;

static int setShmemKey()
{
    /* TODO: pid of simulation, etc. */
    key = DEFAULT_SHMEM_KEY;

    return 0;
}

static int connectSharedScene()
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

static int connectSharedScene()
{
    return 0;
}

#endif /* !BOINC_APPLICATION */

static void sceneInit()
{
    scene->xrot = 0.0f;
    scene->yrot = 0.0f;
    scene->z = -6.0f;
    scene->starsize = STARSIZE;
    scene->fullscreen = FALSE;
    scene->drawaxes = TRUE;
    scene->ntri = NTRI;
    scene->paused = FALSE;
    scene->step = 0;
    scene->mousemode = 1;
    scene->changed = FALSE;
    scene->dt = 300;
    scene->usleepdt = USLEEPDT;
}

int nbodyGLSetup(int* argc, char** argv)
{
    if (connectSharedScene())
        return 1;

    sceneInit();
    glutInit(argc, argv);

    /* prepare rendering */
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutInitWindowPosition(0, 0);
    window = glutCreateWindow("Milkyway@Home N-body");
    glutDisplayFunc(drawGLScene);

    if (scene->fullscreen)
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
    mwFreeA(color);
    color = NULL;
}

