/*
  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
  Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _NBODY_TYPES_H_
#define _NBODY_TYPES_H_

/* Body and cell data structures are used to represent the tree.  During
 * tree construction, descendent pointers are stored in the subp arrays:
 *
 *          +-------------------------------------------------------------+
 * root --> | CELL: mass, pos, next, rcrit2, more, subp:[/,o,/,/,/,/,o,/] |
 *          +----------------------------------------------|---------|----+
 *                                                         |         |
 *     +---------------------------------------------------+         |
 *     |                                                             |
 *     |    +--------------------------------------+                 |
 *     +--> | BODY: mass, pos, next, vel, acc, phi |                 |
 *          +--------------------------------------+                 |
 *                                                                   |
 *     +-------------------------------------------------------------+
 *     |
 *     |    +-------------------------------------------------------------+
 *     +--> | CELL: mass, pos, next, rcrit2, more, subp:[o,/,/,o,/,/,o,/] |
 *          +--------------------------------------------|-----|-----|----+
 *                                                      etc   etc   etc
 *
 * After the tree is complete, it is threaded to permit linear force
 * calculation, using the next and more pointers.  The storage used for
 * the subp arrays may be reused to store quadrupole moments.
 *
 *          +-----------------------------------------------+
 * root --> | CELL: mass, pos, next:/, rcrit2, more:o, quad |
 *          +---------------------------------------|-------+
 *                                                  |
 *     +--------------------------------------------+
 *     |
 *     |    +----------------------------------------+
 *     +--> | BODY: mass, pos, next:o, vel, acc, phi |
 *          +-----------------------|----------------+
 *                                  |
 *     +----------------------------+
 *     |
 *     |    +-----------------------------------------------+
 *     +--> | CELL: mass, pos, next:/, rcrit2, more:o, quad |
 *          +---------------------------------------|-------+
 *                                                 etc
 */

#include "nbody_config.h"
#include "milkyway_math.h"
#include "milkyway_extra.h"
#include "milkyway_util.h"
#include "nbody_graphics.h"
#include "nbody_potential_types.h"

#include <lua.h>
#include <time.h>


/* There are bodies and cells. Cells are 0, bodies are nonzero. Bodies
   will be set to 1 or -1 if the body is to be ignored in the final
   likelihood calculation (i.e. a dark matter body)
 */
#define CELL(x) (0)
#define BODY(x) ((x) ? -1 : 1)

#define isBody(x) (((NBodyNode*) (x))->type != 0)
#define isCell(x) (((NBodyNode*) (x))->type == 0)
#define ignoreBody(x) (((NBodyNode*) (x))->type < 0)
#define bodyTypeIsIgnore(x) ((x) < 0)

typedef short body_t;

/* node: data common to BODY and CELL structures. */
typedef struct NBODY_ALIGN _NBodyNode
{
    body_t type;              /* code for node type */
    real mass;                /* total mass of node */
    mwvector pos;             /* position of node */
    struct _NBodyNode* next;  /* link to next force-calc */
} NBodyNode;

#define EMPTY_NODE { 0, 0.0, ZERO_VECTOR, NULL }

#define Type(x) (((NBodyNode*) (x))->type)
#define Mass(x) (((NBodyNode*) (x))->mass)
#define Pos(x)  (((NBodyNode*) (x))->pos)
#define Next(x) (((NBodyNode*) (x))->next)

/* BODY: data structure used to represent particles. */

typedef struct NBODY_ALIGN
{
    NBodyNode bodynode;         /* data common to all nodes */
    mwvector vel;               /* velocity of body */
} Body;

#define EMPTY_BODY { EMPTY_NODE, ZERO_VECTOR }

#define BODY_TYPE "Body"

#define Vel(x)  (((Body*) (x))->vel)

/* CELL: structure used to represent internal nodes of tree. */

#define NSUB (1 << NDIM)        /* subcells per cell */

typedef struct NBODY_ALIGN
{
    NBodyNode cellnode;         /* data common to all nodes */
    real rcrit2;                /* critical c-of-m radius^2 */
    NBodyNode* more;            /* link to first descendent */
    union                       /* shared storage for... */
    {
        NBodyNode* subp[NSUB];  /* descendents of cell */
        mwmatrix quad;          /* quad. moment of cell */
    } stuff;
} NBodyCell;

/* use alternate criteria */
typedef enum
{
    InvalidCriterion = InvalidEnum,
    NewCriterion,  /* FIXME: What is this exactly? Rename it. */
    SW93,
    BH86,
    Exact
} criterion_t;


typedef enum
{
    EXTERNAL_POTENTIAL_DEFAULT,
    EXTERNAL_POTENTIAL_NONE,
    EXTERNAL_POTENTIAL_CUSTOM_LUA
} ExternalPotentialType;


#define Rcrit2(x) (((NBodyCell*) (x))->rcrit2)
#define More(x)   (((NBodyCell*) (x))->more)
#define Subp(x)   (((NBodyCell*) (x))->stuff.subp)
#define Quad(x)   (((NBodyCell*) (x))->stuff.quad)

/* Variables used in tree construction. */

typedef struct NBODY_ALIGN
{
    NBodyCell* root;         /* pointer to root cell */
    real rsize;              /* side-length of root cell */

    unsigned int cellused;   /* count of cells in tree */
    unsigned int maxlevel;   /* count of levels in tree */
    int structureError;
} NBodyTree;


#ifndef _WIN32

typedef struct NBODY_ALIGN
{
    int fd;            /* File descriptor for checkpoint file */
    char* mptr;        /* mmap'd pointer for checkpoint file */
    size_t cpFileSize; /* For checking how big the file should be for expected bodies */
} CheckpointHandle;

#define EMPTY_CHECKPOINT_HANDLE { -1, NULL, 0 }

#else

typedef struct NBODY_ALIGN
{
    HANDLE file;
    HANDLE mapFile;
    char* mptr;
    DWORD cpFileSize;
} CheckpointHandle;

#define EMPTY_CHECKPOINT_HANDLE { INVALID_HANDLE_VALUE, INVALID_HANDLE_VALUE, NULL, 0 }

#endif /* _WIN32 */


/* Mutable state used during an evaluation */
typedef struct NBODY_ALIGN
{
    NBodyTree tree;
    NBodyNode* freecell;   /* list of free cells */
    time_t lastCheckpoint;
    real tnow;
    int nbody;
    Body* bodytab;      /* points to array of bodies */
    mwvector* acctab;   /* Corresponding accelerations of bodies */
    int treeIncest;     /* Tree incest has occured */

    FILE* outFile;            /* file for snapshot output */
    char* checkpointResolved;

    mwvector* orbitTrace;  /* Trail of center of masses for display purposes */
    scene_t* scene;
    int shmId; /* shmid, key when using shmem */
} NBodyState;

#define NBODYSTATE_TYPE "NBodyState"

#define EMPTY_NBODYSTATE { EMPTY_TREE, NULL, 0, 0.0, 0, NULL, NULL, FALSE, NULL, NULL, NULL, NULL, -1 }


typedef struct
{
    real phi;
    real theta;
    real psi;
    real startRaw;
    real endRaw;
    real binSize;
    real center;
} HistogramParams;

#define EMPTY_HISTOGRAM_PARAMS { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#define HISTOGRAM_PARAMS_TYPE "HistogramParams"


typedef struct
{
    int useBin;
    real lambda;
    real err;
    real count;
} HistData;


/* The context tracks settings of the simulation.  It should be set
   once at the beginning of a simulation based on settings, and then
   stays constant for the actual simulation.
 */
typedef struct NBODY_ALIGN
{
    Potential pot;
    ExternalPotentialType potentialType;

    real timestep;
    real timeEvolve;

    real theta;               /* accuracy parameter: 0.0 */
    real eps2;                /* (potential softening parameter)^2 */
    real treeRSize;

    real sunGCDist;
    criterion_t criterion;

    mwbool useQuad;           /* use quadrupole corrections */
    mwbool allowIncest;
    mwbool quietErrors;

    time_t checkpointT;       /* Period to checkpoint when not using BOINC */
    unsigned int freqOut;
    HistogramParams histogramParams;
} NBodyCtx;

#define NBODYCTX_TYPE "NBodyCtx"

/* Negative codes can be nonfatal but useful return statuses.
   Positive can be different hard failures.
 */
typedef enum
{
    NBODY_TREE_INCEST_NONFATAL = -(1 << 2), /* Negative of NBODY_TREE_INCEST */
    NBODY_SUCCESS              = 0 << 0,
    NBODY_ERROR                = 1 << 0,
    NBODY_TREE_STRUCTURE_ERROR = 1 << 1,
    NBODY_TREE_INCEST_FATAL    = 1 << 2,
    NBODY_IO_ERROR             = 1 << 3,
    NBODY_CHECKPOINT_ERROR     = 1 << 4
} NBodyStatus;

#define nbodyStatusIsFatal(x) ((x) > 0)


/* Note: 'type' should first field for all types. */
#define SET_TYPE(x, y) (((Disk*)x)->type = y)
#define NBODY_TYPEOF(x) (((Disk*)x)->type)


#define EMPTY_TREE { NULL, 0.0, 0, 0, FALSE }
#define EMPTY_NBODYCTX { EMPTY_POTENTIAL, EXTERNAL_POTENTIAL_DEFAULT, 0.0, \
                         0.0, 0.0, 0.0, 0.0,                               \
                         0.0, InvalidCriterion,                            \
                         FALSE, FALSE, FALSE,                              \
                         0, 0, EMPTY_HISTOGRAM_PARAMS }

int destroyNBodyState(NBodyState* st);
int detachSharedScene(NBodyState* st);
void setInitialNBodyState(NBodyState* st, const NBodyCtx* ctx, Body* bodies, int nbody);
void cloneNBodyState(NBodyState* st, const NBodyState* oldSt);
int equalNBodyState(const NBodyState* st1, const NBodyState* st2);

void sortBodies(Body* bodies, int nbody);

int equalSpherical(const Spherical* s1, const Spherical* s2);
int equalHalo(const Halo* h1, const Halo* h2);
int equalDisk(const Disk* d1, const Disk* d2);
int equalPotential(const Potential* p1, const Potential* p2);

int equalNBodyCtx(const NBodyCtx* ctx1, const NBodyCtx* ctx2);

int equalHistogramParams(const HistogramParams* hp1, const HistogramParams* hp2);

#endif /* _NBODY_TYPES_H_ */

