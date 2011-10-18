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


#if NBODY_OPENCL
  #include "milkyway_cl.h"
#endif



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
    mwvector pos;             /* position of node */
    struct _NBodyNode* next;  /* link to next force-calc */
    real mass;                /* total mass of node */
    body_t type;              /* code for node type */
} NBodyNode;

#define EMPTY_NODE { ZERO_VECTOR, NULL, 0.0, 0  }

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
    real xx, xy, xz;
    real yy, yz;
    real zz;
} NBodyQuadMatrix;



typedef struct NBODY_ALIGN
{
    NBodyNode cellnode;         /* data common to all nodes */
    real rcrit2;                /* critical c-of-m radius^2 */
    NBodyNode* more;            /* link to first descendent */
    union                       /* shared storage for... */
    {
        NBodyNode* subp[NSUB];  /* descendents of cell */
        NBodyQuadMatrix quad;   /* quad. moment of cell. Unique symmetric matrix components */
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


#if NBODY_OPENCL

typedef struct
{
    cl_mem pos[3];
    cl_mem vel[3];
    cl_mem acc[3];
    cl_mem max[3];
    cl_mem min[3];
    cl_mem masses;

    cl_mem start; /* TODO: We can reuse other buffers with this later to save memory */
    cl_mem count;
    cl_mem child;
    cl_mem sort;

    cl_mem critRadii; /* Used by the alternative cell opening criterion.
                         Unnecessary for BH86. */

    struct
    {
        cl_mem xx, xy, xz;
        cl_mem yy, yz;
        cl_mem zz;
    } quad;

    cl_mem treeStatus;
} NBodyBuffers;


/* 6 used by tree + 1 with quad, 2 used by exact. 1 shared. */
#define NKERNELS 8


typedef struct
{
    cl_kernel boundingBox;
    cl_kernel buildTree;
    cl_kernel summarization;
    cl_kernel sort;
    cl_kernel quadMoments;
    cl_kernel forceCalculation;
    cl_kernel integration;

    /* Used by exact one only */
    cl_kernel forceCalculation_Exact;
} NBodyKernels;

#endif /* NBODY_OPENCL */


typedef struct
{
    size_t factors[7];
    size_t threads[7];
    double timings[7];        /* In a single iteration */
    double chunkTimings[7];   /* Average time per chunk */
    double kernelTimings[7];  /* Running totals */

    size_t global[7];
    size_t local[7];
} NBodyWorkSizes;



/* Mutable state used during an evaluation */
typedef struct NBODY_ALIGN
{
    NBodyTree tree;
    NBodyNode* freecell;      /* list of free cells */
    char* checkpointResolved;
    Body* bodytab;            /* points to array of bodies */
    mwvector* acctab;         /* Corresponding accelerations of bodies */
    mwvector* orbitTrace;     /* Trail of center of masses for display purposes */
    scene_t* scene;

    lua_State** potEvalStates;  /* If using a Lua closure as a potential, the evaluation states.
                                  We need one per thread in the general case. */
    int* potEvalClosures;       /* Lua closure for each state */

    time_t lastCheckpoint;

    unsigned int step;
    int nbody;
    int treeIncest;          /* Tree incest has occured */
    int potentialEvalError;  /* Error occured in calling custom Lua potential */

    int shmId;               /* shmid, key when using shmem */

    mwbool ignoreResponsive;
    mwbool usesExact;
    mwbool usesQuad;
    mwbool dirty;      /* Whether the view of the bodies is consistent with the view in the CL buffers */
    mwbool usesCL;
    mwbool reportProgress;

  #if NBODY_OPENCL
    CLInfo* ci;
    NBodyKernels* kernels;
    NBodyBuffers* nbb;
  #else
    void* kernels;
    void* ci;
    void* nbb;
  #endif /* NBODY_OPENCL */
    NBodyWorkSizes* workSizes;
} NBodyState;

#define NBODYSTATE_TYPE "NBodyState"

#define EMPTY_NBODYSTATE { EMPTY_TREE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, -1, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, NULL, NULL, NULL, NULL }


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
    real eps2;                /* (potential softening parameter)^2 */
    real theta;               /* accuracy parameter: 0.0 */
    real timestep;
    real timeEvolve;
    real treeRSize;
    real sunGCDist;

    criterion_t criterion;
    ExternalPotentialType potentialType;

    mwbool useQuad;           /* use quadrupole corrections */
    mwbool allowIncest;
    mwbool quietErrors;

    time_t checkpointT;       /* Period to checkpoint when not using BOINC */
    unsigned int freqOut;
    unsigned int nStep;

    Potential pot;
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
    NBODY_CHECKPOINT_ERROR     = 1 << 4,
    NBODY_CL_ERROR             = 1 << 5,
    NBODY_CAPABILITY_ERROR     = 1 << 6,
    NBODY_CONSISTENCY_ERROR    = 1 << 7,
    NBODY_UNIMPLEMENTED        = 1 << 8,
    NBODY_UNSUPPORTED          = 1 << 9,
    NBODY_USER_ERROR           = 1 << 10,
    NBODY_PARAM_FILE_ERROR     = 1 << 11,
    NBODY_LUA_POTENTIAL_ERROR  = 1 << 12
} NBodyStatus;

#define nbStatusIsFatal(x) ((x) > 0)
#define nbStatusIsOK(x) ((x) <= 0)
#define nbStatusIsWarning(x) ((x) < 0)

/* Reserve positive numbers for reporting depth > MAXDEPTH. Should match in kernel  */
typedef enum
{
    NBODY_KERNEL_OK                   = 0,
    NBODY_KERNEL_CELL_LEQ_NBODY       = -1,
    NBODY_KERNEL_TREE_INCEST          = -2,
    NBODY_KERNEL_TREE_STRUCTURE_ERROR = -3,
    NBODY_KERNEL_ERROR_OTHER          = -4
} NBodyKernelError;


/* Note: 'type' should first field for all types. */
#define SET_TYPE(x, y) (((Disk*)x)->type = y)
#define NBODY_TYPEOF(x) (((Disk*)x)->type)


#define EMPTY_TREE { NULL, 0.0, 0, 0, FALSE }
#define EMPTY_NBODYCTX { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                  \
                         InvalidCriterion, EXTERNAL_POTENTIAL_DEFAULT,  \
                         FALSE, FALSE, FALSE,                           \
                         0, 0, 0,                                       \
                         EMPTY_POTENTIAL, EMPTY_HISTOGRAM_PARAMS }



#if NBODY_OPENCL
NBodyStatus initCLNBodyState(NBodyState* st, const NBodyCtx* ctx, const CLRequest* clr);
#endif

int destroyNBodyState(NBodyState* st);
int nbDetachSharedScene(NBodyState* st);
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

