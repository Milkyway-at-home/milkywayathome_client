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

#include <lua.h>
#include <time.h>

#ifndef __OPENCL_VERSION__   /* Not compiling CL kernel */
  #if NBODY_OPENCL
    #include <OpenCL/cl.h>
    #include <OpenCL/cl_platform.h>
    #include "build_cl.h"
  #endif /* NBODY_OPENCL */

  #include <stdio.h>
  #ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
	#define VC_EXTRALEAN
    #include <windows.h>
  #endif /* _WIN32 */

#else
  /* FIXME: Remove IO from context?
  These aren't allowed in the kernels, so make these go away with same size type.
  CHECKME: type of HANDLE on windows = ? */
  #define FILE void
  #define HANDLE void*
#endif

#ifndef _MSC_VER
  #define NBODY_ALIGN __attribute__((packed))
#else
  #define NBODY_ALIGN
#endif /* _MSC_VER */

/* There are bodies and cells. Cells are 0, bodies are nonzero. Bodies
   will be set to 1 or -1 if the body is to be ignored in the final
   likelihood calculation (i.e. a dark matter body)
 */
#define CELL(x) (0)
#define BODY(x) ((x) ? -1 : 1)

#define isBody(x) (((Node*) (x))->type != 0)
#define isCell(x) (((Node*) (x))->type == 0)
#define ignoreBody(x) (((Node*) (x))->type < 0)
#define bodyTypeIsIgnore(x) ((x) < 0)

typedef short body_t;

/* node: data common to BODY and CELL structures. */
typedef struct NBODY_ALIGN _Node
{
    body_t type;            /* code for node type */
    real mass;              /* total mass of node */
    mwvector pos;           /* position of node */
    struct _Node* next;     /* link to next force-calc */
} Node;

#define EMPTY_NODE { 0, NAN, EMPTY_MWVECTOR, NULL }

#define Type(x) (((Node*) (x))->type)
#define Mass(x) (((Node*) (x))->mass)
#define Pos(x)  (((Node*) (x))->pos)
#define Next(x) (((Node*) (x))->next)

/* BODY: data structure used to represent particles. */

typedef struct NBODY_ALIGN
{
    Node bodynode;              /* data common to all nodes */
    mwvector vel;               /* velocity of body */
} Body;

#define EMPTY_BODY { EMPTY_NODE, EMPTY_MWVECTOR }

#define BODY_TYPE "Body"

#define Vel(x)  (((Body*) (x))->vel)

/* CELL: structure used to represent internal nodes of tree. */

#define NSUB (1 << NDIM)        /* subcells per cell */

typedef struct NBODY_ALIGN
{
    Node cellnode;              /* data common to all nodes */
    real rcrit2;                /* critical c-of-m radius^2 */
    Node* more;                 /* link to first descendent */
    union                       /* shared storage for... */
    {
        Node* subp[NSUB];       /* descendents of cell */
        mwmatrix quad;         /* quad. moment of cell */
    } stuff;
} Cell;

#define InvalidEnum (-1)
typedef int generic_enum_t;  /* A general enum type. */

/* use alternate criteria */
typedef enum
{
    InvalidCriterion = InvalidEnum,
    NewCriterion,  /* FIXME: What is this exactly? Rename it. */
    Exact,
    BH86,
    SW93
} criterion_t;

#define _SPHERICAL 0

typedef enum
{
    InvalidSpherical   = InvalidEnum,
    SphericalPotential = _SPHERICAL
} spherical_t;

/* Spherical potential */
typedef struct NBODY_ALIGN
{
    spherical_t type;
    real mass;
    real scale;
} Spherical;

#define SPHERICAL_TYPE "Spherical"


/* Can't get the enum value in preprocessor, so do this */
#define _MN_DISK 0
#define _EXP_DISK 1


/* Supported disk models */
typedef enum
{
    InvalidDisk       = InvalidEnum,
    MiyamotoNagaiDisk = _MN_DISK,
    ExponentialDisk   = _EXP_DISK
} disk_t;

typedef struct NBODY_ALIGN
{
    disk_t type;
    real mass;         /* disk mass */
    real scaleLength;  /* "a" for M-N, "b" for exp disk */
    real scaleHeight;  /* unused for exponential disk. "b" for Miyamoto-Nagai disk */
} Disk;

#define DISK_TYPE "Disk"

/* Supported halo models */

/* Can't get the enum value in preprocessor, so do this */
#define _LOG_HALO 0
#define _NFW_HALO 1
#define _TRIAXIAL_HALO 2
typedef enum
{
    InvalidHalo     = InvalidEnum,
    LogarithmicHalo = _LOG_HALO,
    NFWHalo         = _NFW_HALO,
    TriaxialHalo    = _TRIAXIAL_HALO
} halo_t;

typedef struct NBODY_ALIGN
{
    halo_t type;
    real vhalo;         /* common to all 3 halos */
    real scaleLength;   /* common to all 3 halos */
    real flattenZ;      /* used by logarithmic and triaxial */
    real flattenY;      /* used by triaxial */
    real flattenX;      /* used by triaxial */
    real triaxAngle;    /* used by triaxial */

    real c1;           /* Constants calculated for triaxial from other params */
    real c2;        /* TODO: Lots more stuff could be cached, but should be done less stupidly */
    real c3;
} Halo;

#define HALO_TYPE "Halo"

typedef struct NBODY_ALIGN
{
    Spherical sphere[1];  /* 1 for now, flexibility can be added later */
    Disk disk;
    Halo halo;
    void* rings;         /* reserved for future use */
} Potential;

#define POTENTIAL_TYPE "Potential"


#define Rcrit2(x) (((Cell*) (x))->rcrit2)
#define More(x)   (((Cell*) (x))->more)
#define Subp(x)   (((Cell*) (x))->stuff.subp)
#define Quad(x)   (((Cell*) (x))->stuff.quad)

/* Variables used in tree construction. */

typedef struct NBODY_ALIGN
{
    Cell* root;     /* pointer to root cell */
    real rsize;     /* side-length of root cell */

    unsigned int cellused;   /* count of cells in tree */
    unsigned int maxlevel;   /* count of levels in tree */
    mwbool firstcall;
} Tree;

typedef struct NBODY_ALIGN
{
    mwvector position;
    mwvector velocity;
} InitialConditions;

#define EMPTY_INITIAL_CONDITIONS { EMPTY_MWVECTOR, EMPTY_MWVECTOR }

#define INITIAL_CONDITIONS_TYPE "InitialConditions"


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

#if NBODY_OPENCL && !defined(__OPENCL_VERSION__)

typedef struct
{
    cl_mem acc;
    cl_mem bodies;
    cl_mem nbctx;
    cl_mem root;
} NBodyCLMem;

#define EMPTY_NBODY_CL_MEM { NULL, NULL, NULL, NULL }

#endif /* NBODY_OPENCL && !defined(__OPENCL_VERSION__)) */


#ifndef __OPENCL_VERSION__

/* Mutable state used during an evaluation */
typedef struct NBODY_ALIGN
{
    Tree tree;
    Node* freecell;   /* list of free cells */
    unsigned int outputTime;
    time_t lastCheckpoint;
    real tnow;
    unsigned int nbody;
    Body* bodytab;      /* points to array of bodies */
    mwvector* acctab;   /* Corresponding accelerations of bodies */

  #if NBODY_OPENCL
    CLInfo ci;
    NBodyCLMem cm;
  #endif /* NBODY_OPENCL */
} NBodyState;

#define NBODYSTATE_TYPE "NBodyState"

#if NBODY_OPENCL
  #define EMPTY_NBODYSTATE { EMPTY_TREE, NULL, 0, NAN, NULL, NULL, EMPTY_CL_INFO, EMPTY_NBODY_CL_MEM }
#else
  #define EMPTY_NBODYSTATE { EMPTY_TREE, NULL, 0, 0, NAN, 0, NULL, NULL }
#endif /* NBODY_OPENCL */


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

#define EMPTY_HISTOGRAM_PARAMS { NAN, NAN, NAN, NAN, NAN, NAN, NAN }
#define HISTOGRAM_PARAMS_TYPE "HistogramParams"


typedef struct
{
    int useBin;
    real lambda;
    real err;
    real count;
} HistData;


#endif /* __OPENCL_VERSION__ */

/* The context tracks settings of the simulation.  It should be set
   once at the beginning of a simulation based on settings, and then
   stays constant for the actual simulation.
 */
typedef struct NBODY_ALIGN
{
    Potential pot;

    real timestep;
    real timeEvolve;

    real theta;               /* accuracy parameter: 0.0 */
    real eps2;                /* (potential softening parameter)^2 */
    real treeRSize;

    real sunGCDist;
    criterion_t criterion;

    mwbool useQuad;           /* use quadrupole corrections */
    mwbool allowIncest;

    time_t checkpointT;       /* Period to checkpoint when not using BOINC */
    unsigned int freqOut;
    FILE* outfile;            /* file for snapshot output */
    HistogramParams histogramParams;
} NBodyCtx;

#define NBODYCTX_TYPE "NBodyCtx"


/* Note: 'type' should first field for all types. */
#define SET_TYPE(x, y) (((Disk*)x)->type = y)
#define NBODY_TYPEOF(x) (((Disk*)x)->type)


/* Useful initializers */
#define EMPTY_SPHERICAL { InvalidSpherical, NAN, NAN }
#define EMPTY_DISK { InvalidDisk, NAN, NAN, NAN }
#define EMPTY_HALO { InvalidHalo, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN }
#define EMPTY_POTENTIAL { {EMPTY_SPHERICAL}, EMPTY_DISK, EMPTY_HALO, NULL }

#define EMPTY_TREE { NULL, NAN, 0, 0, TRUE }
#define EMPTY_NBODYCTX { EMPTY_POTENTIAL, NAN, NAN,                       \
                         NAN, NAN, NAN,                                   \
                         NAN, InvalidCriterion,                           \
                         FALSE, FALSE,                                    \
                         0, 0, NULL, EMPTY_HISTOGRAM_PARAMS }


int destroyNBodyCtx(NBodyCtx* ctx);
void destroyNBodyState(NBodyState* st);
void setInitialNBodyState(NBodyState* st, const NBodyCtx* ctx, Body* bodies, unsigned int nbody);
void cloneNBodyState(NBodyState* st, const NBodyState* oldSt, const unsigned int nbody);

void sortBodies(Body* bodies, unsigned int nbody);

#ifndef __OPENCL_VERSION__  /* No function pointers allowed in kernels */
/* Acceleration functions for a given potential */
typedef mwvector (*SphericalAccel) (const Spherical*, const mwvector);
typedef mwvector (*HaloAccel) (const Halo*, const mwvector);
typedef mwvector (*DiskAccel) (const Disk*, const mwvector);

/* Generic potential function */
typedef mwvector (*AccelFunc) (const void*, const mwvector);

#endif /* __OPENCL_VERSION__ */

#endif /* _NBODY_TYPES_H_ */

