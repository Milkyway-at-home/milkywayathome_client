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

#if NBODY_OPENCL && !defined(bool)
    typedef int bool;
#elif !defined(bool)
    typedef short int bool;
#endif /* NBODY_OPENCL && !defined(bool) */

#ifndef TRUE
  #define TRUE  1
  #define FALSE 0
#endif

#ifndef _MSC_VER
  #define NBODY_ALIGN __attribute__((packed))
#else
  #define NBODY_ALIGN
#endif /* _MSC_VER */

/*
typedef enum
{
    BODY,
    CELL
} body_t; */
#define BODY 01
#define CELL 02

typedef short body_t;

/* node: data common to BODY and CELL structures. */
typedef struct NBODY_ALIGN _node
{
    body_t type;            /* code for node type */
    real mass;              /* total mass of node */
    mwvector pos;           /* position of node */
    struct _node* next;     /* link to next force-calc */
} node, *nodeptr;

#define Type(x) (((nodeptr) (x))->type)
#define Mass(x) (((nodeptr) (x))->mass)
#define Pos(x)  (((nodeptr) (x))->pos)
#define Next(x) (((nodeptr) (x))->next)

/* BODY: data structure used to represent particles. */

typedef struct NBODY_ALIGN
{
    node bodynode;              /* data common to all nodes */
    mwvector vel;               /* velocity of body */
 } body, *bodyptr;

#define Body    body

#define Vel(x)  (((bodyptr) (x))->vel)

/* CELL: structure used to represent internal nodes of tree. */

#define NSUB (1 << NDIM)        /* subcells per cell */

typedef struct NBODY_ALIGN
{
    node cellnode;              /* data common to all nodes */
    real rcrit2;                /* critical c-of-m radius^2 */
    nodeptr more;               /* link to first descendent */
    union                       /* shared storage for... */
    {
        nodeptr subp[NSUB];     /* descendents of cell */
        mwmatrix quad;         /* quad. moment of cell */
    } stuff;
} cell, *cellptr;

/* use alternate criteria */
typedef enum
{
    NEWCRITERION,  /* FIXME: What is this exactly? Rename it. */
    EXACT,
    BH86,
    SW93
} criterion_t;

#define _SPHERICAL 0

typedef enum
{
    SphericalPotential = _SPHERICAL
} spherical_t;

/* Spherical potential */
typedef struct NBODY_ALIGN
{
    spherical_t type;
    real mass;
    real scale;
} Spherical;


/* Can't get the enum value in preprocessor, so do this */
#define _MN_DISK 0
#define _EXP_DISK 1


/* Supported disk models */
typedef enum
{
    MiyamotoNagaiDisk = _MN_DISK,
    ExponentialDisk   = _EXP_DISK
} disk_t;

typedef struct NBODY_ALIGN
{
    disk_t type;
    real mass;          /* disk mass */
    real scale_length;  /* "a" for M-N, "b" for exp disk */
    real scale_height;  /* unused for exponential disk. "b" for Miyamoto-Nagai disk */
} Disk;

/* Supported halo models */

/* Can't get the enum value in preprocessor, so do this */
#define _LOG_HALO 0
#define _NFW_HALO 1
#define _TRIAXIAL_HALO 2
typedef enum
{
    LogarithmicHalo = _LOG_HALO,
    NFWHalo         = _NFW_HALO,
    TriaxialHalo    = _TRIAXIAL_HALO
} halo_t;

typedef struct NBODY_ALIGN
{
    halo_t type;
    real vhalo;         /* common to all 3 halos */
    real scale_length;  /* common to all 3 halos */
    real flattenZ;      /* used by logarithmic and triaxial */
    real flattenY;      /* used by triaxial */
    real flattenX;      /* used by triaxial */
    real triaxAngle;    /* used by triaxial */

    real c1;           /* Constants calculated for triaxial from other params */
    real c2;        /* TODO: Lots more stuff could be cached, but should be done less stupidly */
    real c3;
} Halo;

typedef struct NBODY_ALIGN
{
    Spherical sphere[1];  /* 1 for now, flexibility can be added later */
    Disk disk;
    Halo halo;
    void* rings;         /* reserved for future use */
} Potential;


#define Rcrit2(x) (((cellptr) (x))->rcrit2)
#define More(x)   (((cellptr) (x))->more)
#define Subp(x)   (((cellptr) (x))->stuff.subp)
#define Quad(x)   (((cellptr) (x))->stuff.quad)

/* Variables used in tree construction. */

typedef struct NBODY_ALIGN
{
    cellptr root;   /* pointer to root cell */
    real rsize;     /* side-length of root cell */

    unsigned int cellused;   /* count of cells in tree */
    unsigned int maxlevel;   /* count of levels in tree */
} Tree;

typedef struct NBODY_ALIGN
{
    mwvector position;     /* (x, y, z) if cartesian / useGalC, otherwise (l, b, r) */
    mwvector velocity;
    bool useGalC;
    bool useRadians;
    bool reverseOrbit;     /* Do a reverse orbit. Otherwise, take initial conditions as they are */
} InitialConditions;


typedef struct NBODY_ALIGN
{
    int useFitParams;
    real modelMass;
    real modelRadius;
    real reverseOrbitTime;
    real simulationTime;
} FitParams;

#define EMPTY_FIT_PARAMS { FALSE, NAN, NAN, NAN, NAN }


typedef enum
{
    DwarfModelPlummer,
    DwarfModelKing,
    DwarfModelDehnen
} dwarf_model_t;

typedef struct NBODY_ALIGN
{
    dwarf_model_t type;
    int nbody;
    real time_orbit;

    /* calculated depending on model */
    real timestep;
    real orbit_timestep;

    /* model parameters */
    real mass;
    real scale_radius;

    bool ignoreFinal;
    InitialConditions initialConditions;
} DwarfModel;

#ifndef _WIN32

typedef struct NBODY_ALIGN
{
    int fd;            /* File descriptor for checkpoint file */
    char* mptr;        /* mmap'd pointer for checkpoint file */
} CheckpointHandle;

#define EMPTY_CHECKPOINT_HANDLE { -1, NULL }

#else

typedef struct NBODY_ALIGN
{
    HANDLE file;
    HANDLE mapFile;
    char* mptr;
} CheckpointHandle;

#define EMPTY_CHECKPOINT_HANDLE { INVALID_HANDLE_VALUE, INVALID_HANDLE_VALUE, NULL }

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
    nodeptr freecell;   /* list of free cells */
    real tout;
    real tnow;
    bodyptr bodytab;    /* points to array of bodies */
    mwvector* acctab;   /* Corresponding accelerations of bodies */

  #if NBODY_OPENCL
    CLInfo ci;
    NBodyCLMem cm;
  #endif /* NBODY_OPENCL */
} NBodyState;

#if NBODY_OPENCL
  #define EMPTY_STATE { EMPTY_TREE, NULL, NAN, NAN, NULL, NULL, EMPTY_CL_INFO, EMPTY_NBODY_CL_MEM }
#else
  #define EMPTY_STATE { EMPTY_TREE, NULL, NAN, NAN, NULL, NULL }
#endif /* NBODY_OPENCL */

#endif /* __OPENCL_VERSION__ */

/* The context tracks settings of the simulation.  It should be set
   once at the beginning of a simulation based on settings, and then
   stays constant for the actual simulation.
 */
typedef struct NBODY_ALIGN
{
    Potential pot;

    DwarfModel* models;       /* dwarf models */
    unsigned int modelNum;    /* Number of models */
    unsigned int nbody;       /* Total number of bodies in all models */

    real timestep;
    real time_evolve;
    real orbit_timestep;
    real time_orbit;

    char* headline;           /* message describing calculation */
    const char* outfilename;  /* output */
    const char* histogram;
    const char* histout;
    FILE* outfile;            /* file for snapshot output */

    real freqout;
    real theta;               /* accuracy parameter: 0.0 */
    real eps2;                /* (potential softening parameter)^2 */

    real tree_rsize;
    real sunGCDist;
    criterion_t criterion;
    long seed;                /* random number seed */
    bool usequad;             /* use quadrupole corrections */
    bool allowIncest;
    bool outputCartesian;     /* print (x,y,z) instead of (l, b, r) */
    bool outputBodies;
    bool outputHistogram;

    const char* cp_filename;
    char cp_resolved[1024];
} NBodyCtx;

typedef int generic_enum_t;  /* A general enum type. */


/* Note: 'type' should first field for all types. */
#define SET_TYPE(x, y) (((Disk*)x)->type = y)
#define NBODY_TYPEOF(x) (((Disk*)x)->type)


/* Useful initializers */
#define EMPTY_SPHERICAL { 0, NAN, NAN }
#define EMPTY_DISK { 0, NAN, NAN, NAN }
#define EMPTY_HALO { 0, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN }
#define EMPTY_POTENTIAL { {EMPTY_SPHERICAL}, EMPTY_DISK, EMPTY_HALO, NULL }
#define EMPTY_MODEL { 0, 0, NAN, NAN, NAN, NAN, NAN, NAN, NAN }

#define EMPTY_TREE { NULL, NAN, 0, 0 }
#define EMPTY_VECTOR { NAN, NAN, NAN, NAN }
#define EMPTY_CTX { EMPTY_POTENTIAL, NULL, 0, 0, NAN, NAN, NAN, NAN, \
                    NULL, NULL, NULL, NULL, NULL,                    \
                    NAN, NAN, NAN, NAN, NAN, 0, 0,                   \
                    FALSE, FALSE, FALSE, FALSE, FALSE, NULL, "" }

#ifndef __OPENCL_VERSION__  /* No function pointers allowed in kernels */
/* Acceleration functions for a given potential */
typedef mwvector (*SphericalAccel) (const Spherical*, const mwvector);
typedef mwvector (*HaloAccel) (const Halo*, const mwvector);
typedef mwvector (*DiskAccel) (const Disk*, const mwvector);

/* Generic potential function */
typedef mwvector (*AccelFunc) (const void*, const mwvector);

#endif /* __OPENCL_VERSION__ */

#endif /* _NBODY_TYPES_H_ */

