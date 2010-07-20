/* ************************************************************************** */
/* DEFS.H: include file for hierarchical force calculation routines.  The */
/* definitions in this file are needed for load.c and grav.c; this file */
/* does not provide definitions for other parts of the N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

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


#ifndef __OPENCL_VERSION__   /* Not compiling CL kernel */
  #if NBODY_OPENCL
    #include <OpenCL/cl.h>
    #include <OpenCL/cl_platform.h>
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

#ifndef  DOUBLEPREC
  typedef float real, *realptr;
#else
  typedef double real, *realptr;
#endif /* DOUBLEPREC */

#define NDIM 3

/* Note: vectorptr is NOT the same as vector*.  By using real* as
   vectorptr, we can do nice things to avoid pointer aliasing and
   copying in various places. Use vector* only for mapping over an
   array of vectors.
 */

#if NBODY_OPENCL || defined(__OPENCL_VERSION__)
  #ifdef __OPENCL_VERSION__ /* In the kernel */
    #ifndef  DOUBLEPREC
      typedef float4 real4, *real4ptr;
    #else
      typedef double4 real4, *real4ptr;
    #endif /* DOUBLEPREC */
  #else
    #ifndef bool
      typedef int bool;
    #endif

    #ifndef  DOUBLEPREC
      typedef cl_float4 real4, *real4ptr;
    #else
      typedef cl_double4 real4, *real4ptr;
    #endif /* DOUBLEPREC */
  #endif /* __OPENCL_VERSION__ */

  typedef real4 vector;
  typedef real* vectorptr;

  typedef real4 matrix[NDIM];
  #define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
#else
  #ifndef bool
    typedef short int bool;
  #endif

  typedef real _vector_new __attribute__ ((vector_size (sizeof(real) * 4)));

  typedef union
  {
    _vector_new _v;
    real v[4];
  } vector_new;

  typedef real vector[NDIM];

  typedef real matrix[NDIM][NDIM];
  typedef real* vectorptr;

//#define ZERO_VECTOR { { 0.0, 0.0, 0.0, 0.0 } }
#endif /* NBODY_OPENCL */

#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }
#define ZERO_VECTOR { 0.0, 0.0, 0.0 }
/*
#define L(x) ((((vector)x).v)[0])
#define B(x) ((((vector)x).v)[1])
#define R(x) ((((vector)x).v)[2])
*/

#define XP(x) ((((vector_new)(x)).v)[0])
#define YP(x) ((((vector_new)(x)).v)[1])
#define ZP(x) ((((vector_new)(x)).v)[2])
#define VP(x) ((vector_new)(x)._v)



#define L(x) (((vectorptr) (x))[0])
#define B(x) (((vectorptr) (x))[1])
#define R(x) (((vectorptr) (x))[2])

#define X(x) (((vectorptr) (x))[0])
#define Y(x) (((vectorptr) (x))[1])
#define Z(x) (((vectorptr) (x))[2])


#ifndef TRUE
  #define TRUE  1
  #define FALSE 0
#endif

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
typedef struct _node
{
    body_t type;            /* code for node type */
    real mass;              /* total mass of node */
    vector pos;             /* position of node */
    struct _node* next;     /* link to next force-calc */
} node, *nodeptr;

#define Type(x) (((nodeptr) (x))->type)
#define Mass(x) (((nodeptr) (x))->mass)
#define Pos(x)  (((nodeptr) (x))->pos)
#define Next(x) (((nodeptr) (x))->next)

/* BODY: data structure used to represent particles. */

typedef struct
{
    node bodynode;              /* data common to all nodes */
    vector vel;                 /* velocity of body */
    vector acc;                 /* acceleration of body */
 } body, *bodyptr;

#define Body    body

#define Vel(x)  (((bodyptr) (x))->vel)

/* CELL: structure used to represent internal nodes of tree. */

#define NSUB (1 << NDIM)        /* subcells per cell */

typedef struct
{
    node cellnode;              /* data common to all nodes */
    real rcrit2;                /* critical c-of-m radius^2 */
    nodeptr more;               /* link to first descendent */
    union                       /* shared storage for... */
    {
        nodeptr subp[NSUB];     /* descendents of cell */
        matrix quad;            /* quad. moment of cell */
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
typedef struct
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

typedef struct
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

typedef struct
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

typedef struct
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

typedef struct
{
    cellptr root;   /* pointer to root cell */
    real rsize;     /* side-length of root cell */

    unsigned int cellused;   /* count of cells in tree */
    unsigned int maxlevel;   /* count of levels in tree */
} Tree;

typedef struct
{
    vector position;     /* (x, y, z) if cartesian / useGalC, otherwise (l, b, r) */
    vector velocity;
    bool useGalC;
    bool useRadians;
} InitialConditions;


typedef struct
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

typedef struct
{
    dwarf_model_t type;
    int nbody;
    real time_dwarf;
    real time_orbit;

    /* calculated depending on model */
    real timestep;
    real orbit_timestep;
    real eps;            /* potential softening parameter */

    /* model parameters */
    real mass;
    real scale_radius;
} DwarfModel;

#ifndef _WIN32

typedef struct
{
    int fd;            /* File descriptor for checkpoint file */
    char* mptr;        /* mmap'd pointer for checkpoint file */
    const char* filename;
} CheckpointHandle;

#define EMPTY_CHECKPOINT_HANDLE { -1, NULL, NULL }

#else

typedef struct
{
    HANDLE file;
    HANDLE mapFile;
    char* mptr;
    const char* filename;
} CheckpointHandle;

#define EMPTY_CHECKPOINT_HANDLE { INVALID_HANDLE_VALUE, INVALID_HANDLE_VALUE, NULL, NULL }

#endif /* _WIN32 */

#if NBODY_OPENCL && !defined(__OPENCL_VERSION__)

/* Host OpenCL stuff */
typedef struct
{
    cl_device_id dev;
    cl_device_type devType;
    unsigned int devCount;
    cl_context clctx;
    cl_command_queue queue;
    cl_program prog;
    cl_kernel kern;
} NBodyCLInfo;

#define EMPTY_NBODY_CL_INFO { -1, -1, 0, NULL, NULL, NULL, NULL }

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
typedef struct
{
    Tree tree;
    real tout;
    real tnow;
    bodyptr bodytab;    /* points to array of bodies */
    vector* acctab;     /* Corresponding accelerations of bodies */

  #if NBODY_OPENCL
    NBodyCLInfo ci;
    NBodyCLMem cm;
  #endif /* NBODY_OPENCL */
} NBodyState;

#if NBODY_OPENCL
  #define EMPTY_STATE { EMPTY_TREE, NAN, NAN, NULL, NULL, EMPTY_NBODY_CL_INFO, EMPTY_NBODY_CL_MEM }
#else
  #define EMPTY_STATE { EMPTY_TREE, NAN, NAN, NULL, NULL }
#endif /* NBODY_OPENCL */

#endif /* __OPENCL_VERSION__ */

/* The context tracks settings of the simulation.  It should be set
   once at the beginning of a simulation based on settings, and then
   stays constant for the actual simulation.
 */
typedef struct
{
    Potential pot;
    DwarfModel model;         /* dwarf model */
    char* headline;           /* message describing calculation */

    const char* outfilename;  /* output */
    const char* histogram;
    const char* histout;
    FILE* outfile;            /* file for snapshot output */
    CheckpointHandle cp;

    real freqout;
    real theta;               /* accuracy parameter: 0.0 */
    real tree_rsize;
    real sunGCDist;
    criterion_t criterion;
    long seed;                /* random number seed */
    bool usequad;             /* use quadrupole corrections */
    bool allowIncest;
    bool outputCartesian;     /* print (x,y,z) instead of (l, b, r) */
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
#define EMPTY_VECTOR { NAN, NAN, NAN }
#define EMPTY_INITIAL_CONDITIONS { EMPTY_VECTOR, EMPTY_VECTOR, FALSE, FALSE }
#define EMPTY_CTX { EMPTY_POTENTIAL, EMPTY_MODEL, NULL, NULL, NULL, NULL, NULL, \
                     EMPTY_CHECKPOINT_HANDLE, NAN, NAN, NAN, NAN, 0, 0, FALSE, FALSE, FALSE }

#ifndef __OPENCL_VERSION__  /* No function pointers allowed in kernels */
/* Acceleration functions for a given potential */
typedef void (*SphericalAccel) (vectorptr restrict, const Spherical*, const vectorptr restrict);
typedef void (*HaloAccel) (vectorptr restrict, const Halo*, const vectorptr restrict);
typedef void (*DiskAccel) (vectorptr restrict, const Disk*, const vectorptr restrict);

/* Generic potential function */
typedef void (*AccelFunc) (vectorptr restrict, const void*, const vectorptr restrict);

#endif /* __OPENCL_VERSION__ */

#endif /* _NBODY_TYPES_H_ */

