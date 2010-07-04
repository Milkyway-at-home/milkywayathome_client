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

#include "stdinc.h"
#include "real.h"
#include "vectmath.h"

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

/*
typedef enum
{
    BODY,
    CELL
} body_t; */

#define BODY 01
#define CELL 02

typedef short body_t;

/* NODE: data common to BODY and CELL structures. */

typedef struct _node
{
    body_t type;             /* code for node type */
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
#define Acc(x)  (((bodyptr) (x))->acc)

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

typedef enum
{
    SphericalPotential
} spherical_t;

/* Spherical potential */
typedef struct
{
    spherical_t type;
    real mass;
    real scale;
} Spherical;


/* Supported disk models */
typedef enum
{
    MiyamotoNagaiDisk,
    ExponentialDisk
} disk_t;

typedef struct
{
    disk_t type;
    real mass;          /* disk mass */
    real scale_length;  /* "a" for M-N, "b" for exp disk */
    real scale_height;  /* unused for exponential disk. "b" for Miyamoto-Nagai disk */
} Disk;

/* Supported halo models */
typedef enum
{
    LogarithmicHalo,
    NFWHalo,
    TriaxialHalo
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

    int cellused;   /* count of cells in tree */
    int maxlevel;   /* count of levels in tree */
} Tree;

typedef struct
{
    vector position;     /* (x, y, z) if cartesian / useGalC, otherwise (l, b, r) */
    vector velocity;
    bool useGalC;
    bool useRadians;
} InitialConditions;

#define L(x) (((vectorptr) (x))[0])
#define B(x) (((vectorptr) (x))[1])
#define R(x) (((vectorptr) (x))[2])

#define X(x) (((vectorptr) (x))[0])
#define Y(x) (((vectorptr) (x))[1])
#define Z(x) (((vectorptr) (x))[2])

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


/* The context tracks settings of the simulation.  It should be set
   once at the beginning of a simulation based on settings, and then
   stays constant for the actual simulation.
 */
typedef struct
{
    Potential pot;
    DwarfModel model;         /* dwarf model */
    char* headline;           /* message describing calculation */

    const char* outfilename;  /* filename for snapshot output */
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

/* Mutable state used during an evaluation */
typedef struct
{
    Tree tree;
    real tout;
    real tnow;
    bodyptr bodytab;    /* points to array of bodies */
} NBodyState;

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
#define EMPTY_CTX { EMPTY_POTENTIAL, EMPTY_MODEL, NULL, NULL, NULL, EMPTY_CHECKPOINT_HANDLE, NAN, NAN, NAN, NAN, 0, 0, FALSE, FALSE, FALSE }
#define EMPTY_TREE { NULL, NAN, 0, 0 }
#define EMPTY_STATE { EMPTY_TREE, NAN, NAN, NULL}
#define EMPTY_VECTOR { NAN, NAN, NAN }
#define EMPTY_INITIAL_CONDITIONS { EMPTY_VECTOR, EMPTY_VECTOR, FALSE, FALSE }


/* Acceleration functions for a given potential */
typedef void (*SphericalAccel) (vectorptr restrict, const Spherical*, const vectorptr restrict);
typedef void (*HaloAccel) (vectorptr restrict, const Halo*, const vectorptr restrict);
typedef void (*DiskAccel) (vectorptr restrict, const Disk*, const vectorptr restrict);

/* Generic potential function */
typedef void (*AccelFunc) (vectorptr restrict, const void*, const vectorptr restrict);


/* PROC, IPROC: pointers to procedures and integer-valued functions. */
typedef void (*proc)();
typedef int (*iproc)();

#endif /* _NBODY_TYPES_H_ */

