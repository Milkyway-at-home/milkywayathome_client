/* ************************************************************************** */
/* DEFS.H: include file for hierarchical force calculation routines.  The */
/* definitions in this file are needed for load.c and grav.c; this file */
/* does not provide definitions for other parts of the N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _DEFS_H_
#define _DEFS_H_

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

/* NODE: data common to BODY and CELL structures. */

typedef struct _node
{
    short type;             /* code for node type */
    real mass;              /* total mass of node */
    vector pos;             /* position of node */
    struct _node* next;     /* link to next force-calc */
} node, *nodeptr;

#define Type(x) (((nodeptr) (x))->type)
#define Mass(x) (((nodeptr) (x))->mass)
#define Pos(x)  (((nodeptr) (x))->pos)
#define Next(x) (((nodeptr) (x))->next)

/*  * BODY: data structure used to represent particles.
 */

#define BODY 01                 /* type code for bodies */

typedef struct
{
    node bodynode;              /* data common to all nodes */
    vector vel;                 /* velocity of body */
    vector acc;                 /* acceleration of body */
    real phi;                   /* potential at body */
} body, *bodyptr;

#define Body    body

#define Vel(x)  (((bodyptr) (x))->vel)
#define Acc(x)  (((bodyptr) (x))->acc)
#define Phi(x)  (((bodyptr) (x))->phi)

/* CELL: structure used to represent internal nodes of tree. */

#define CELL 02                 /* type code for cells */

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
    BH86,
    SW93
} criterion_t;

typedef enum
{
    SphericalPotential = 1 << 1
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
    MiaymotoNagaiDisk = 1 << 2,
    ExponentialDisk   = 1 << 3
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
    LogarithmicHalo = 1 << 4,
    NFWHalo         = 1 << 5,
    TriaxialHalo    = 1 << 6
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

/* Parameters and results for gravitational calculation. */

typedef struct
{
    real PluMass, r0;
    real lstart, bstart, Rstart;
    real XC, YC, ZC;
    real VXC, VYC, VZC;
    real Xinit, Yinit, Zinit;
    real VXinit, VYinit, VZinit;
    real orbittstop, dtorbit;
    real sunGCDist;
} NBodyParams;

typedef struct
{
    bool useGalC;
    vector position;     /* (x, y, z) if cartesian / useGalC, otherwise (r, l, b) */
    vector velocity;
} InitialConditions;

typedef enum
{
    DwarfModelPlummer = 1 << 1,
    DwarfModelKing    = 1 << 2,
    DwarfModelDehnen  = 1 << 3
} dwarf_model_t;

typedef struct
{
    dwarf_model_t type;
    real timestep;        /* calculated from other information */
    real time_dwarf;
    real time_orbit;
    real mass;
    real scale_radius;
} DwarfModel;


/* The context tracks settings of the simulation.  It should be set
   once at the beginning of a simulation based on settings, and then
   stays constant for the actual simulation.
 */
typedef struct
{
    /* Actual parameters */
    Potential pot;
    DwarfModel model;   /* dwarf model */
    int nbody;

    /* TODO: these should go away */
    real tstop;
    real dtout;
    real freq;
    real freqout;

    /* Simulation settings */
    bool usequad;          /* use quadrupole corrections */
    bool allowIncest;
    criterion_t criterion;         /* bh86 or sw93 */
    real theta;     /* accuracy parameter: 0.0 => exact */
    real eps;       /* potential softening parameter */

    /* Utilitarian type information */
    int seed;             /* random number seed */
    FILE* outfile;        /* file for snapshot output */
    char* outfilename;    /* filename for snapshot output */
    char* headline;       /* message describing calculation */
} NBodyCtx;

/* Mutable state used during an evaluation */
typedef struct
{
    int n2bterm;    /* number 2-body of terms evaluated */
    int nbcterm;    /* num of body-cell terms evaluated */
    int nstep;      /* number of time-steps */

    real tout;
    real tnow;
    bodyptr bodytab;      /* points to array of bodies */
} NBodyState;

typedef int generic_enum_t;  /* A general enum type. */

/* Note: 'type' should first field for all types. */
#define SET_TYPE(x, y) (((Disk*)x)->type = y)
#define NBODY_TYPEOF(x) (((Disk*)x)->type)


#define fail(msg, ...) { fprintf(stderr, msg, ##__VA_ARGS__);  \
                         exit(EXIT_FAILURE); }



/* Useful initializers */
#define EMPTY_SPHERICAL { 0, NAN, NAN }
#define EMPTY_DISK { 0, NAN, NAN, NAN }
#define EMPTY_HALO { 0, NAN, NAN, NAN, NAN, NAN, NAN }
#define EMPTY_POTENTIAL { {EMPTY_SPHERICAL}, EMPTY_DISK, EMPTY_HALO, NULL }
#define EMPTY_MODEL { 0, NAN, NAN, NAN, NAN }
#define EMPTY_CTX { EMPTY_POTENTIAL,  EMPTY_MODEL, 0, NAN, NAN, NAN, NAN, FALSE, FALSE, 0, NAN, NAN, 0, NULL, NULL, NULL }

#define EMPTY_TREE { NULL, NAN, 0, 0 }

#define EMPTY_STATE { 0, 0, 0, NAN, NAN, NULL }
#define EMPTY_INITIAL_CONDITIONS { FALSE, { NAN, NAN, NAN }, { NAN, NAN, NAN } }

/* Utility routines used in load.c and grav.c.  These are defined in
 * util.c, which must be compiled with the same choice of precision.
 */

real cputime(void);              /* return elapsed CPU time */
void* allocate(int);             /* allocate and zero memory */
real distv(vector, vector);      /* distance between vectors */
void error(char*, ...);          /* report error and exit */
void eprintf(char*, ...);        /* printf to error FILE* */
int compare (const void* a, const void* b);     /* comparison function used in chisq */
float chisq();                  /* likelihood calculator */

#endif /* _DEFS_H_ */

