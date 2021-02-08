/*
  Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
  Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright (c) 2016-2018 Siddhartha Shelton

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
#define idBody(x) (((NBodyNode*) (x))->id)
#define bodyTypeIsIgnore(x) ((x) < 0)

typedef short body_t;

/* node: data common to BODY and CELL structures. */
typedef struct MW_ALIGN_TYPE _NBodyNode
{
    mwvector pos;             /* position of node */
    struct _NBodyNode* next;  /* link to next force-calc */
    real mass;                /* total mass of node */
    body_t type;              /* code for node type */
    unsigned int id;          /* body id */
} NBodyNode;

#define EMPTY_NODE { ZERO_VECTOR, NULL, 0.0, 0, 0  }

#define Type(x) (((NBodyNode*) (x))->type)
#define Mass(x) (((NBodyNode*) (x))->mass)
#define Pos(x)  (((NBodyNode*) (x))->pos)
#define Next(x) (((NBodyNode*) (x))->next)

/* BODY: data structure used to represent particles. */

typedef struct MW_ALIGN_TYPE
{
    NBodyNode bodynode;         /* data common to all nodes */
    mwvector vel;               /* velocity of body */
} Body;

#define EMPTY_BODY { EMPTY_NODE, ZERO_VECTOR }

#define BODY_TYPE "Body"

#define Vel(x)  (((Body*) (x))->vel)

/* CELL: structure used to represent internal nodes of tree. */

#define NSUB (1 << NDIM)        /* subcells per cell */

typedef struct MW_ALIGN_TYPE
{
    real xx, xy, xz;
    real yy, yz;
    real zz;
} NBodyQuadMatrix;



typedef struct MW_ALIGN_TYPE
{
    NBodyNode cellnode;         /* data common to all nodes */
    real rcrit2;                /* critical c-of-m radius^2 */
    NBodyNode* more;            /* link to first descendent */
    union MW_ALIGN_V(16)        /* shared storage for... */
    {
        NBodyNode* subp[NSUB];  /* descendents of cell */
        NBodyQuadMatrix quad;   /* quad. moment of cell. Unique symmetric matrix components */
    } stuff;
} NBodyCell;

/* use alternate criteria */
typedef enum
{
    InvalidCriterion = InvalidEnum,
    TreeCode, 
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

typedef struct MW_ALIGN_TYPE
{
    NBodyCell* root;         /* pointer to root cell */
    real rsize;              /* side-length of root cell */

    unsigned int cellUsed;   /* count of cells in tree */
    unsigned int maxDepth;   /* count of levels in tree */
    int structureError;
} NBodyTree;

#define EMPTY_TREE { NULL, 0.0, 0, 0, 0 }


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

    /* Just valid non-aliasing, read only buffers.
     * Used as dummy arguments for kernel arguments we don't need
     * depending on the specific simulation options since you can't set
     * a kernel argument to null.
     *
     * The AMD runtime will rebuild the kernel with noalias if it
     * thinks an argument aliases, so use this to avoid that.
     *
     * This is probably less of a mess than the work required to
     * account for fewer / a different order of arguments in different
     * cases.
     */
    cl_mem dummy[18];
} NBodyBuffers;


/* 8 used by tree + 1 with quad, 2 used by exact. 1 shared. */
#define NKERNELS 11


typedef struct
{
    cl_kernel boundingBox;
    cl_kernel buildTreeClear;
    cl_kernel buildTree;
    cl_kernel summarizationClear;
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
    size_t factors[8];
    size_t threads[8];
    real timings[8];        /* In a single iteration */
    real chunkTimings[8];   /* Average time per chunk */
    real kernelTimings[8];  /* Running totals */

    size_t global[8];
    size_t local[8];
} NBodyWorkSizes;


typedef struct
{
    int useBin;
    unsigned int rawCount;
    real lambda;
    real beta;
    real variable;  // can be customized to any of the different histograms
    real err;       // the error in variable
    real sum;       // used by beta & vel disp
    real sq_sum;    // used by beta & vel disp
    real outliersRemoved;   
    
} HistData;


typedef struct
{
    /* Center Point */
    real phi;
    real theta;
    real psi;

    /* Lambda */
    real lambdaStart;
    real lambdaEnd;
    unsigned int lambdaBins;

    /* Beta */
    real betaStart;
    real betaEnd;
    unsigned int betaBins;

} HistogramParams;

#define EMPTY_HISTOGRAM_PARAMS { 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0 }
#define HISTOGRAM_PARAMS_TYPE "HistogramParams"

typedef struct
{
    unsigned int lambdaBins;
    unsigned int betaBins;
    unsigned int totalNum;
    unsigned int totalSimulated;
    int hasRawCounts;
    HistogramParams params;
    real massPerParticle;

    /* This is used as a variable length struct. Do not add any fields
     * after data. */
    HistData data[1];
} NBodyHistogram;

/* class to hold all of the different histograms */
typedef struct
{
    // data input must ALWAYS be in the following order:
    // 0: normalized counts
    // 1: beta disp
    // 2: vel disp
    // 3: vlos avg
    // 4: beta avg
    // 5: dist avg
    mwbool usage[6];
    NBodyHistogram* histograms[6]; 
} MainStruct;

/* Mutable state used during an evaluation */
typedef struct MW_ALIGN_TYPE
{
    NBodyTree tree;
    NBodyNode* freeCell;      /* list of free cells */
    char* checkpointResolved;
    Body* bodytab;            /* points to array of bodies */
    mwvector* acctab;         /* Corresponding accelerations of bodies */
    mwvector* orbitTrace;     /* Trail of center of masses for display purposes */
    scene_t* scene;

    mwvector* shiftByLMC;      /* Accelerations on MW from LMC */
    mwvector* LMCpos;        /* Position of LMC */
    mwvector* LMCvel;        /* Velocity of LMC */

    lua_State** potEvalStates;  /* If using a Lua closure as a potential, the evaluation states.
                                   We need one per thread in the general case. */
    int* potEvalClosures;       /* Lua closure for each state */

    size_t nOrbitTrace;         /* Number of items in orbitTrace */
    size_t nShiftLMC;          /* Size of LMC Acceleration array */
    time_t lastCheckpoint;

    unsigned int step;
    int nbody;
    int effNBody;            /* Sometimes needed rounded up number of bodies. >= nbody are just padding */
    int treeIncest;          /* Tree incest has occured */
    int potentialEvalError;  /* Error occured in calling custom Lua potential */

    unsigned int maxDepth;   /* Maximum depth before overflow. Used for CL version */
    
    real bestLikelihood;           /* new parameter for best likelihood eval */
    real bestLikelihood_EMD;       /* EMD component of likelihood */
    real bestLikelihood_Mass;      /* Mass component of likelihood */
    real bestLikelihood_Beta;      /* Beta disp component of likelihood */
    real bestLikelihood_Vel;       /* Velocity disp component of likelihood */
    real bestLikelihood_BetaAvg;   /* Beta avg component of likelihood */
    real bestLikelihood_VelAvg;    /* Velocity avg component of likelihood */
    real bestLikelihood_Dist;      /* Distance component of likelihood */
    real bestLikelihood_time;      /* to store the evolve time at which the best likelihood occurred */
    int bestLikelihood_count;      /* count of how many times the likelihood improved */
    mwbool useVelDisp;             /* whether or not to use the vel disp comparison */
    mwbool useBetaDisp;            /* whether or not to use the beta disp comparison */
    mwbool useBetaComp;            /* whether or not to use the avg beta comparison */
    mwbool useVlos;                /* whether or not to use the avg vlos comparison */
    mwbool useDist;                /* whether or not to use the avg distance comparison */

    mwbool ignoreResponsive;
    mwbool usesExact;
    mwbool usesQuad;
    mwbool usesConsistentMemory;
    mwbool dirty;      /* Whether the view of the bodies is consistent with the view in the CL buffers */
    mwbool usesCL;
    mwbool useCLCheckpointing;
    mwbool reportProgress;

    
    real previousForwardTime;   //used to calibrate bar time
    
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

#define EMPTY_NBODYSTATE { EMPTY_TREE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, \
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, \
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 0, NULL,\
NULL, NULL, NULL}


/* The context tracks settings of the simulation.  It should be set
   once at the beginning of a simulation based on settings, and then
   stays constant for the actual simulation.
 */
typedef struct MW_ALIGN_TYPE
{
    real eps2;                /* (potential softening parameter)^2 */
    real theta;               /* accuracy parameter: 0.0 */
    real timestep;
    real timeEvolve;
    real timeBack;
    real treeRSize;
    real sunGCDist;

    real b;     /* orbital parameters */
    real r;
    real vx;
    real vy;
    real vz;

    criterion_t criterion;
    ExternalPotentialType potentialType;
    
    mwbool Nstep_control;     /* manually control how many timesteps simulation runs */
    mwbool useBestLike;       /* use best likelihood return code */
    mwbool useVelDisp;        /* use the velocity dispersion comparison calc */
    mwbool useBetaDisp;       /* use the beta dispersion comparison calc */
    mwbool useBetaComp;       /* use the beta average comparison calc */
    mwbool useVlos;           /* use the line of sight velocity comparison calc */
    mwbool useDist;           /* use the average distance comparison calc */
    mwbool MultiOutput;       /* whether to have algorithm put out multiple outputs */
    
    mwbool useQuad;           /* use quadrupole corrections */
    mwbool allowIncest;
    mwbool quietErrors;
    
    real BestLikeStart;       /* after what portion of the sim should the calc start */
    real OutputFreq;          /* frequency of writing outputs */
    
    real BetaSigma;           /* sigma cutoff for the outlier rejection for the bin beta dispersions */ 
    real VelSigma;            /* sigma cutoff for the outlier rejection for the bin vel dispersions */ 
    real DistSigma;            /* sigma cutoff for the outlier rejection for the bin dists dispersions */ 
    real IterMax;             /* number of times to apply outlier rejection with sigma cutoff */ 
    real BetaCorrect;         /* correction factor for correcting the distribution after outlier rejection */
    real VelCorrect;          /* correction factor for correcting the distribution after outlier rejection */
    real DistCorrect;          /* correction factor for correcting the distribution after outlier rejection */
    mwbool LMC;

    real LMCmass;              /* Mass of LMC */
    real LMCscale;             /* Scale radius of LMC */
    mwbool LMCDynaFric;        /* LMC Dynamical Friction switch */

    unsigned int calibrationRuns; //for calibrating time-dependent potentials

    real Ntsteps;              /* number of time steps to run when manual control is on */
    time_t checkpointT;        /* Period to checkpoint when not using BOINC */
    unsigned int nStep;

    Potential pot;


} NBodyCtx;

#define NBODYCTX_TYPE "NBodyCtx"
#define EMPTY_NBODYCTX { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                  \
                         InvalidCriterion, EXTERNAL_POTENTIAL_DEFAULT,                                \
                         FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, \
                         0, 0, 0, 0, 0, 0, 0, 0, 0, FALSE,                                           \
                         0, 0, FALSE, 0,                                                                \
                         0, 0, 0,                                                                     \
                         EMPTY_POTENTIAL}

/* Negative codes can be nonfatal but useful return statuses.
   Positive can be different hard failures.
 */
typedef enum
{
    NBODY_SUCCESS              = 0 << 0,
    NBODY_ERROR                = 1 << 1,
    NBODY_ASSERTION_FAILURE    = 1 << 2,
    NBODY_TREE_STRUCTURE_ERROR = 1 << 3,
    NBODY_TREE_INCEST_FATAL    = 1 << 4,
    NBODY_IO_ERROR             = 1 << 5,
    NBODY_CHECKPOINT_ERROR     = 1 << 6,
    NBODY_CL_ERROR             = 1 << 7,
    NBODY_CAPABILITY_ERROR     = 1 << 8,
    NBODY_CONSISTENCY_ERROR    = 1 << 9,
    NBODY_UNIMPLEMENTED        = 1 << 10,
    NBODY_UNSUPPORTED          = 1 << 11,
    NBODY_USER_ERROR           = 1 << 12,
    NBODY_PARAM_FILE_ERROR     = 1 << 13,
    NBODY_LUA_POTENTIAL_ERROR  = 1 << 14,
    NBODY_LIKELIHOOD_ERROR     = 1 << 15,
    NBODY_MAX_DEPTH_ERROR      = 1 << 16,
    NBODY_CELL_OVERFLOW_ERROR  = 1 << 17,
    NBODY_GRAPHICS_TIMEOUT     = 1 << 18,
    NBODY_GRAPHICS_DEAD        = 1 << 19,
    NBODY_RESERVED_ERROR_1     = 1 << 20,
    NBODY_RESERVED_ERROR_2     = 1 << 21,
    NBODY_RESERVED_ERROR_3     = 1 << 22,
    NBODY_RESERVED_ERROR_4     = 1 << 23,

    /* Warnings */
    NBODY_TREE_INCEST_NONFATAL = 1 << 24,
    NBODY_RESERVED_WARNING_1   = 1 << 25,
    NBODY_RESERVED_WARNING_2   = 1 << 26,
    NBODY_RESERVED_WARNING_3   = 1 << 27,
    NBODY_RESERVED_WARNING_4   = 1 << 28,
    NBODY_RESERVED_WARNING_5   = 1 << 29,
    NBODY_RESERVED_WARNING_6   = 1 << 30,
    NBODY_RESERVED_WARNING_7   = 1 << 31
} NBodyStatus;

#define NBODY_STATUS_ALL_WARNINGS (NBODY_TREE_INCEST_NONFATAL | NBODY_RESERVED_WARNING_1 | NBODY_RESERVED_WARNING_2 | NBODY_RESERVED_WARNING_3 | NBODY_RESERVED_WARNING_4 | NBODY_RESERVED_WARNING_5 | NBODY_RESERVED_WARNING_6 | NBODY_RESERVED_WARNING_7)

/* Right now there is only one warning */
#define nbStatusIsFatal(x) (((x) & ~NBODY_STATUS_ALL_WARNINGS) != 0)
#define nbStatusIsWarning(x) (((x) & NBODY_STATUS_ALL_WARNINGS) != 0)

/* Reserve positive numbers for reporting depth > MAXDEPTH. Should match in kernel  */
typedef enum
{
    NBODY_KERNEL_OK                   = 0,
    NBODY_KERNEL_CELL_OVERFLOW        = -1,
    NBODY_KERNEL_TREE_INCEST          = -2,
    NBODY_KERNEL_TREE_STRUCTURE_ERROR = -3,
    NBODY_KERNEL_ERROR_OTHER          = -4
} NBodyKernelError;


/* Note: 'type' should first field for all types. */
#define SET_TYPE(x, y) (((Disk*)x)->type = y)
#define NBODY_TYPEOF(x) (((Disk*)x)->type)



typedef enum
{
    NBODY_INVALID_METHOD = -1,
    NBODY_EMD,
    NBODY_ORIG_CHISQ,
    NBODY_ORIG_ALT,
    NBODY_CHISQ_ALT,
    NBODY_POISSON,
    NBODY_KOLMOGOROV,
    NBODY_KULLBACK_LEIBLER,
    NBODY_SAHA
} NBodyLikelihoodMethod;

typedef struct MW_ALIGN_TYPE
{
    mwvector revOrbitPos;
    mwvector revOrbitVel;
    mwvector revOrbitLMCPos;
    mwvector revOrbitLMCVel;
    real revOrbitdt;
    real LMCmass;
    real revOrbitTstop;
    Potential pot;
    real previousForwardTime;

}SingleParticleOrbitParams;

NBodyStatus nbInitCL(NBodyState* st, const NBodyCtx* ctx, const CLRequest* clr);
NBodyStatus nbInitNBodyStateCL(NBodyState* st, const NBodyCtx* ctx);

int destroyNBodyState(NBodyState* st);
int nbDetachSharedScene(NBodyState* st);
void setLMCShiftArray(NBodyState* st, mwvector* shiftArray, size_t shiftSize);
void setLMCPosVel(NBodyState* st, mwvector* PosArray, mwvector* VelArray);
void setInitialNBodyState(NBodyState* st, const NBodyCtx* ctx, Body* bodies, int nbody);
void setRandomLMCNBodyState(NBodyState* st, int nShift, dsfmt_t* dsfmtState);
void cloneNBodyState(NBodyState* st, const NBodyState* oldSt);
int equalNBodyState(const NBodyState* st1, const NBodyState* st2);

void sortBodies(Body* bodies, int nbody);

int equalSpherical(const Spherical* s1, const Spherical* s2);
int equalHalo(const Halo* h1, const Halo* h2);
int equalDisk(const Disk* d1, const Disk* d2);

// int equalDwarf(const Dwarf* h1, const Dwarf* h2);

int equalPotential(const Potential* p1, const Potential* p2);

int equalNBodyCtx(const NBodyCtx* ctx1, const NBodyCtx* ctx2);

int equalHistogramParams(const HistogramParams* hp1, const HistogramParams* hp2);


#endif /* _NBODY_TYPES_H_ */

