/*
 * Copyright (c) 2010 The University of Texas at Austin
 * Copyright (c) 2010 Dr. Martin Burtscher
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
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

#ifndef DOUBLEPREC
  #error Precision not defined
#endif


#if !BH86 && !SW93 && !NEWCRITERION && !EXACT
  #error Opening criterion not set
#endif

/* Reserve positive numbers for reporting depth > MAXDEPTH */
typedef enum
{
    NBODY_KERNEL_OK                   = 0,
    NBODY_KERNEL_CELL_LEQ_NBODY       = -1,
    NBODY_KERNEL_TREE_INCEST          = -2,
    NBODY_KERNEL_TREE_STRUCTURE_ERROR = -3,
    NBODY_KERNEL_ERROR_OTHER          = -4
} NBodyKernelError;

#if DOUBLEPREC
  #if cl_amd_fp64
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
  #elif cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
  #else
    #error Missing double precision extension
  #endif
#endif /* DOUBLEPREC */

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable



#if 0
#pragma OPENCL EXTENSION cl_amd_printf : enable
#define assert(x) if (!(x)) { printf("Assertion failed[%d]: %s\n", __LINE__, #x); }
#define dprintf(fmt, ...) printf("Line %d: " fmt, __LINE__, __VA_ARGS__)

#define ddprintf(fmt, ...) printf("[%d][%d][%d]: Line %d: " fmt, (int) get_global_id(0), (int) get_group_id(0), (int) get_local_id(0),__LINE__, __VA_ARGS__)
#define __BARRIER(type) \
    dprintf("\t__BARRIER[%d] hit by thread %d\n", __LINE__, get_local_id(0)); \
  barrier(type);

#else
#define assert(x)
#define dprintf(fmt, ...)
#define __BARRIER(type) barrier((type))
#define ddprintf(fmt, ...)
#endif


#if DOUBLEPREC
typedef double real;
#else
typedef float real;
#endif /* DOUBLEPREC */


/* FIXME: This needs to be the same as on the host */
typedef struct __attribute__((aligned))
{
    volatile real radius;
    volatile int bottom;
    volatile int maxDepth;
    volatile int errorCode;
    volatile unsigned int blkCnt;
} TreeStatus;

#define NSUB 8

typedef struct
{
    volatile real f[32];
    volatile int i[64];
} Debug;

#define isBody(n) ((n) < NBODY)
#define isCell(n) ((n) >= NBODY)


typedef __global volatile real* restrict RVPtr;
typedef __global volatile int* restrict IVPtr;

//typedef __global real* restrict RVPtr;
//typedef __global int* restrict IVPtr;


/* All kernels will use the same parameters for now */
#define NBODY_KERNEL(name) name(                        \
    RVPtr _posX, RVPtr _posY, RVPtr _posZ,              \
    RVPtr _velX, RVPtr _velY, RVPtr _velZ,              \
    RVPtr _accX, RVPtr _accY, RVPtr _accZ,              \
                                                        \
    RVPtr _maxX, RVPtr _maxY, RVPtr _maxZ,              \
    RVPtr _minX, RVPtr _minY, RVPtr _minZ,              \
                                                        \
    RVPtr _mass,                                        \
    IVPtr _start, IVPtr _count,                         \
    IVPtr _child, IVPtr _sort,                          \
    __global volatile TreeStatus* restrict _treeStatus, \
                                                        \
    int step,                                           \
    RVPtr _critRadii,                                   \
    __global volatile Debug* _debug                     \
    )


__attribute__ ((reqd_work_group_size(THREADS1, 1, 1)))
__kernel void NBODY_KERNEL(boundingBox)
{
    __local volatile real minX[THREADS1], minY[THREADS1], minZ[THREADS1];
    __local volatile real maxX[THREADS1], maxY[THREADS1], maxZ[THREADS1];

    int i = (int) get_local_id(0);
    if (i == 0)
    {
        minX[0] = _posX[0];
        minY[0] = _posY[0];
        minZ[0] = _posZ[0];
    }
    barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

    /* initialize with valid data (in case #bodies < #threads) */
    minX[i] = maxX[i] = minX[0];
    minY[i] = maxY[i] = minY[0];
    minZ[i] = maxZ[i] = minZ[0];

    int inc = get_local_size(0) * get_num_groups(0);
    int j = i + get_group_id(0) * get_local_size(0); // = get_global_id(0);
    while (j < NBODY) /* Scan bodies */
    {
        real tmp = _posX[j];
        minX[i] = min(minX[i], tmp);
        maxX[i] = max(maxX[i], tmp);

        tmp = _posY[j];
        minY[i] = min(minY[i], tmp);
        maxY[i] = max(maxY[i], tmp);

        tmp = _posZ[j];
        minZ[i] = min(minZ[i], tmp);
        maxZ[i] = max(maxZ[i], tmp);

        j += inc;  /* Move on to next body */
    }

    /* Reduction in shared memory */
    j = get_local_size(0) >> 1;
    while (j > 0)
    {
        barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
        if (i < j)
        {
            minX[i] = min(minX[i], minX[i + j]);
            minY[i] = min(minY[i], minY[i + j]);
            minZ[i] = min(minZ[i], minZ[i + j]);

            maxX[i] = max(maxX[i], maxX[i + j]);
            maxY[i] = max(maxY[i], maxY[i + j]);
            maxZ[i] = max(maxZ[i], maxZ[i + j]);
        }

        j >>= 1;
    }

    if (i == 0)
    {
        /* Write block result to global memory */
        j = get_group_id(0);

        _minX[j] = minX[0];
        _minY[j] = minY[0];
        _minZ[j] = minZ[0];

        _maxX[j] = maxX[0];
        _maxY[j] = maxY[0];
        _maxZ[j] = maxZ[0];
        mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

        inc = get_num_groups(0) - 1;
        if (inc == atom_inc(&_treeStatus->blkCnt))
        {
            /* I'm the last block, so combine all block results */
            for (j = 0; j <= inc; ++j)
            {
                minX[0] = min(minX[0], _minX[j]);
                minY[0] = min(minY[0], _minY[j]);
                minZ[0] = min(minZ[0], _minZ[j]);

                maxX[0] = max(maxX[0], _maxX[j]);
                maxY[0] = max(maxY[0], _maxY[j]);
                maxZ[0] = max(maxZ[0], _maxZ[j]);
            }

            /* Compute radius */
            real tmpR = max(maxX[0] - minX[0], maxY[0] - minY[0]);
            real radius = 0.5 * max(tmpR, maxZ[0] - minZ[0]);

            _treeStatus->radius = radius;

          #if NEWCRITERION || SW93
            _critRadii[NNODE] = radius;
          #endif


            /* Create root node */
            _mass[NNODE] = -1.0;
            _start[NNODE] = 0;
            _posX[NNODE] = 0.5 * (minX[0] + maxX[0]);
            _posY[NNODE] = 0.5 * (minY[0] + maxY[0]);
            _posZ[NNODE] = 0.5 * (minZ[0] + maxZ[0]);

            #pragma unroll NSUB
            for (int k = 0; k < NSUB; ++k)
                _child[NSUB * NNODE + k] = -1;

            _treeStatus->bottom = NNODE;
            _treeStatus->blkCnt = 0;  /* If this isn't 0'd for next time, everything explodes */
        }
    }
}

__attribute__ ((reqd_work_group_size(THREADS2, 1, 1)))
__kernel void NBODY_KERNEL(buildTree)
{
    __local real radius, rootX, rootY, rootZ;
    __local volatile int deadCount;

    if (get_local_id(0) == 0)
    {
        /* Cache root data */
        radius = _treeStatus->radius;
        rootX = _posX[NNODE];
        rootY = _posY[NNODE];
        rootZ = _posZ[NNODE];

        deadCount = 0;
    }
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    int localMaxDepth = 1;
    int skip = 1;
    int inc = get_local_size(0) * get_num_groups(0);
    int i = get_global_id(0);


    bool dead = false;
    while (1)  /* while (i < NBODY) */
    {

        /* Ugly hackery to prevent conditional barrier() for when some
         * items have another body and others don't */
        if (i >= NBODY && !dead)
        {
            dead = true;
            (void) atom_inc(&deadCount);
        }

        /* Wait for other wavefronts to finish loading */
        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

        if (deadCount == THREADS2)
            break;

        if (dead)
            continue;

        real r;
        real px, py, pz;
        int j, n, depth;

        if (skip != 0)
        {
            /* New body, so start traversing at root */
            skip = 0;

            px = _posX[i];
            py = _posY[i];
            pz = _posZ[i];
            n = NNODE;
            depth = 1;
            r = radius;

            /* Determine which child to follow */
            j = 0;
            if (rootX <= px)
                j = 1;
            if (rootY <= py)
                j += 2;
            if (rootZ <= pz)
                j += 4;
        }

        int ch = _child[NSUB * n + j];
        while (ch >= NBODY)  /* Follow path to leaf cell */
        {
            n = ch;
            ++depth;
            r *= 0.5;

            /* Determine which child to follow */
            j = 0;
            if (_posX[n] <= px)
                j = 1;
            if (_posY[n] <= py)
                j += 2;
            if (_posZ[n] <= pz)
                j += 4;
            ch = _child[NSUB * n + j];
        }

        if (ch != -2) /* Skip if child pointer is locked and try again later */
        {
            int locked = NSUB * n + j;
            if (ch == atom_cmpxchg(&_child[locked], ch, -2)) /* Try to lock */
            {
                if (ch == -1)
                {
                    /* If null, just insert the new body */
                    _child[locked] = i;
                }
                else  /* There already is a body in this position */
                {
                    int patch = -1;
                    /* Create new cell(s) and insert the old and new body */
                    do
                    {
                        ++depth;

                        int cell = atom_dec(&_treeStatus->bottom) - 1;
                        if (cell <= NBODY)
                        {
                            _treeStatus->errorCode = NBODY_KERNEL_CELL_LEQ_NBODY;
                            _treeStatus->bottom = NNODE;
                        }
                        patch = max(patch, cell);

                        _mass[cell] = -1.0;
                        _start[cell] = -1;

                      #if SW93 || NEWCRITERION
                        _critRadii[cell] = r;  /* Save cell size */
                      #endif /* SW93 || NEWCRITERION */


                        real nx = _posX[n];
                        real ny = _posY[n];
                        real nz = _posZ[n];

                        r *= 0.5;

                        real x = nx + (px < nx ? -r : r);
                        real y = ny + (py < ny ? -r : r);
                        real z = nz + (pz < nz ? -r : r);

                        _posX[cell] = x;
                        _posY[cell] = y;
                        _posZ[cell] = z;

                        #pragma unroll NSUB
                        for (int k = 0; k < NSUB; ++k)
                        {
                            _child[NSUB * cell + k] = -1;
                        }

                        if (patch != cell)
                        {
                            _child[NSUB * n + j] = cell;
                        }

                        j = 0;
                        if (x <= _posX[ch])
                            j = 1;
                        if (y <= _posY[ch])
                            j += 2;
                        if (z <= _posZ[ch])
                            j += 4;

                        _child[NSUB * cell + j] = ch;

                        /* The AMD compiler reorders the next read
                         * from _child, which then reads the old/wrong
                         * value when the children are the same without this.
                         */
                        mem_fence(CLK_GLOBAL_MEM_FENCE);

                        n = cell;
                        j = 0;
                        if (x <= px)
                            j = 1;
                        if (y <= py)
                            j += 2;
                        if (z <= pz)
                            j += 4;

                        ch = _child[NSUB * n + j];
                        /* Repeat until the two bodies are different children */
                    }
                    while (ch >= 0);

                    _child[NSUB * n + j] = i;
                    mem_fence(CLK_GLOBAL_MEM_FENCE);
                    _child[locked] = patch;
                }
                mem_fence(CLK_GLOBAL_MEM_FENCE);
                localMaxDepth = max(depth, localMaxDepth);
                i += inc;  /* move on to next body */
                skip = 1;
            }
        }
    }

    atom_max(&_treeStatus->maxDepth, localMaxDepth);
}

/* Used by sw93 */
inline real bmax2Inc(real cmPos, real pPos, real psize)
{
    real dmin = cmPos - (pPos - 0.5 * psize);         /* dist from 1st corner */
    real tmp = max(dmin, psize - dmin);
    return tmp * tmp;      /* sum max distance^2 */
}

inline bool checkTreeDim(real cmPos, real pPos, real halfPsize)
{
    return (cmPos < pPos - halfPsize || cmPos > pPos + halfPsize);
}

__attribute__ ((reqd_work_group_size(THREADS3, 1, 1)))
__kernel void NBODY_KERNEL(summarization)
{
    __local int bottom;
    __local volatile int child[NSUB * THREADS3];

    if (get_local_id(0) == 0)
    {
        bottom = _treeStatus->bottom;
    }
    barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

    int inc = get_local_size(0) * get_num_groups(0);
    int k = (bottom & (-WARPSIZE)) + get_global_id(0);  /* Align to warp size */
    if (k < bottom)
        k += inc;

    int missing = 0;
    while (k <= NNODE) /* Iterate over all cells assigned to thread */
    {
        real m, cm, px, py, pz;
        int cnt, ch;

        if (missing == 0)
        {
            /* New cell, so initialize */
            cm = px = py = pz = 0.0;
            cnt = 0;
            int j = 0;
            #pragma unroll NSUB
            for (int i = 0; i < NSUB; ++i)
            {
                ch = _child[NSUB * k + i];
                if (ch >= 0)
                {
                    if (i != j)
                    {
                        /* Move children to front (needed later for speed) */
                        _child[NSUB * k + i] = -1;
                        _child[NSUB * k + j] = ch;
                    }
                    child[THREADS3 * missing + get_local_id(0)] = ch; /* Cache missing children */
                    m = _mass[ch];
                    ++missing;
                    if (m >= 0.0)
                    {
                        /* Child is ready */
                        --missing;
                        if (ch >= NBODY) /* Count bodies (needed later) */
                        {
                            cnt += _count[ch] - 1;
                        }

                        /* Add child's contribution */
                        cm += m;
                        px += _posX[ch] * m;
                        py += _posY[ch] * m;
                        pz += _posZ[ch] * m;
                    }
                    ++j;
                }
            }
            mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE); /* Only for performance */
            cnt += j;
        }

        if (missing != 0)
        {
            do
            {
                /* poll missing child */
                ch = child[THREADS3 * (missing - 1) + get_local_id(0)];
                m = _mass[ch];
                if (m >= 0.0)
                {
                    /* child is now ready */
                    --missing;
                    if (ch >= NBODY) /* Is a cell */
                    {
                        /* count bodies (needed later) */
                        cnt += _count[ch] - 1;
                    }
                    /* add child's contribution */
                    cm += m;
                    px += _posX[ch] * m;
                    py += _posY[ch] * m;
                    pz += _posZ[ch] * m;
                }
                /* repeat until we are done or child is not ready */
            }
            while ((m >= 0.0) && (missing != 0));
        }

        if (missing == 0)
        {
            /* All children are ready, so store computed information */
            _count[k] = cnt;
            real cx = _posX[k];  /* Load geometric center */
            real cy = _posY[k];
            real cz = _posZ[k];

          #if SW93 || NEWCRITERION
            real psize = _critRadii[k]; /* Get saved size (half cell = radius) */
          #endif

            m = 1.0 / cm;
            px *= m; /* Scale up to position */
            py *= m;
            pz *= m;

            /* Calculate opening criterion if necessary */
          #if SW93
            real bmax2 = bmax2Inc(px, cx, psize);
            bmax2 += bmax2Inc(py, cy, psize);
            bmax2 += bmax2Inc(pz, cz, psize);
            real rc2 = bmax2 / (THETA * THETA);
          #elif NEWCRITERION
            real dx = px - cx;  /* Find distance from center of mass to geometric center */
            real dy = py - cy;
            real dz = pz - cz;
            real dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));

            real rc = (psize / THETA) + dr;

            real rc2 = rc * rc;
          #endif /* SW93 */

          #if SW93 || NEWCRITERION
            /* We don't have the size of the cell for the others, but really still should check */
            bool xTest = checkTreeDim(px, cx, psize);
            bool yTest = checkTreeDim(py, cy, psize);
            bool zTest = checkTreeDim(pz, cz, psize);
            bool structureCheck = xTest || yTest || zTest;
            if (structureCheck)
            {
                _treeStatus->errorCode = NBODY_KERNEL_TREE_STRUCTURE_ERROR;
            }
          #endif /* SW93 || NEWCRITERION */

            _posX[k] = px;
            _posY[k] = py;
            _posZ[k] = pz;

          #if SW93 || NEWCRITERION
            _critRadii[k] = rc2;
          #endif

            write_mem_fence(CLK_GLOBAL_MEM_FENCE); /* Make sure data is visible before setting mass */
            _mass[k] = cm;
            write_mem_fence(CLK_GLOBAL_MEM_FENCE);

            k += inc;  /* Move on to next cell */
        }
    }
}

__attribute__ ((reqd_work_group_size(THREADS4, 1, 1)))
__kernel void NBODY_KERNEL(sort)
{
    __local int bottoms;

    if (get_local_id(0) == 0)
    {
        bottoms = _treeStatus->bottom;
    }
    barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

    int bottom = bottoms;
    int dec = get_local_size(0) * get_num_groups(0);
    int k = NNODE + 1 - dec + get_global_id(0);

    while (k >= bottom) /* Iterate over all cells assigned to thread */
    {
        int start = _start[k];
        if (start >= 0)
        {
            #pragma unroll NSUB
            for (int i = 0; i < NSUB; ++i)
            {
                int ch = _child[NSUB * k + i];
                if (ch >= NBODY)         /* Child is a cell */
                {
                    _start[ch] = start;  /* Set start ID of child */
                    start += _count[ch]; /* Add #bodies in subtree */
                }
                else if (ch >= 0)        /* Child is a body */
                {
                    _sort[start] = ch;  /* Record body in sorted array */
                    ++start;
                }
            }
            k -= dec;  /* Move on to next cell */
        }
    }
}

#if WARPSIZE <= 0
  #error Invalid warp size
#endif

/* These were problems when being lazy and writing it */
#if (THREADS5 / WARPSIZE) <= 0
  #error (THREADS5 / WARPSIZE) must be > 0
#elif (MAXDEPTH * THREADS5 / WARPSIZE) <= 0
  #error (MAXDEPTH * THREADS5 / WARPSIZE) must be > 0
#endif


/* OpenCL is missing thread voting functions.
   This should be equivalent roughtly to CUDA's __all() with the conditions
 * A barrier should be unnecessary here since
 * all the threads in a wavefront should be
 * forced to run simulatenously. This is not
 * over the workgroup, but the actual
 * wavefront.
 * CHECKME: I'm not entirely sure if separate ones needed for each wavefront in a workgroup
 */
inline int forceAllPredicate(__local volatile int allBlock[THREADS1], int warpId, int cond)
{

    allBlock[get_local_id(0)] = cond;

    /* Relies on underlying wavefronts (not whole workgroup)
       executing in lockstep to not require barrier */
    int predicate = 1;
    for (int x = 0; x < WARPSIZE; ++x)
    {
        predicate &= allBlock[WARPSIZE * warpId + x];
    }

    /* For exact, this could always just return false */

    return predicate;
}


__attribute__ ((reqd_work_group_size(THREADS6, 1, 1)))
__kernel void NBODY_KERNEL(velocityIntegration)
{
    int inc = get_local_size(0) * get_num_groups(0);

    /* Iterate over all bodies assigned to thread */
    for (int i = (int) get_global_id(0); i < NBODY; i += inc)
    {
        real dvx = _accX[i] * (0.5 * TIMESTEP);
        real dvy = _accY[i] * (0.5 * TIMESTEP);
        real dvz = _accZ[i] * (0.5 * TIMESTEP);

        real vhx = _velX[i] + dvx;
        real vhy = _velY[i] + dvy;
        real vhz = _velZ[i] + dvz;

        _posX[i] += vhx * TIMESTEP;
        _posY[i] += vhy * TIMESTEP;
        _posZ[i] += vhz * TIMESTEP;

        mem_fence(CLK_GLOBAL_MEM_FENCE);

        _velX[i] = vhx;
        _velY[i] = vhy;
        _velZ[i] = vhz;
    }
}

__attribute__ ((reqd_work_group_size(THREADS5, 1, 1)))
__kernel void NBODY_KERNEL(forceCalculation)
{
    __local int maxDepth;
    __local real rootCritRadius;
    __local int volatile ch[THREADS5 / WARPSIZE];
    __local int volatile pos[MAXDEPTH * THREADS5 / WARPSIZE], node[MAXDEPTH * THREADS5 / WARPSIZE];
    __local volatile real nx[THREADS5 / WARPSIZE], ny[THREADS5 / WARPSIZE], nz[THREADS5 / WARPSIZE];
    __local volatile real nm[THREADS5 / WARPSIZE];

    __local real dq[MAXDEPTH * THREADS5 / WARPSIZE]; /* Used by BH86 and Exact */

    /* Used by the fake thread voting function.
       We rely on the lockstep behaviour of warps/wavefronts to avoid using a barrier
     */
    __local volatile int allBlock[THREADS5];

    /* Excess threads will "die", however their slots in the
     * fake warp vote are still counted, but not set,
     * resulting in garbage in the last few votes. Make sure
     * that dead threads can't prevent a successful vote.
     *
     * Barrier should not be necessary when used, since
     * communication should happen on wavefront/warp level
     */
    allBlock[get_local_id(0)] = 1;

    if (get_local_id(0) == 0)
    {
        maxDepth = _treeStatus->maxDepth;
        real rootSize = _treeStatus->radius;

      #if SW93 || NEWCRITERION
        rootCritRadius = _critRadii[NNODE];
      #endif


      #if BH86
        real rc = rootSize / THETA;
        /* Precompute values that depend only on tree level */
        dq[0] = rc * rc;
        for (int i = 1; i < maxDepth; ++i)
        {
            dq[i] = 0.25 * dq[i - 1];
        }
      #elif EXACT
        real rc = 2.0 * rootSize;
        /* Just fill dq to simplify things. This shouldn't really ever be used anyway */
        for (int i = 0; i < maxDepth; ++i)
        {
            dq[i] = rc * rc;
        }
      #endif /* BH86 */


        if (maxDepth > MAXDEPTH)
        {
            _treeStatus->errorCode = maxDepth;
        }
    }
    barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

    if (maxDepth <= MAXDEPTH)
    {
        /* Figure out first thread in each warp */
        int base = get_local_id(0) / WARPSIZE;
        int sbase = base * WARPSIZE;
        int j = base * MAXDEPTH;
        int diff = get_local_id(0) - sbase; /* Index in warp */

      #if BH86 || EXACT
        /* Make multiple copies to avoid index calculations later */
        if (diff < MAXDEPTH)
        {
            dq[diff + j] = dq[diff];
        }
        barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
      #endif /* BH86 || EXACT */

        /* iterate over all bodies assigned to thread */
        for (int k = get_global_id(0); k < NBODY; k += get_local_size(0) * get_num_groups(0))
        {
            int i = _sort[k];  /* Get permuted index */

            /* Cache position info */
            real px = _posX[i];
            real py = _posY[i];
            real pz = _posZ[i];

            real ax = 0.0;
            real ay = 0.0;
            real az = 0.0;

            /* Initialize iteration stack, i.e., push root node onto stack */
            int depth = j;
            if (get_local_id(0) == sbase)
            {
                node[j] = NNODE;
                pos[j] = 0;

              #if SW93 || NEWCRITERION
                dq[j] = rootCritRadius;
              #endif
            }
            mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

            bool skipSelf = false;
            while (depth >= j)
            {
                /* Stack is not empty */
                while (pos[depth] < 8)
                {
                    int n;
                    /* Node on top of stack has more children to process */
                    if (get_local_id(0) == sbase)
                    {
                        /* I'm the first thread in the warp */
                        n = _child[NSUB * node[depth] + pos[depth]]; /* Load child pointer */
                        pos[depth]++;
                        ch[base] = n; /* Cache child pointer */
                        if (n >= 0)
                        {
                            /* Cache position and mass */
                            nx[base] = _posX[n];
                            ny[base] = _posY[n];
                            nz[base] = _posZ[n];
                            nm[base] = _mass[n];
                        }
                    }
                    mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

                    /* All threads retrieve cached data */
                    n = ch[base];
                    if (n >= 0)
                    {
                        real dx = nx[base] - px;
                        real dy = ny[base] - py;
                        real dz = nz[base] - pz;
                        real rSq = (dx * dx) + (dy * dy) + (dz * dz); /* Compute distance squared */

                        /* Check if all threads agree that cell is far enough away (or is a body) */
                        if (isBody(n) || forceAllPredicate(allBlock, base, rSq >= dq[depth]))
                        {
                            if (n != i) /* Skip self interaction */
                            {
                                real r = sqrt(rSq + EPS2); /* Compute distance with softening */
                                real ai = nm[base] / (r * r * r);
                                ax += ai * dx;
                                ay += ai * dy;
                                az += ai * dz;
                            }

                            if (isBody(n) && n == i) /* Watch for tree incest */
                            {
                                skipSelf = true;
                            }
                        }
                        else
                        {
                            /* Push cell onto stack */
                            ++depth;
                            if (get_local_id(0) == sbase)
                            {
                                node[depth] = n;
                                pos[depth] = 0;

                              #if SW93 || NEWCRITERION
                                dq[depth] = _critRadii[n];
                              #endif
                            }
                            mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
                        }
                    }
                    else
                    {
                        /* Early out because all remaining children are also zero */
                        depth = max(j, depth - 1);
                    }
                }
                --depth;  /* Done with this level */
            }

            if (step > 0)
            {
                _velX[i] += (ax - _accX[i]) * (0.5 * TIMESTEP);
                _velY[i] += (ay - _accY[i]) * (0.5 * TIMESTEP);
                _velZ[i] += (az - _accZ[i]) * (0.5 * TIMESTEP);
            }


            /* Save computed acceleration */
            _accX[i] = ax;
            _accY[i] = ay;
            _accZ[i] = az;

            if (!skipSelf)
            {
                _treeStatus->errorCode = NBODY_KERNEL_TREE_INCEST;
            }


            /* In case this thread is done with bodies and others in the wavefront aren't */
            allBlock[get_local_id(0)] = 1;
        }
    }

}

__attribute__ ((reqd_work_group_size(THREADS6, 1, 1)))
__kernel void NBODY_KERNEL(integration)
{
    int inc = get_local_size(0) * get_num_groups(0);

    /* Iterate over all bodies assigned to thread */
    for (int i = (int) get_global_id(0); i < NBODY; i += inc)
    {
        real dvx = _accX[i] * (0.5 * TIMESTEP);
        real dvy = _accY[i] * (0.5 * TIMESTEP);
        real dvz = _accZ[i] * (0.5 * TIMESTEP);

        real vhx = _velX[i] + dvx;
        real vhy = _velY[i] + dvy;
        real vhz = _velZ[i] + dvz;

        _posX[i] += vhx * TIMESTEP;
        _posY[i] += vhy * TIMESTEP;
        _posZ[i] += vhz * TIMESTEP;

        _velX[i] = vhx + dvx;
        _velY[i] = vhy + dvy;
        _velZ[i] = vhz + dvz;
    }
}

