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

#if DOUBLEPREC
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable


#if DOUBLEPREC
typedef double real;
#else
typedef float real;
#endif /* DOUBLEPREC */

/* FIXME: This needs to be the same as on the host */
typedef struct __attribute__((aligned))
{
    real radius;
    int bottom;
    int maxDepth;
    int errorCode;
    unsigned int blkCnt;
} TreeStatus;

#define NSUB 8

/* All kernels will use the same parameters for now */

#define NBODY_KERNEL(name) name(                                                                 \
    __global real* restrict _posX, __global real* restrict _posY, __global real* restrict _posZ, \
    __global real* restrict _velX, __global real* restrict _velY, __global real* restrict _velZ, \
    __global real* restrict _accX, __global real* restrict _accY, __global real* restrict _accZ, \
                                                                                                 \
    __global real* restrict _maxX, __global real* restrict _maxY, __global real* restrict _maxZ, \
    __global real* restrict _minX, __global real* restrict _minY, __global real* restrict _minZ, \
                                                                                                 \
    __global real* restrict _mass,                                                               \
    __global int* restrict _start, __global int* restrict _count, __global int* restrict _child, \
    __global int* restrict _sort,                                                                \
    __global TreeStatus* restrict _treeStatus,                                                   \
                                                                                                 \
    int step                                                                                     \
    )

__kernel void NBODY_KERNEL(boundingBox)
{
    __local real minX[THREADS1], minY[THREADS1], minZ[THREADS1];
    __local real maxX[THREADS1], maxY[THREADS1], maxZ[THREADS1];

    int i = (int) get_local_id(0);
    if (i == 0)
    {
        minX[0] = _posX[0];
        minY[0] = _posY[0];
        minZ[0] = _posZ[0];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    /* initialize with valid data (in case #bodies < #threads) */
    minX[i] = maxX[i] = minX[0];
    minY[i] = maxY[i] = minY[0];
    minZ[i] = maxZ[i] = minZ[0];

    int inc = get_local_size(0) * get_num_groups(0);
    int j = i + get_group_id(0) * get_local_size(0);

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
        barrier(CLK_LOCAL_MEM_FENCE);
        if (i < j)
        {
            minX[i] = min(minX[i], minZ[i + j]);
            minY[i] = min(minY[i], minZ[i + j]);
            minZ[i] = min(minZ[i], minZ[i + j]);

            maxX[i] = max(maxX[i], maxX[i + j]);
            maxX[i] = max(maxY[i], maxY[i + j]);
            maxX[i] = max(maxZ[i], maxZ[i + j]);
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

        mem_fence(CLK_GLOBAL_MEM_FENCE); /* CHECKME: equivalent to CUDA's __threadfence()? */

        inc = get_num_groups(0) - 1;

        /* CHECKME: CUDA atomicInc seems to have different behaviour when old >= val */
        if (inc == atom_add(&_treeStatus->blkCnt, inc))
        {
            /* I'm the last block so combine all block results */
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
            real tmp = max(maxX[0] - minX[0], maxY[0] - minY[0]);
            _treeStatus->radius = 0.5 * max(tmp, maxZ[0] - minZ[0]);

            /* Create root node */
            j = NNODE;
            _mass[j] = -1.0;
            _start[j] = 0;
            _posX[j] = 0.5 * (minX[0] + maxX[0]);
            _posY[j] = 0.5 * (minY[0] + maxY[0]);
            _posZ[j] = 0.5 * (minZ[0] + maxZ[0]);

            #pragma unroll NSUB
            for (i = 0; i < NSUB; i++)
                _child[NSUB * j + i] = -1;

            _treeStatus->bottom = j;
        }
    }
}

__kernel void NBODY_KERNEL(buildTree)
{
    __local real radius, rootX, rootY, rootZ;

    int i = (int) get_local_id(0);
    if (i == 0)
    {
        /* Cache root data */
        radius = _treeStatus->radius;
        rootX = _posX[NNODE];
        rootY = _posY[NNODE];
        rootZ = _posZ[NNODE];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    int localMaxDepth = 1;
    int skip = 1;

    int inc = get_local_size(0) * get_num_groups(0);
    i += get_group_id(0) * get_local_size(0);

    while (i < NBODY) /* Iterate over all bodies assigned to thread */
    {
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
            j = 0;

            /* Determine which child to follow */
            if (rootX < px)
                j = 1;
            if (rootY < py)
                j += 2;
            if (rootZ < pz)
                j += 4;
        }

        int ch = _child[NSUB * n + j];

        /* Follow path to leaf cell */
        while (ch >= NBODY)
        {
            n = ch;
            ++depth;
            r *= 0.5;
            j = 0;

            /* Determine which child to follow */
            if (rootX < px)
                j = 1;
            if (rootY < py)
                j += 2;
            if (rootZ < pz)
                j += 4;

            ch = _child[NSUB * n + j];
        }

        if (ch != -2) /* Skip if child pointer is locked and try again later */
        {
            int locked = NSUB * n + j;

            /* Try to lock */
            if (ch == atom_cmpxchg(&_child[locked], ch, -2))
            {
                if (ch == -1)
                {
                    /* If null, just insert the new menu */
                    _child[locked] = i;
                }
                else  /* There is already a body in this position */
                {
                    int patch = -1;

                    /* Create new cell(s) and insert the old and new body */
                    do
                    {
                        ++depth;

                        int cell = atom_sub(&_treeStatus->bottom, 1) - 1;
                        if (cell <= NBODY)
                        {
                            _treeStatus->errorCode = 1;
                            _treeStatus->bottom = NNODE;
                        }
                        patch = max(patch, cell);

                        real x = (j & 1) * r;
                        real y = ((j >> 1) & 1) * r;
                        real z = ((j >> 2) & 1) * r;
                        r *= 0.5;

                        _mass[cell] = -1.0;
                        _start[cell] = -1;

                        x = _posX[cell] = _posX[n] - r + x;
                        y = _posY[cell] = _posY[n] - r + y;
                        z = _posZ[cell] = _posZ[n] - r + z;

                        #pragma unroll NSUB
                        for (int k = 0; k < NSUB; ++k)
                        {
                            _child[NSUB * n + j] = cell;
                        }

                        j = 0;
                        if (x < _posX[ch])
                            j = 1;

                        if (y < _posY[ch])
                            j += 2;

                        if (z < _posZ[ch])
                            j += 4;

                        _child[NSUB * cell + j] = ch;

                        n = cell;
                        j = 0;
                        if (rootX < px)
                            j = 1;
                        if (rootY < py)
                            j += 2;
                        if (rootZ < pz)
                            j += 4;

                        ch = _child[NSUB * n + j];
                        /* Repeat until the two bodies are different children */
                    }
                    while (ch >= 0);

                    _child[NSUB * n + j] = i;
                    mem_fence(CLK_GLOBAL_MEM_FENCE); /* CHECKME: same as __threadfence()? */
                    _child[locked] = patch;
                }
                mem_fence(CLK_GLOBAL_MEM_FENCE);

                localMaxDepth = max(depth, localMaxDepth);
                i += inc; /* Move on to next body */
                skip = 1;
            }
        }

        barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
    }

    atom_max(&_treeStatus->maxDepth, localMaxDepth);
}

__kernel void NBODY_KERNEL(summarization)
{
    __local int bottom;
    __local int child[NSUB * THREADS3];

    int i = (int) get_local_id(0);
    if (i == 0)
    {
        bottom = _treeStatus->bottom;
    }
    barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

    int inc = get_local_size(0) * get_num_groups(0);
    int k = (bottom & (-WARPSIZE)) + get_global_id(0);  /* align to warp size */
    if (k < bottom)
        k += inc;

    int missing = 0;
    while (k <= NNODE) /* Iterate over all cells assigned to thread */
    {
        real cm, px, py, pz;
        real m;
        int j;
        int cnt;

        if (missing == 0)
        {
            /* New cell, so initialize */
            cm = px = py = pz = 0.0;
            cnt = j = 0;

            #pragma unroll NSUB
            for (i = 0; i < NSUB; ++i)
            {
                int ch = _child[NSUB * k + i];
                if (ch >= 0)
                {
                    if (i != j)
                    {
                        /* Move children to front (needed later for speed) */
                        _child[NSUB * k + i] = -1;
                        _child[NSUB * k + j] = ch;
                    }
                    child[missing * THREADS3 + get_local_id(0)] = ch; /* Cache missing children */
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
            cnt += j;
        }

        if (missing != 0)
        {
            do
            {
                /* Poll missing children */
                int ch = child[(missing - 1) * THREADS3 + get_local_id(0)];
                m = _mass[ch];
                if (m >= 0.0)
                {
                    /* Child is now ready */
                    --missing;
                    if (ch >= NBODY)
                    {
                        /* Count bodies (needed later) */
                        cnt += _count[ch] - 1;
                    }

                    /* Add child's contribution */
                    cm += m;
                    px += _posX[ch] * m;
                    py += _posY[ch] * m;
                    pz += _posZ[ch] * m;
                }
                /* Repeat until we are done or child is not ready */
            }
            while ((m >= 0.0) && (missing != 0));
        }

        if (missing == 0)
        {
            /* All children are ready, so store computed information */
            _count[k] = cnt;
            m = 1.0 / cm;
            _posX[k] = m * px;
            _posY[k] = m * py;
            _posZ[k] = m * pz;
            mem_fence(CLK_GLOBAL_MEM_FENCE); /* CHECKME: equivalent to CUDA's __threadfence()? */
            _mass[k] = cm;
            k += inc;
        }
    }
}

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

    /* Iterate over all cells assigned to thread */
    while (k >= bottom)
    {
        int start = _start[k];
        if (start >= 0)
        {
            #pragma unroll NSUB
            for (int i = 0; i < NSUB; ++i)
            {
                int ch = _child[NSUB * k + i];
                if (ch >= NBODY)
                {
                    /* Child is a cell */
                    _start[ch] = start; /* Set start ID of child */
                    start += _count[ch]; /* Add #bodies in subtree */
                }
                else if (ch >= 0)
                {
                    /* Child is a body */
                    _sort[start] = ch; /* Record body in sorted array */
                    ++start;
                }
            }
            k -= dec; /* Next cell */
        }
        barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
    }
}

__kernel void NBODY_KERNEL(forceCalculation)
{
    __local int maxDepth;
    __local int ch[THREADS5 / WARPSIZE];
    __local int pos[MAXDEPTH * THREADS5 / WARPSIZE];
    __local int node[MAXDEPTH * THREADS5 / WARPSIZE];
    __local real dq[MAXDEPTH * THREADS5 / WARPSIZE];
    __local real nx[THREADS5 / WARPSIZE], ny[THREADS5 / WARPSIZE], nz[THREADS5 / WARPSIZE];
    __local real nm[THREADS5 / WARPSIZE];

    if (get_local_id(0) == 0)
    {
        /* Precompute values that only depend on tree level */
        maxDepth = _treeStatus->maxDepth;
        real rootSize = _treeStatus->radius;

        real rc;
      #if BH86
        rc = rootSize / THETA;
      #elif SW93
        #error SW93 unimplemented
      #elif NEWCRITERION
        #error NewCriterion unimplemented
      #elif EXACT
        rc = 2.0 * rootSize;
      #else
        #error No opening criterion defined
      #endif /* BH86 */

        dq[0] = rc * rc;
        for (int i = 1; i < maxDepth; ++i)
        {
            dq[i] = 0.25 * dq[i - 1];
        }

        if (maxDepth > MAXDEPTH)
        {
            _treeStatus->errorCode = maxDepth;
        }
    }
    barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

    if (maxDepth <= MAXDEPTH)
    {
        int base = get_local_id(0) / WARPSIZE;
        int sbase = base * WARPSIZE;
        int j = base * MAXDEPTH;

        int diff = get_local_id(0) - sbase;

        /* Make multiple copies to avoid index calculations later */
        if (diff < MAXDEPTH)
        {
            dq[diff + j] = dq[diff];
        }
        barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

        /* Iterate ove all bodies assigned to thread */
        for (int k = get_global_id(0); k < NBODY; k += get_local_size(0) * get_num_groups(0))
        {
            int i = _sort[k]; /* Get permuted index */

            /* Cache position info */
            real px = _posX[i];
            real py = _posY[i];
            real pz = _posZ[i];

            real ax = 0.0;
            real ay = 0.0;
            real az = 0.0;

            /* Initialize iteration stack, i.e. push root node onto stack */
            int depth = j;

            if (sbase == get_local_id(0))
            {
                node[j] = NNODE;
                pos[j] = 0;
            }
            mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

            while (depth >= j)
            {
                /* Stack is not empty */
                while (pos[depth] < 8)
                {
                    int n;
                    /* Node on top of stack has more children to process */
                    if (sbase == get_local_id(0))
                    {
                        /* I'm the first thread in the warp */
                        n = _child[NSUB * node[depth] + pos[depth]]; /* Load child pointer */
                        pos[depth]++;
                        ch[base] = n; /* Cache child pointer */
                        if (n >= 0)
                        {
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

                        /* Compute distance squared */
                        real tmp = (dx * dx) + (dy * dy) + (dz * dz);

                        /* check if all threads agree that cell is far enough away (or is a body) */
                        /* if ((n < nbodiesd) || __all(tmp >= dq[depth])) */
                        /* CHECKME: What exactly is the purpose of the __all here?
                           OpenCL is missing thread voting functions
                         */
                        if (n < NBODY)
                        {
                            if (n != i)
                            {
                                tmp = sqrt(tmp + EPS2);
                                tmp = tmp * tmp * tmp;

                                tmp = nm[base] / tmp;
                                ax = tmp * dx;
                                ay = tmp * dy;
                                az = tmp * dz;
                            }
                        }
                        else
                        {
                            /* Push cell onto stack */
                            ++depth;
                            if (sbase == get_local_id(0))
                            {
                                node[depth] = n;
                                pos[depth] = 0;
                            }
                            mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
                        }
                    }
                    else
                    {
                        /* Early out because all remaining children are also zero */
                        depth = max(j, depth - 1);
                    }
                    --depth; /* Done with this level */
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
            }
        }
    }
}

__kernel void NBODY_KERNEL(integration)
{
    int inc = get_local_size(0) * get_num_groups(0);

    /* Iterate over all bodies assigned to thread */
    for (int i = get_global_id(0); i < NBODY; i += inc)
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

