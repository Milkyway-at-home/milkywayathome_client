/*///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//*/

/*
    Partially based on Yossi Rubner code:
    =========================================================================
    emd.c

    Last update: 3/14/98

    An implementation of the Earth Movers Distance.
    Based of the solution for the Transportation problem as described in
    "Introduction to Mathematical Programming" by F. S. Hillier and
    G. J. Lieberman, McGraw-Hill, 1990.

    Copyright (C) 1998 Yossi Rubner
    Computer Science Department, Stanford University
    E-Mail: rubner@cs.stanford.edu   URL: http://vision.stanford.edu/~rubner
    ==========================================================================
*/

/* likelihood Scaling:
 * Copyright (c) 2016-2018 Siddhartha Shelton

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
#include "milkyway_util.h"
#include "nbody_emd.h"
#include "nbody_defaults.h"
#include "nbody_types.h"
#include "nbody_mass.h"

#define MAX_ITERATIONS 2500
#define EMD_INF   ((real)1.0e20)
#define EMD_EPS   ((real)1.0e-5)
#define EMD_INVALID NAN

typedef enum
{
    EMD_DIST_L1,
    EMD_DIST_L2,
    EMD_DIST_C
} EMDDistanceType;

typedef real (*EMDDistanceFunction)(const real* a, const real* b, void* user_param);

/* EMDNode1D is used for lists, representing 1D sparse array */
typedef struct EMDNode1D
{
    real val;
    struct EMDNode1D* next;
} EMDNode1D;

/* EMDNode2D is used for lists, representing 2D sparse matrix */
typedef struct EMDNode2D
{
    real val;
    struct EMDNode2D* next[2];  /* next row & next column */
    int i, j;
} EMDNode2D;


typedef struct EMDState
{
    int ssize;
    int dsize;

    real** cost;
    EMDNode2D* _x;
    EMDNode2D* end_x;
    EMDNode2D* enter_x;
    char** is_x;

    EMDNode2D** rows_x;
    EMDNode2D** cols_x;

    EMDNode1D* u;
    EMDNode1D* v;

    int* idx1;
    int* idx2;

    /* find_loop buffers */
    EMDNode2D** loop;
    char* is_used;

    /* russel buffers */
    real* s;
    real* d;
    real** delta;

    real weight, max_cost;
    char* buffer;
} EMDState;


/* static function declaration */
static size_t emdAllocateStateBuffer(EMDState* state, int size1, int size2, int dims)
{
    size_t bufferSize;

    /* calculate buffer size */
    bufferSize = (size1 + 1) * (size2 + 1) * (sizeof(real) +   /* cost */
                 sizeof(char) +       /* is_x */
                 sizeof(real)) +     /* delta matrix */
                 (size1 + size2 + 2) * (sizeof(EMDNode2D) +   /* _x */
                                        sizeof(EMDNode2D*) +  /* cols_x & rows_x */
                                        sizeof(EMDNode1D) +   /* u & v */
                                        sizeof(real) +      /* s & d */
                                        sizeof(int) + sizeof(EMDNode2D*)) +    /* idx1 & idx2 */
                 (size1 + 1) * (sizeof(real*) + sizeof(char*) +     /* rows pointers for */
                                sizeof(real*)) + 256;               /*  cost, is_x and delta */

    if (bufferSize < dims * 2 * sizeof(real))
    {
        bufferSize = dims * 2 * sizeof(real);
    }


    state->buffer = mwCalloc(bufferSize, sizeof(char));

    return bufferSize;
}

static void emdReleaseEMD(EMDState* state)
{
    free(state->buffer);
}


/****************************************************************************************\
*                                  standard  metrics                                     *
\****************************************************************************************/
static real emdDistL1(const real* x, const real* y, void* user_param)
{
    int i;
    int dims = (int)(size_t)user_param;
    real s = 0.0;

    for (i = 0; i < dims; i++)
    {
        real t = x[i] - y[i];

        s += mw_fabs(t);
    }

    return (real) s;
}

static real emdDistL2(const real* x, const real* y, void* user_param)
{
    int i;
    int dims = (int)(size_t)user_param;
    real s = 0.0;

    for (i = 0; i < dims; i++)
    {
        real t = x[i] - y[i];

        s += t * t;
    }

    return mw_sqrt((real) s);
}

static real emdDistC(const real* x, const real* y, void* user_param)
{
    int i;
    int dims = (int)(size_t)user_param;
    real s = 0.0;

    for (i = 0; i < dims; i++)
    {
        real t = fabs(x[i] - y[i]);

        if (s < t)
        {
            s = t;
        }
    }

    return (real) s;
}

static int emdFindBasicVariables(real** cost, char** is_x,
                                 EMDNode1D* u, EMDNode1D* v, int ssize, int dsize)
{
    int i, j, found;
    int u_cfound, v_cfound;
    EMDNode1D u0_head, u1_head, *cur_u, *prev_u;
    EMDNode1D v0_head, v1_head, *cur_v, *prev_v;

    /* initialize the rows list (u) and the columns list (v) */
    u0_head.next = u;

    for (i = 0; i < ssize; i++)
    {
        u[i].next = u + i + 1;
    }

    u[ssize - 1].next = 0;
    u1_head.next = 0;

    v0_head.next = ssize > 1 ? v + 1 : 0;

    for (i = 1; i < dsize; i++)
    {
        v[i].next = v + i + 1;
    }

    v[dsize - 1].next = 0;
    v1_head.next = 0;

    /* there are ssize+dsize variables but only ssize+dsize-1 independent equations,
    so set v[0]=0 */
    v[0].val = 0;
    v1_head.next = v;
    v1_head.next->next = 0;

    /* loop until all variables are found */
    u_cfound = v_cfound = 0;

    while (u_cfound < ssize || v_cfound < dsize)
    {
        found = 0;

        if (v_cfound < dsize)
        {
            /* loop over all marked columns */
            prev_v = &v1_head;

            for (found |= (cur_v = v1_head.next) != NULL; cur_v != NULL; cur_v = cur_v->next)
            {
                real cur_v_val = cur_v->val;

                j = (int)(cur_v - v);
                /* find the variables in column j */
                prev_u = &u0_head;

                for (cur_u = u0_head.next; cur_u != NULL;)
                {
                    i = (int)(cur_u - u);

                    if (is_x[i][j])
                    {
                        /* compute u[i] */
                        cur_u->val = cost[i][j] - cur_v_val;
                        /* ...and add it to the marked list */
                        prev_u->next = cur_u->next;
                        cur_u->next = u1_head.next;
                        u1_head.next = cur_u;
                        cur_u = prev_u->next;
                    }
                    else
                    {
                        prev_u = cur_u;
                        cur_u = cur_u->next;
                    }
                }

                prev_v->next = cur_v->next;
                v_cfound++;
            }
        }

        if (u_cfound < ssize)
        {
            /* loop over all marked rows */
            prev_u = &u1_head;

            for (found |= (cur_u = u1_head.next) != NULL; cur_u != NULL; cur_u = cur_u->next)
            {
                real cur_u_val = cur_u->val;
                real* _cost;
                char* _is_x;

                i = (int)(cur_u - u);
                _cost = cost[i];
                _is_x = is_x[i];
                /* find the variables in rows i */
                prev_v = &v0_head;

                for (cur_v = v0_head.next; cur_v != NULL;)
                {
                    j = (int)(cur_v - v);

                    if (_is_x[j])
                    {
                        /* compute v[j] */
                        cur_v->val = _cost[j] - cur_u_val;
                        /* ...and add it to the marked list */
                        prev_v->next = cur_v->next;
                        cur_v->next = v1_head.next;
                        v1_head.next = cur_v;
                        cur_v = prev_v->next;
                    }
                    else
                    {
                        prev_v = cur_v;
                        cur_v = cur_v->next;
                    }
                }

                prev_u->next = cur_u->next;
                u_cfound++;
            }
        }

        if (!found)
        {
            return -1;
        }
    }

    return 0;
}

static void emdAddBasicVariable(EMDState* state,
                                int min_i, int min_j,
                                EMDNode1D* prev_u_min_i, EMDNode1D* prev_v_min_j, EMDNode1D* u_head)
{
    real temp;
    EMDNode2D* end_x = state->end_x;

    if (state->s[min_i] < state->d[min_j] + state->weight * EMD_EPS)
    {
        /* supply exhausted */
        temp = state->s[min_i];
        state->s[min_i] = 0;
        state->d[min_j] -= temp;
    }
    else                        /* demand exhausted */
    {
        temp = state->d[min_j];
        state->d[min_j] = 0;
        state->s[min_i] -= temp;
    }

    /* x(min_i,min_j) is a basic variable */
    state->is_x[min_i][min_j] = 1;

    end_x->val = temp;
    end_x->i = min_i;
    end_x->j = min_j;
    end_x->next[0] = state->rows_x[min_i];
    end_x->next[1] = state->cols_x[min_j];
    state->rows_x[min_i] = end_x;
    state->cols_x[min_j] = end_x;
    state->end_x = end_x + 1;

    /* delete supply row only if the empty, and if not last row */
    if (state->s[min_i] == 0.0 && u_head->next->next != NULL)
    {
        prev_u_min_i->next = prev_u_min_i->next->next;    /* remove row from list */
    }
    else
    {
        prev_v_min_j->next = prev_v_min_j->next->next;    /* remove column from list */
    }
}

static void emdRussel(EMDState* state)
{
    int i, j, min_i = -1, min_j = -1;
    real min_delta, diff;

    EMDNode1D u_head;
    EMDNode1D* cur_u;
    EMDNode1D* prev_u;

    EMDNode1D v_head;
    EMDNode1D* cur_v;
    EMDNode1D* prev_v;

    EMDNode1D* prev_u_min_i = NULL;
    EMDNode1D* prev_v_min_j = NULL;
    EMDNode1D* remember;

    EMDNode1D* u = state->u;
    EMDNode1D* v = state->v;

    int ssize = state->ssize;
    int dsize = state->dsize;
    real eps = EMD_EPS * state->max_cost;
    real** cost = state->cost;
    real** delta = state->delta;

    /* initialize the rows list (ur), and the columns list (vr) */
    u_head.next = u;

    for (i = 0; i < ssize; i++)
    {
        u[i].next = u + i + 1;
    }

    u[ssize - 1].next = 0;

    v_head.next = v;

    for (i = 0; i < dsize; i++)
    {
        v[i].val = -EMD_INF;
        v[i].next = v + i + 1;
    }

    v[dsize - 1].next = 0;

    /* find the maximum row and column values (ur[i] and vr[j]) */
    for (i = 0; i < ssize; i++)
    {
        real u_val = -EMD_INF;
        real* cost_row = cost[i];

        for (j = 0; j < dsize; j++)
        {
            real temp = cost_row[j];

            if (u_val < temp)
            {
                u_val = temp;
            }

            if (v[j].val < temp)
            {
                v[j].val = temp;
            }
        }

        u[i].val = u_val;
    }

    /* compute the delta matrix */
    for (i = 0; i < ssize; i++)
    {
        real u_val = u[i].val;
        real* delta_row = delta[i];
        real* cost_row = cost[i];

        for (j = 0; j < dsize; j++)
        {
            delta_row[j] = cost_row[j] - u_val - v[j].val;
        }
    }

    /* find the basic variables */
    do
    {
        /* find the smallest delta[i][j] */
        min_i = -1;
        min_delta = EMD_INF;
        prev_u = &u_head;

        for (cur_u = u_head.next; cur_u != NULL; cur_u = cur_u->next)
        {
            real* delta_row;

            i = (int)(cur_u - u);
            delta_row = delta[i];

            prev_v = &v_head;

            for (cur_v = v_head.next; cur_v != NULL; cur_v = cur_v->next)
            {
                j = (int)(cur_v - v);

                if (min_delta > delta_row[j])
                {
                    min_delta = delta_row[j];
                    min_i = i;
                    min_j = j;
                    prev_u_min_i = prev_u;
                    prev_v_min_j = prev_v;
                }

                prev_v = cur_v;
            }

            prev_u = cur_u;
        }

        if (min_i < 0)
        {
            break;
        }

        /* add x[min_i][min_j] to the basis, and adjust supplies and cost */
        remember = prev_u_min_i->next;
        emdAddBasicVariable(state, min_i, min_j, prev_u_min_i, prev_v_min_j, &u_head);

        /* update the necessary delta[][] */
        if (remember == prev_u_min_i->next)    /* line min_i was deleted */
        {
            for (cur_v = v_head.next; cur_v != NULL; cur_v = cur_v->next)
            {
                j = (int)(cur_v - v);

                if (cur_v->val == cost[min_i][j])      /* column j needs updating */
                {
                    real max_val = -EMD_INF;

                    /* find the new maximum value in the column */
                    for (cur_u = u_head.next; cur_u != NULL; cur_u = cur_u->next)
                    {
                        real temp = cost[cur_u - u][j];

                        if (max_val < temp)
                        {
                            max_val = temp;
                        }
                    }

                    /* if needed, adjust the relevant delta[*][j] */
                    diff = max_val - cur_v->val;
                    cur_v->val = max_val;

                    if (mw_fabs(diff) < eps)
                    {
                        for (cur_u = u_head.next; cur_u != NULL; cur_u = cur_u->next)
                        {
                            delta[cur_u - u][j] += diff;
                        }
                    }
                }
            }
        }
        else                    /* column min_j was deleted */
        {
            for (cur_u = u_head.next; cur_u != NULL; cur_u = cur_u->next)
            {
                i = (int)(cur_u - u);

                if (cur_u->val == cost[i][min_j])      /* row i needs updating */
                {
                    real max_val = -EMD_INF;

                    /* find the new maximum value in the row */
                    for (cur_v = v_head.next; cur_v != NULL; cur_v = cur_v->next)
                    {
                        real temp = cost[i][cur_v - v];

                        if (max_val < temp)
                        {
                            max_val = temp;
                        }
                    }

                    /* if needed, adjust the relevant delta[i][*] */
                    diff = max_val - cur_u->val;
                    cur_u->val = max_val;

                    if (mw_fabs(diff) < eps)
                    {
                        for (cur_v = v_head.next; cur_v != NULL; cur_v = cur_v->next)
                        {
                            delta[i][cur_v - v] += diff;
                        }
                    }
                }
            }
        }
    }
    while (u_head.next != NULL || v_head.next != NULL);
}

/************************************************************************************\
*          initialize structure, allocate buffers and generate initial golution      *
\************************************************************************************/
static int emdInitEMD(const real* signature1, int size1,
                      const real* signature2, int size2,
                      int dims, EMDDistanceFunction dist_func, void* user_param,
                      const real* cost, int cost_step,
                      EMDState* state, real* lower_bound)
{
    real s_sum = 0.0, d_sum = 0.0, diff;
    int i, j;
    int ssize = 0;
    int dsize = 0;
    int equal_sums = 1;
    size_t buffer_size;
    real max_cost = 0.0;
    char* buffer;
    char* buffer_end;

    memset(state, 0, sizeof(*state));
    assert(cost_step % sizeof(real) == 0);
    cost_step /= sizeof(real);

    buffer_size = emdAllocateStateBuffer(state, size1, size2, dims);
    buffer = state->buffer;
    buffer_end = buffer + buffer_size;

    state->idx1 = (int*) buffer;
    buffer += (size1 + 1) * sizeof(int);

    state->idx2 = (int*) buffer;
    buffer += (size2 + 1) * sizeof(int);

    state->s = (real*) buffer;
    buffer += (size1 + 1) * sizeof(real);

    state->d = (real*) buffer;
    buffer += (size2 + 1) * sizeof(real);

    /* sum up the supply and demand */
    for (i = 0; i < size1; i++)
    {
        real weight = signature1[i * (dims + 1)];

        if (weight > 0.0)
        {
            s_sum += weight;
            state->s[ssize] = weight;
            state->idx1[ssize++] = i;

        }
        else if (weight < 0.0)
        {
            mw_printf("Weight out of range\n");
            return -1;
        }
    }

    for (i = 0; i < size2; i++)
    {
        real weight = signature2[i * (dims + 1)];

        if (weight > 0.0)
        {
            d_sum += weight;
            state->d[dsize] = weight;
            state->idx2[dsize++] = i;
        }
        else if (weight < 0.0)
        {
            mw_printf("Weight out of range\n");
            return -1;
        }
    }

    if (ssize == 0 || dsize == 0)
    {
        mw_printf("ssize or dsize out of range\n");
        return -1;
    }

    /* if supply different than the demand, add a zero-cost dummy cluster */
    diff = s_sum - d_sum;

    if (mw_fabs(diff) >= EMD_EPS * s_sum)
    {
        equal_sums = 0;

        if (diff < 0)
        {
            state->s[ssize] = -diff;
            state->idx1[ssize++] = -1;
        }
        else
        {
            state->d[dsize] = diff;
            state->idx2[dsize++] = -1;
        }
    }

    state->ssize = ssize;
    state->dsize = dsize;
    state->weight = s_sum > d_sum ? s_sum : d_sum;

    if (lower_bound && equal_sums)     /* check lower bound */
    {
        int sz1 = size1 * (dims + 1);
        int sz2 = size2 * (dims + 1);
        real lb = 0.0;

        real* xs = (real*) buffer;
        real* xd = xs + dims;

        memset(xs, 0, dims * sizeof(xs[0]));
        memset(xd, 0, dims * sizeof(xd[0]));

        for (j = 0; j < sz1; j += dims + 1)
        {
            real weight = signature1[j];

            for (i = 0; i < dims; i++)
            {
                xs[i] += signature1[j + i + 1] * weight;
            }
        }

        for (j = 0; j < sz2; j += dims + 1)
        {
            real weight = signature2[j];

            for (i = 0; i < dims; i++)
            {
                xd[i] += signature2[j + i + 1] * weight;
            }
        }

        lb = dist_func(xs, xd, user_param) / state->weight;
        i = *lower_bound <= lb;
        *lower_bound = lb;

        if (i)
        {
            return 1;
        }
    }

    /* assign pointers */
    state->is_used = (char*) buffer;
    /* init delta matrix */
    state->delta = (real**) buffer;
    buffer += ssize * sizeof(real*);

    for (i = 0; i < ssize; i++)
    {
        state->delta[i] = (real*) buffer;
        buffer += dsize * sizeof(real);
    }

    state->loop = (EMDNode2D**) buffer;
    buffer += (ssize + dsize + 1) * sizeof(EMDNode2D*);

    state->_x = state->end_x = (EMDNode2D*) buffer;
    buffer += (ssize + dsize) * sizeof(EMDNode2D);

    /* init cost matrix */
    state->cost = (real**) buffer;
    buffer += ssize * sizeof(real*);

    /* compute the distance matrix */
    for (i = 0; i < ssize; i++)
    {
        int ci = state->idx1[i];

        state->cost[i] = (real*) buffer;
        buffer += dsize * sizeof(real);

        if (ci >= 0)
        {
            for (j = 0; j < dsize; j++)
            {
                int cj = state->idx2[j];

                if (cj < 0)
                {
                    state->cost[i][j] = 0;
                }
                else
                {
                    real val;

                    if (dist_func)
                    {
                        val = dist_func(signature1 + ci * (dims + 1) + 1,
                                        signature2 + cj * (dims + 1) + 1,
                                        user_param);
                    }
                    else
                    {
                        assert(cost);
                        val = cost[cost_step * ci + cj];
                    }

                    state->cost[i][j] = val;

                    if (max_cost < val)
                    {
                        max_cost = val;
                    }
                }
            }
        }
        else
        {
            for (j = 0; j < dsize; j++)
            {
                state->cost[i][j] = 0;
            }
        }
    }

    state->max_cost = max_cost;

    memset(buffer, 0, buffer_end - buffer);

    state->rows_x = (EMDNode2D**) buffer;
    buffer += ssize * sizeof(EMDNode2D*);

    state->cols_x = (EMDNode2D**) buffer;
    buffer += dsize * sizeof(EMDNode2D*);

    state->u = (EMDNode1D*) buffer;
    buffer += ssize * sizeof(EMDNode1D);

    state->v = (EMDNode1D*) buffer;
    buffer += dsize * sizeof(EMDNode1D);

    /* init is_x matrix */
    state->is_x = (char**) buffer;
    buffer += ssize * sizeof(char*);

    for (i = 0; i < ssize; i++)
    {
        state->is_x[i] = buffer;
        buffer += dsize;
    }

    assert(buffer <= buffer_end);

    emdRussel(state);

    state->enter_x = (state->end_x)++;
    return 0;
}

static real emdIsOptimal(real** cost, char** is_x,
                          EMDNode1D* u, EMDNode1D* v, int ssize, int dsize, EMDNode2D* enter_x)
{
    real delta, min_delta = EMD_INF;
    int i, j;
    int min_i = 0;
    int min_j = 0;

    /* find the minimal cij-ui-vj over all i,j */
    for (i = 0; i < ssize; i++)
    {
        real u_val = u[i].val;
        real* _cost = cost[i];
        char* _is_x = is_x[i];

        for (j = 0; j < dsize; j++)
        {
            if (!_is_x[j])
            {
                delta = _cost[j] - u_val - v[j].val;

                if (min_delta > delta)
                {
                    min_delta = delta;
                    min_i = i;
                    min_j = j;
                }
            }
        }
    }

    enter_x->i = min_i;
    enter_x->j = min_j;

    return min_delta;
}

static int emdFindLoop(EMDState* state)
{
    int i;
    int steps = 1;
    EMDNode2D* new_x;
    EMDNode2D** loop = state->loop;
    EMDNode2D* enter_x = state->enter_x;
    EMDNode2D* _x = state->_x;
    char* is_used = state->is_used;

    memset(is_used, 0, state->ssize + state->dsize);

    new_x = loop[0] = enter_x;
    is_used[enter_x - _x] = 1;
    steps = 1;

    do
    {
        if ((steps & 1) == 1)
        {
            /* find an unused x in the row */
            new_x = state->rows_x[new_x->i];

            while (new_x != NULL && is_used[new_x - _x])
            {
                new_x = new_x->next[0];
            }
        }
        else
        {
            /* find an unused x in the column, or the entering x */
            new_x = state->cols_x[new_x->j];

            while (new_x != NULL && is_used[new_x - _x] && new_x != enter_x)
            {
                new_x = new_x->next[1];
            }

            if (new_x == enter_x)
            {
                break;
            }
        }

        if (new_x != NULL)        /* found the next x */
        {
            /* add x to the loop */
            loop[steps++] = new_x;
            is_used[new_x - _x] = 1;
        }
        else                    /* didn't find the next x */
        {
            /* backtrack */
            do
            {
                i = steps & 1;
                new_x = loop[steps - 1];

                do
                {
                    new_x = new_x->next[i];
                }
                while (new_x != NULL && is_used[new_x - _x]);

                if (new_x == NULL)
                {
                    is_used[loop[--steps] - _x] = 0;
                }
            }
            while (new_x == NULL && steps > 0);

            is_used[loop[steps - 1] - _x] = 0;
            loop[steps - 1] = new_x;
            is_used[new_x - _x] = 1;
        }
    }
    while (steps > 0);

    return steps;
}

static mwbool emdNewSolution(EMDState* state)
{
    int i, j;
    real min_val = EMD_INF;
    int steps;
    EMDNode2D head;
    EMDNode2D* cur_x;
    EMDNode2D* next_x;
    EMDNode2D* leave_x = NULL;
    EMDNode2D* enter_x = state->enter_x;
    EMDNode2D** loop = state->loop;

    /* enter the new basic variable */
    i = enter_x->i;
    j = enter_x->j;
    state->is_x[i][j] = 1;
    enter_x->next[0] = state->rows_x[i];
    enter_x->next[1] = state->cols_x[j];
    enter_x->val = 0;
    state->rows_x[i] = enter_x;
    state->cols_x[j] = enter_x;

    /* find a chain reaction */
    steps = emdFindLoop(state);

    if (steps == 0)
    {
        return FALSE;
    }

    /* find the largest value in the loop */
    for (i = 1; i < steps; i += 2)
    {
        real temp = loop[i]->val;

        if (min_val > temp)
        {
            leave_x = loop[i];
            min_val = temp;
        }
    }

    if (!leave_x)
    {
        return FALSE;
    }

    /* update the loop */
    for (i = 0; i < steps; i += 2)
    {
        real temp0 = loop[i]->val + min_val;
        real temp1 = loop[i + 1]->val - min_val;

        loop[i]->val = temp0;
        loop[i + 1]->val = temp1;
    }

    /* remove the leaving basic variable */
    i = leave_x->i;
    j = leave_x->j;
    state->is_x[i][j] = 0;

    head.next[0] = state->rows_x[i];
    cur_x = &head;

    while ((next_x = cur_x->next[0]) != leave_x)
    {
        cur_x = next_x;
        assert(cur_x);
    }

    cur_x->next[0] = next_x->next[0];
    state->rows_x[i] = head.next[0];

    head.next[1] = state->cols_x[j];
    cur_x = &head;

    while ((next_x = cur_x->next[1]) != leave_x)
    {
        cur_x = next_x;
        assert(cur_x);
    }

    cur_x->next[1] = next_x->next[1];
    state->cols_x[j] = head.next[1];

    /* set enter_x to be the new empty slot */
    state->enter_x = leave_x;

    return TRUE;
}

static int emdIterateSolution(EMDState* state)
{
    int result;
    real min_delta;
    real eps = EMD_EPS * state->max_cost;

    /* if ssize = 1 or dsize = 1 then we are done, else ... */
    if (state->ssize > 1 && state->dsize > 1)
    {
        int itr;

        for (itr = 1; itr < MAX_ITERATIONS; itr++)
        {
            /* find basic variables */
            result = emdFindBasicVariables(state->cost, state->is_x,
                                           state->u, state->v, state->ssize, state->dsize);
            if (result < 0)
            {
                break;
            }

            /* check for optimality */
            min_delta = emdIsOptimal(state->cost, state->is_x,
                                     state->u, state->v,
                                     state->ssize, state->dsize, state->enter_x);

            if (min_delta == EMD_INF)
            {
                mw_printf("Iteration didn't converge");
                return 1;
            }

            /* if no negative deltamin, we found the optimal solution */
            if (min_delta >= -eps)
            {
                break;
            }

            /* improve solution */
            if (!emdNewSolution(state))
            {
                mw_printf("Iteration didn't converge");
                return 1;
            }
        }
    }

    return 0;
}

static void emdPrintFlowMatrix(const real* flow, int size1, int size2)
{
    int i, j;
    const int flowStep = 1;

    if (!flow)
    {
        mw_printf("Empty flow\n");
        return;
    }

    for (i = 0; i < size1; ++i)
    {
        for (j = 0; j < size2; ++j)
        {
            mw_printf("Flow[ % d][ % d] = % f\n", i, j, flow[flowStep * j + i]);
        }
    }
}

static EMDDistanceFunction nbMetricDistanceFunction(EMDDistanceType distType)
{
    switch (distType)
    {
    case EMD_DIST_L1:
        return emdDistL1;
        break;

    case EMD_DIST_L2:
        return emdDistL2;
        break;

    case EMD_DIST_C:
        return emdDistC;
        break;

    default:
        mw_panic("Bad or unsupported metric type");
        return NULL;
    }
}

static real emdComputeTotalFlow(EMDState* state, real* flow)
{
    EMDNode2D* xp = NULL;
    real totalCost = 0.0;
    const int flowStep = 1;

    for (xp = state->_x; xp < state->end_x; xp++)
    {
        int ci, cj;
        real val = xp->val;
        int i = xp->i;
        int j = xp->j;

        if (xp == state->enter_x)
        {
            continue;
        }

        ci = state->idx1[i];
        cj = state->idx2[j];

        if (ci >= 0 && cj >= 0)
        {
            totalCost += (real) val * state->cost[i][j];

            if (flow)
            {
                flow[flowStep * ci + cj] = val;
            }
        }
    }
    return totalCost;
}

/* The main function */
real emdCalc(const real* RESTRICT signature_arr1,
              const real* RESTRICT signature_arr2,
              unsigned int size1,
              unsigned int size2,
              real* RESTRICT lower_bound)
{
    EMDState state;
    real emd = (real) EMD_INVALID;
    real totalCost = 0.0;
    int result = 0;
    EMDDistanceFunction dist_func = NULL;
    const EMDDistanceType dist_type = EMD_DIST_L2;
    const mwbool debugFlow = FALSE;
    real* flow = NULL;
    const int dims = 2; /* We have 2 dimensions, lambda and beta */
    void* user_param = (void*) dims;

    memset(&state, 0, sizeof(state));

    dist_func = nbMetricDistanceFunction(dist_type);
    result = emdInitEMD(signature_arr1, size1,
                        signature_arr2, size2,
                        dims, dist_func, user_param,
                        NULL, 0,
                        &state, lower_bound);

    if (result > 0 && lower_bound)
    {
        emd = *lower_bound;
        emdReleaseEMD(&state);
        return emd;
    }
    else if (result < 0)
    {
        emdReleaseEMD(&state);
        return (real) EMD_INVALID;
    }

    if (debugFlow)
    {
        flow = mwCalloc(size1 * size2, sizeof(real));
    }

    if (!emdIterateSolution(&state))
    {
        //mw_printf("emdComputeTotalFlow ACCESSED!\n");
        totalCost = emdComputeTotalFlow(&state, flow);
        //mw_printf("totalCost = %.15f\n", totalCost);
        //mw_printf("weight = %.15f\n", state.weight);
        emd = (real)(totalCost / state.weight);
        //mw_printf("emd = %.15f\n", emd);
    }

    if (debugFlow)
    {
        emdPrintFlowMatrix(flow, size1, size2);
    }

    free(flow);
    emdReleaseEMD(&state);
    return emd;
}


real nbWorstCaseEMD(const NBodyHistogram* hist)
{
    //(This makes no sense to be defined this way now that histograms are not normalized.
    //  return fabs(hist->data[0].lambda - hist->data[hist->nBin - 1].lambda);
    return DEFAULT_WORST_CASE;
}

real nbMatchEMD(const MainStruct* data, const MainStruct* histogram)
{
    // all the histograms have the same lambda/betaBins info
    // but histogram 0 is the one that contains the counts information
    NBodyHistogram* first_data = (data->histograms[0]);
    NBodyHistogram* first_hist = (histogram->histograms[0]);

    unsigned int lambdaBins = first_data->lambdaBins;
    unsigned int betaBins = first_data->betaBins;
    unsigned int bins = lambdaBins * betaBins;
    unsigned int nSim_uncut = first_hist->totalNum;
    unsigned int nSim = nSim_uncut;
    unsigned int nData = first_data->totalNum;
    real histMass = first_hist->massPerParticle;
    real dataMass = first_data->massPerParticle;
    unsigned int i;
    unsigned int j;
    unsigned int rawCount;
    WeightPos* hist;
    WeightPos* dat;
    real EMDStart;
    real EMDEnd;
    unsigned int rangeCount;
    unsigned int rangeBins;
    unsigned int totalRangeCount; // Total counts in all EMD Ranges, used for weighting the EMD
    real emd;
    real likelihood;

    if (first_data->lambdaBins != first_hist->lambdaBins || first_data->betaBins != first_hist->betaBins)
    {
        /* FIXME?: We could have mismatched histogram sizes, but I'm
        * not sure what to do with ignored bins and
        * renormalization */
        return NAN;
    }

    if (nSim == 0 || nData == 0)
    {
        /* If the histogram is totally empty, it is worse than the worst case */
        return INFINITY;
    }

    if (histMass <= 0.0 || dataMass <= 0.0)
    {
        /*In order to calculate likelihood the masses are necessary*/
        return NAN;
    }

    if (first_data->params.nRange%2 != 0)
    {
        mw_printf("Error reading EMD calculation ranges: odd number of lambda values recieved \n");
        return NAN;
    }

    emd = 0;
    totalRangeCount = 0;
    if(first_data->params.nRange == 0)
    {
        // mw_printf("No EMD Ranges defined, using full histogram\n");
        first_data->params.nRange = 2;
        first_data->params.EMDRange[0] = first_data->data[0].lambda;
        first_data->params.EMDRange[1] = first_data->data[bins - 1].lambda;
    }
    for(i = 0; i < first_data->params.nRange; i = i + 2)
    {
        /*Renormalize simulated hist to given EMD Range*/
        EMDStart = first_data->params.EMDRange[i];
        EMDEnd = first_data->params.EMDRange[i+1];
        // mw_printf("Using EMD Range: {%f,%f}\n", EMDStart, EMDEnd);
        if(EMDStart > EMDEnd)
        {
            mw_printf("Error reading EMD calculation ranges: EMDStart > EMDEnd \n");
            return NAN;
        }
        rangeCount = 0;
        rangeBins = 0;
        for(j = 0; j < bins; j++) /*Sum up total counts of bins in emd range*/
        {
            if(first_hist->data[j].lambda >= EMDStart && first_hist->data[j].lambda <= EMDEnd && first_data->data[j].useBin)
            {
                rangeCount += mw_round(first_hist->data[j].variable * nSim_uncut);
                rangeBins += 1;
            }
        }
        hist = mwCalloc(rangeBins, sizeof(WeightPos)); /*Create histogram emdCalc can use*/
        rangeBins = 0;
        for(j = 0; j < bins; j++)
        {
            if(first_hist->data[j].lambda >= EMDStart && first_hist->data[j].lambda <= EMDEnd && first_data->data[j].useBin)
            {
                hist[rangeBins].weight = (real) first_hist->data[j].variable * nSim_uncut / (1.0*rangeCount);
                hist[rangeBins].lambda = (real) first_hist->data[j].lambda;
                hist[rangeBins].beta = (real) first_hist->data[j].beta;
                rangeBins += 1;
            }
        }

        /*Repeat for the input data hist*/
        rangeCount = 0;
        rangeBins = 0;
        for(j = 0; j < bins; j++) /*Sum up total counts of bins in emd range*/
        {
            if(first_data->data[j].lambda >= EMDStart && first_data->data[j].lambda <= EMDEnd && first_data->data[j].useBin)
            {
                rangeCount += mw_round(first_data->data[j].variable * nData);
                rangeBins += 1;
            }
        }
        dat = mwCalloc(rangeBins, sizeof(WeightPos)); /*Create histogram emdCalc can use*/
        totalRangeCount += rangeCount;
        rangeBins = 0;
        for(j = 0; j < bins; j++)
        {
            if(first_data->data[j].lambda >= EMDStart && first_data->data[j].lambda <= EMDEnd && first_data->data[j].useBin)
            {
                dat[rangeBins].weight = (real) first_data->data[j].variable * nData / (1.0*rangeCount);
                dat[rangeBins].lambda = (real) first_data->data[j].lambda;
                dat[rangeBins].beta = (real) first_data->data[j].beta;
                rangeBins += 1;
            }
        }

        /*Calculate EMD, each range is weighted by the % of counts in that range in data hist*/
        emd += emdCalc((const real*) dat, (const real*) hist, rangeBins, rangeBins, NULL) * (rangeCount/(1.0*nData));  //temporary weighting by total counts
        free(hist);
        free(dat);
    }     

    emd = emd * (nData / totalRangeCount); //correct weighting to only consider counts in EMD ranges
    emd *= 1.0e9;
    emd = mw_round(emd);
    emd *= 1.0e-9;
    
    if (emd > 50.0)
    {
        /* emd's max value is 50 */
        return NAN;
    }

    /* This calculates the likelihood as the combination of the
    * probability distribution and (1.0 - emd / max_dist) */

    real EMDComponent = 1.0 - emd / 50.0;
    
    /* the 300 is there to add weight to the EMD component */
    likelihood = 300.0 * mw_log(EMDComponent);
    
    /* the emd is a negative. returning a positive value */
    return -likelihood;
}

