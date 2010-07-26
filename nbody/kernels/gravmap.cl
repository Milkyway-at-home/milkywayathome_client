/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
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

#include "nbody_types.h"

/* FIXME: I don't belong here */
#define sqr(x) ((x) * (x))
#define cube(x) ((x) * (x) * (x))
#define ZVEC ((vector) 0.0)

#ifndef NULL
  #define NULL ((void*) 0)
#endif


#if SPHERICALTYPE == _SPHERICAL

inline vector sphericalAccel(const Spherical* sph, const vector pos)
{
    real r = length(pos);
    return -sph->mass * pos / (r * sqr(sph->scale + r));
}

#else
  #error "No spherical potential function chosen"
#endif /* SPHERICAL */

#if DISKTYPE == _MN_DISK

inline vector diskAccel(const Disk* disk, const vector pos)
{
    const real a   = disk->scale_length;
    const real b   = disk->scale_height;
    const real zp  = rsqrt( sqr(pos.z) + sqr(b) );
    const real azp = a + zp;
    const real tmp = sqr(pos.x) + sqr(pos.y) + sqr(azp);
    const real rth = sqrt(cube(tmp));

    vector acc = ZVEC;

    acc.x = pos.x / rth;
    acc.y = pos.y / rth;
    acc.z = pos.z * azp / (zp * rth);

    return -disk->mass * acc;
}

#elif DISKTYPE == _EXPONENTIAL_DISK

inline vector diskAccel(const Disk* disk, const vector pos)
{
    const real b        = disk->scale_length;
    const real r        = length(pos);
    const real expPiece = exp(-r / b) * (r + b) / b;
    const real factor   = disk->mass * (expPiece - 1) / cube(r);
    return pos * factor;
}

#else
  #error "No disk acceleration function chosen"
#endif /* MIYAMOTO_NAGAI_DISK */

#if HALOTYPE == _LOG_HALO
inline vector haloAccel(const Halo* halo, const vector pos)
{
    const real tvsqr = -2.0 * sqr(halo->vhalo);
    const real qsqr  = sqr(halo->flattenZ);
    const real d     = halo->scale_length;
    const real zsqr  = sqr(pos.z);

    const real arst  = sqr(d) + sqr(pos.x) + sqr(pos.y);
    const real denom = arst + zsqr / qsqr;

    vector acc = ZVEC;

    acc.x = pos.x / denom;
    acc.y = pos.y / denom;
    acc.z = pos.z / (qsqr * arst + zsqr);

    return tvsqr * acc;
}

#elif HALOTYPE == _NFW_HALO
inline vector nfwHaloAccel(const Halo* halo, const vector pos)
{
    const real r  = length(pos);
    const real a  = halo->scale_length;
    const real ar = a + r;
    const real c  = a * sqr(halo->vhalo) * (r - ar * rlog1p(r / a)) / (0.216 * cube(r) * ar);
    return c * pos;
}

#elif HALOTYPE == _TRIAXIAL_HALO

/* CHECKME: Seems to have precision related issues for a small number of cases for very small qy */
inline vector haloAccel(const Halo* h, const vector pos)
{
    /* TODO: More things here can be cached */
    const real qzs      = sqr(h->flattenZ);
    const real rhalosqr = sqr(h->scale_length);
    const real vsqr     = -sqr(h->vhalo);
    const vector possqr = dot(pos, pos);

    const real arst  = rhalosqr + (h->c1 * possqr.x) + (h->c3 * pos.x * pos.y) + (h->c2 * possqr.y);
    const real arst2 = arst + possqr.z / qzs;

    vector acc = ZVEC;

    acc.x = ( 2 * h->c1 * pos.x + h->c3 * pos.y ) / arst2;
    acc.y = ( 2 * h->c2 * pos.y + h->c3 * pos.x ) / arst2;
    acc.z = 2 * pos.z / (qzs * arst + possqr.z);
    acc *= vsqr;
}

#else
  #error "No halo acceleration chosen"
#endif /* HALOTYPE == _LOG_HALO */

inline vector externalAcc(const NBodyCtx* ctx, const vector r)
{
    return   sphericalAccel(&ctx->pot.sphere[0], r)
           + haloAccel(&ctx->pot.halo, r)
           + diskAccel(&ctx->pot.disk, r);
}

typedef struct
{
    vector pos0;       /* point to evaluate field */
    vector acc0;       /* resulting acceleration */
    vector dr;         /* vector from q to pos0 */
    bodyptr pskip;     /* skip in force evaluation */
    cellptr qmem;      /* data shared with gravsub */
    real drsq;         /* squared distance to pos0 */
} ForceEvalState;

#define EMPTY_FORCE_EVAL_STATE { ZVEC, ZVEC, ZVEC, NULL, NULL, 0.0 }

/* subdivp: decide if cell q is too close to accept as a single
 * term. Also sets qmem, dr, and drsq for use by gravsub.
 */
inline static bool subdivp(ForceEvalState* fest, cellptr q)
{
    fest->dr = Pos(q) - fest->pos0;         /* compute displacement */
    fest->drsq = dot(fest->dr, fest->dr);   /* and find dist squared */
    fest->qmem = q;                         /* remember we know them */
    return (fest->drsq < Rcrit2(q));        /* apply standard rule */
}

/* gravsub: compute contribution of node q to gravitational field at
 * point pos0, and add to running totals phi0 and acc0.
 */
inline static void gravsub(const NBodyCtx* ctx, ForceEvalState* fest, nodeptr q)
{
    real drab, phii, mor3;
    real4 ai, quaddr;
    real dr5inv, phiquad, drquaddr;

    if (q != (nodeptr) fest->qmem)                    /* cant use memorized data? */
    {
        fest->dr = Pos(q) - fest->pos0;               /* then compute sep. */
        fest->drsq = dot(fest->dr, fest->dr);         /* and sep. squared */
    }
    fest->drsq += sqr(ctx->model.eps);                /* use standard softening */
    drab = sqrt(fest->drsq);
    phii = Mass(q) / drab;
    mor3 = phii / fest->drsq;
    ai = fest->dr * mor3;
    fest->acc0 += ai;                     /* ... and to total accel. */
#if 0 /* TODO: matrix stuff */
    if (ctx->usequad && Type(q) == CELL)             /* if cell, add quad term */
    {
        dr5inv = 1.0 / (sqr(fest->drsq) * drab); /* form dr^-5 */
        MULMV(quaddr, Quad(q), fest->dr);        /* form Q * dr */
        drquaddr = dot(fest->dr, quaddr);        /* form dr * Q * dr */
        phiquad = -0.5 * dr5inv * drquaddr;      /* get quad. part of phi */
        phiquad = 5.0 * phiquad / fest->drsq;    /* save for acceleration */
        MULVS(ai, fest->dr, phiquad);            /* components of acc. */
        INCSUBV(fest->acc0, ai);                 /* increment */
        INCMULVS(quaddr, dr5inv);
        fest->acc0 -= quaddr;                    /* acceleration */
    }
#endif
}

/* treescan: iterative routine to do force calculation, starting with
 * node q, which is typically the t.root cell. Watches for tree
 * incest.
 */
inline static bool treescan(const NBodyCtx* ctx,
                            ForceEvalState* fest,
                            nodeptr q)
{
    bool skipself = FALSE;

    while (q != NULL)               /* while not at end of scan */
    {
        if (Type(q) == CELL &&                /* is node a cell and... */
            subdivp(fest, (cellptr) q))       /* too close to accept? */
            q = More(q);            /* follow to next level */
        else                    /* else accept this term */
        {
            if (q == (nodeptr) fest->pskip)    /* self-interaction? */
                skipself = TRUE;               /* then just skip it */
            else                               /* not self-interaction */
                gravsub(ctx, fest, q);         /* so compute gravity */
            q = Next(q);            /* follow next link */
        }
    }

    return skipself;
}

/* hackGrav: evaluate gravitational field on body p; checks to be
 * sure self-interaction was handled correctly if intree is true.
 */
__kernel void gravMap(__const NBodyCtx* ctx,
                      __const nodeptr root,
                      __const bodyptr bs,
                      __global vector* a,
                      const unsigned int count)
{
    size_t id = get_global_id(0);

    bodyptr p = &bs[id];

#if 0
    a[id] = (root == NULL) ? 42.0 : 9000.0;
    return;

    /* Causes things to go invalid, so passing the tree doesn't work
     * as expected */
    nodeptr q = More(root);
    a[id] = (Type(q) == CELL) ? 1337.0 : 55.0;
    return;
#endif

    static bool treeincest = FALSE;     /* tree-incest occured */
    bool skipself          = FALSE;     /* self-interaction skipped */

    ForceEvalState fest = EMPTY_FORCE_EVAL_STATE;

    bool intree = Mass(p) > 0.0;

    fest.pskip = p;             /* exclude p from f.c. */
    fest.pos0  = Pos(p);        /* set field point */

    /* watch for tree-incest */
    skipself = treescan(ctx, &fest, root);         /* scan tree from t.root */
    #if 0
    if (intree && !skipself)            /* did tree-incest occur? */
    {
        if (!ctx->allowIncest) /* treat as catastrophic? */
            fail("hackgrav: tree-incest detected\n");
        if (!treeincest)           /* for the first time? */
            warn("\n[hackgrav: tree-incest detected]\n");
        treeincest = TRUE;          /* don't repeat warning */
    }
    #endif

    /* Adding the external potential */
    fest.acc0 += externalAcc(ctx, Pos(p));

    /* TODO: Sharing */
    a[id] = fest.acc0;
}

