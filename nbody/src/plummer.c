/* ************************************************************************** */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "nbody_priv.h"

#define MFRAC  0.999                /* cut off 1-MFRAC of mass */

/* pickshell: pick a random point on a sphere of specified radius. */
static void pickshell(vector vec, real rad)
{
    int k;
    real rsq, rsc;

    do                      /* pick point in NDIM-space */
    {
        for (k = 0; k < NDIM; k++)      /* loop over dimensions */
            vec[k] = xrandom(-1.0, 1.0);        /* pick from unit cube */
        DOTVP(rsq, vec, vec);           /* compute radius squared */
    }
    while (rsq > 1.0);              /* reject if outside sphere */

    rsc = rad / rsqrt(rsq);         /* compute scaling factor */
    INCMULVS(vec, rsc);             /* rescale to radius given */
}

/* generatePlummer: generate Plummer model initial conditions for test
 * runs, scaled to units such that M = -4E = G = 1 (Henon, Hegge,
 * etc).  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
 * 183.
 */
void generatePlummer(const NBodyCtx* ctx, const InitialConditions* ic, NBodyState* st)
{
    bodyptr p, endp;
    real rsc, vsc, r, v, x, y;
    vector scaledrshift = ZERO_VECTOR;
    vector scaledvshift = ZERO_VECTOR;
    vector cmr          = ZERO_VECTOR;
    vector cmv          = ZERO_VECTOR;

    const real rnbody = (real) ctx->model.nbody;
    const real mass   = ctx->model.mass;
    const real mpp    = mass / rnbody;     /* mass per particle */

    // The coordinates to shift the plummer sphere by
    vector rshift = { ic->position[0], ic->position[1], ic->position[2] };
    vector vshift = { ic->velocity[0], ic->velocity[1], ic->velocity[2] };

    printf("Shifting plummer sphere to r = (%f, %f, %f) v = (%f, %f, %f)...\n",
           rshift[0],
           rshift[1],
           rshift[2],
           vshift[0],
           vshift[1],
           vshift[2]);

    rsc = ctx->model.scale_radius;              /* set length scale factor */
    vsc = rsqrt(ctx->model.mass / rsc);         /* and recip. speed scale */

    MULVS(scaledrshift, rshift, rsc);   /* Multiply shift by scale factor */
    MULVS(scaledvshift, vshift, vsc);   /* Multiply shift by scale factor */

    endp = st->bodytab + ctx->model.nbody;
    for (p = st->bodytab; p < endp; ++p)   /* loop over particles */
    {
        Type(p) = BODY;    /* tag as a body */
        Mass(p) = mpp;     /* set masses equal */
        r = 1.0 / rsqrt(rpow(xrandom(0.0, MFRAC), /* pick r in struct units */
                           -2.0 / 3.0) - 1);
        pickshell(Pos(p), rsc * r);     /* pick scaled position */
        INCADDV(Pos(p), rshift);        /* move the position */
        INCADDV(cmr, Pos(p));           /* add to running sum */

        do                      /* select from fn g(x) */
        {
            x = xrandom(0.0, 1.0);      /* for x in range 0:1 */
            y = xrandom(0.0, 0.1);      /* max of g(x) is 0.092 */
        }
        while (y > x * x * rpow(1 - x * x, 3.5)); /* using von Neumann tech */
        v = rsqrt(2.0) * x / rpow(1 + r * r, 0.25); /* find v in struct units */
        pickshell(Vel(p), vsc * v);     /* pick scaled velocity */
        INCADDV(Vel(p), vshift);       /* move the velocity */
        INCADDV(cmv, Vel(p));         /* add to running sum */
    }

    INCDIVVS(cmr, rnbody);      /* normalize cm coords */
    INCDIVVS(cmv, rnbody);
}


