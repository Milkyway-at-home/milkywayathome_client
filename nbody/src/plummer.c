/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
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

#include "nbody_priv.h"
#include "nbody_util.h"
#include "dSFMT.h"

/* pickshell: pick a random point on a sphere of specified radius. */
inline static void pickshell(dsfmt_t* dsfmtState, vector vec, real rad)
{
    unsigned int k;
    real rsq, rsc;

    do                      /* pick point in NDIM-space */
    {
        for (k = 0; k < NDIM; k++)      /* loop over dimensions */
            vec[k] = xrandom(dsfmtState, -1.0, 1.0);        /* pick from unit cube */
        SQRV(rsq, vec);            /* compute radius squared */
    }
    while (rsq > 1.0);              /* reject if outside sphere */

    rsc = rad / mw_sqrt(rsq);         /* compute scaling factor */
    INCMULVS(vec, rsc);             /* rescale to radius given */
}

#if BOINC_APPLICATION

static void printPlummer(vectorptr restrict rshift, vectorptr restrict vshift)
{
    fprintf(DEFAULT_OUTPUT_FILE,
            "<plummer_r> %.14g %.14g %.14g </plummer_r>\n"
            "<plummer_v> %.14g %.14g %.14g </plummer_v>\n",
            X(rshift), Y(rshift), Z(rshift),
            X(vshift), Y(vshift), Z(vshift));
}

#else

static void printPlummer(vectorptr restrict rshift, vectorptr restrict vshift)
{
    printf("Shifting plummer sphere to r = (%.10f, %.10f, %.10f) v = (%.10f, %.10f, %.10f)\n",
           X(rshift), Y(rshift), Z(rshift),
           X(vshift), Y(vshift), Z(vshift));
}

#endif /* BOINC_APPLICATION */

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

    dsfmt_t dsfmtState;
    real rnd;

    const real rnbody = (real) ctx->model.nbody;
    const real mass   = ctx->model.mass;
    const real mpp    = mass / rnbody;     /* mass per particle */

    // The coordinates to shift the plummer sphere by
    vector rshift = { X(ic->position), Y(ic->position), Z(ic->position) };
    vector vshift = { X(ic->velocity), Y(ic->velocity), Z(ic->velocity) };

    dsfmt_init_gen_rand(&dsfmtState, ctx->seed);

    printPlummer(rshift, vshift);

    rsc = ctx->model.scale_radius;              /* set length scale factor */
    vsc = mw_sqrt(ctx->model.mass / rsc);         /* and recip. speed scale */

    MULVS(scaledrshift, rshift, rsc);   /* Multiply shift by scale factor */
    MULVS(scaledvshift, vshift, vsc);   /* Multiply shift by scale factor */

    endp = st->bodytab + ctx->model.nbody;
    for (p = st->bodytab; p < endp; ++p)   /* loop over particles */
    {
        Type(p) = BODY;    /* tag as a body */
        Mass(p) = mpp;     /* set masses equal */

        /* returns [0, 1) */
        rnd = (real) dsfmt_genrand_close_open(&dsfmtState);

        /* pick r in struct units */

        r = 1.0 / mw_sqrt(mw_pow(rnd, -2.0 / 3.0) - 1.0);

        pickshell(&dsfmtState, Pos(p), rsc * r);     /* pick scaled position */
        INCADDV(Pos(p), rshift);        /* move the position */
        INCADDV(cmr, Pos(p));           /* add to running sum */

        do                      /* select from fn g(x) */
        {
            x = xrandom(&dsfmtState, 0.0, 1.0);      /* for x in range 0:1 */
            y = xrandom(&dsfmtState, 0.0, 0.1);      /* max of g(x) is 0.092 */
        }   /* using von Neumann tech */
        while (y > -cube(x - 1.0) * sqr(x) * cube(x + 1.0) * mw_sqrt(1.0 - sqr(x)));

        v = M_SQRT2 * x / mw_sqrt(mw_sqrt(1.0 + sqr(r)));   /* find v in struct units */
        pickshell(&dsfmtState, Vel(p), vsc * v);        /* pick scaled velocity */
        INCADDV(Vel(p), vshift);       /* move the velocity */
        INCADDV(cmv, Vel(p));         /* add to running sum */
    }

    INCDIVVS(cmr, rnbody);      /* normalize cm coords */
    INCDIVVS(cmv, rnbody);
}


