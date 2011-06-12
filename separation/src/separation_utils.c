/*
Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
and Rensselaer Polytechnic Institute.

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "separation.h"
#include "milkyway_math.h"
#include "separation_utils.h"
#include "milkyway_util.h"


/* Determine the rotation matrix needed to transform f into t.  The result is an
   array of 9 elements representing a (flattened) 3x3 matrix.
   Adapted from information at http://www.flipcode.com/documents/matrfaq.html */
void get_transform(mwmatrix mat, const mwvector f, const mwvector t)
{
    real angle, sin_a;
    real x, y, z, w;
    mwvector axis;

    axis = mw_crossv(f, t);
    mw_normalize(axis);

    angle = mw_vecangle(f, t);
    sin_a = mw_sin(0.5 * angle);

    x = X(axis) * sin_a;
    y = Y(axis) * sin_a;
    z = Z(axis) * sin_a;
    w = mw_cos(0.5 * angle);

    X(mat[0]) = 1.0 - 2.0 * (y * y + z * z);
    Y(mat[0]) =       2.0 * (x * y - z * w);
    Z(mat[0]) =       2.0 * (x * z + y * w);

    X(mat[1]) =       2.0 * (x * y + z * w);
    Y(mat[1]) = 1.0 - 2.0 * (x * x + z * z);
    Z(mat[1]) =       2.0 * (y * z - x * w);

    X(mat[2]) =       2.0 * (x * z - y * w);
    Y(mat[2]) =       2.0 * (y * z + x * w);
    Z(mat[2]) = 1.0 - 2.0 * (x * x + y * y);
}


/* Transform v by applying the rotation matrix mat */
/* apply coordinate transformations to the given point */
mwvector transform_point(const AstronomyParameters* ap,
                         mwvector point,
                         const mwmatrix cmat,
                         mwvector xsun)
{
    real mcutoff = 11.0;
    mwvector logPoint = xyz_mag(ap, point, mcutoff);

    mw_incsubv(logPoint, xsun);

    return mw_mulmv(cmat, logPoint); /* do transform */
}

static dsfmt_t dsfmtState;

/* Initialize seed for prob_ok; time based seed */
void prob_ok_init(uint32_t seed, int setSeed)
{
    seed = setSeed ? seed : (uint32_t) time(NULL);

    dsfmt_init_gen_rand(&dsfmtState, seed);
}

/* FIXME: WTF? */
/* FIXME: lack of else leads to possibility of returned garbage */
/* determines if star with prob p should be separated into stream */
int prob_ok(StreamStats* ss, int n)
{
    int s_ok = 0;
    real r;
    real step1, step2, step3;

    r = mwXrandom(&dsfmtState, 0.0, 1.0);

    switch (n)
    {
        case 1:
            if (r > ss[0].sprob)
                s_ok = 0;
            else
                s_ok = 1;
            break;
        case 2:
            step1 = ss[0].sprob + ss[1].sprob;
            if (r > step1)
                s_ok = 0;
            else if (r < ss[0].sprob)
                s_ok = 1;
            else if (r > ss[0].sprob && r <= step1)
                s_ok = 2;
            break;
        case 3:
            step1 = ss[0].sprob + ss[1].sprob;
            step2 = ss[0].sprob + ss[1].sprob + ss[2].sprob;
            if (r > step2)
                s_ok = 0;
            else if (r < ss[0].sprob)
                s_ok = 1;
            else if (r > ss[0].sprob && r <= step1)
                s_ok = 2;
            else if (r > step1 && r <= step2)
                s_ok = 3;
            /* CHECKME: else? */
            break;
        case 4:
            step1 = ss[0].sprob + ss[1].sprob;
            step2 = ss[0].sprob + ss[1].sprob + ss[2].sprob;
            step3 = ss[0].sprob + ss[1].sprob + ss[2].sprob + ss[3].sprob;
            if (r > step3)
                s_ok = 0;
            else if (r <= ss[0].sprob)
                s_ok = 1;
            else if (r > ss[0].sprob && r <= step1)
                s_ok = 2;
            else if (r > step1 && r <= step2)
                s_ok = 3;
            else if (r > step2 && r <= step3)
                s_ok = 4;
            break;
        default:
            mw_panic("Too many streams (%d > %d) to separate using current code\n", n, 4);
    }

    return s_ok;
}

SeparationResults* newSeparationResults(unsigned int numberStreams)
{
    unsigned int i;
    SeparationResults* p;

    p = mwCalloc(1, sizeof(SeparationResults));
    p->streamIntegrals = mwCalloc(numberStreams, sizeof(real));
    p->streamLikelihoods = mwCalloc(numberStreams, sizeof(real));

    p->backgroundIntegral = NAN;
    p->backgroundLikelihood = NAN;
    p->likelihood = NAN;

    for (i = 0; i < numberStreams; ++i)
    {
        p->streamIntegrals[i] = NAN;
        p->streamLikelihoods[i] = NAN;
    }

    return p;
}

void freeSeparationResults(SeparationResults* p)
{
    free(p->streamIntegrals);
    free(p->streamLikelihoods);
    free(p);
}

int checkSeparationResults(const SeparationResults* results, unsigned int numberStreams)
{
    unsigned int i;
    int rc = 0;

    rc |= !isfinite(results->likelihood);
    rc |= !isfinite(results->backgroundIntegral);
    rc |= !isfinite(results->backgroundLikelihood);

    for (i = 0; i < numberStreams; ++i)
    {
        rc |= !isfinite(results->streamIntegrals[i]);
        rc |= !isfinite(results->streamLikelihoods[i]);
    }

    if (rc)
        warn("Non-finite result\n");

    return rc;
}

