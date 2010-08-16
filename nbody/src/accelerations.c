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

#include <stdlib.h>
#include "nbody_priv.h"
#include "accelerations.h"

/* CHECKME: order of operations and effect on precision, and where can
 * we share divisions and such */


/* Pure functions are the best ones */

void sphericalAccel(vectorptr restrict acc, const Spherical* sph, const vectorptr restrict pos)
{
    real r;
    ABSV(r, pos);
    const real tmp = sph->scale + r;
    MULVS(acc, pos, -sph->mass / (r * sqr(tmp)));
}

/* gets negative of the acceleration vector of this disk component */
void miyamotoNagaiDiskAccel(vectorptr restrict acc, const Disk* disk, const vectorptr restrict pos)
{
    const real a   = disk->scale_length;
    const real b   = disk->scale_height;
    const real zp  = rsqrt( sqr(Z(pos)) + sqr(b) );
    const real azp = a + zp;

    const real rp  = sqr(X(pos)) + sqr(Y(pos)) + sqr(azp);
    const real rth = rsqrt(cube(rp));  /* rp ^ (3/2) */

    X(acc) = -disk->mass * X(pos) / rth;
    Y(acc) = -disk->mass * Y(pos) / rth;
    Z(acc) = -disk->mass * Z(pos) * azp / (zp * rth);
}

void exponentialDiskAccel(vectorptr restrict acc, const Disk* disk, const vectorptr restrict pos)
{
    const real b = disk->scale_length;
    real r;
    ABSV(r, pos);

    const real expPiece = rexp(-r / b) * (r + b) / b;
    const real factor   = disk->mass * (expPiece - 1.0) / cube(r);
    MULVS(acc, pos, factor);
}

void logHaloAccel(vectorptr restrict acc, const Halo* halo, const vectorptr restrict pos)
{
    const real tvsqr = -2.0 * sqr(halo->vhalo);
    const real qsqr  = sqr(halo->flattenZ);
    const real d     = halo->scale_length;
    const real zsqr  = sqr(Z(pos));

    const real arst  = sqr(d) + sqr(X(pos)) + sqr(Y(pos));
    const real denom = (zsqr / qsqr) +  arst;

    X(acc) = tvsqr * X(pos) / denom;
    Y(acc) = tvsqr * Y(pos) / denom;
    Z(acc) = tvsqr * Z(pos) / ((qsqr * arst) + zsqr);
}

void nfwHaloAccel(vectorptr restrict acc, const Halo* halo, const vectorptr restrict pos)
{
    real r;
    ABSV(r, pos);
    const real a  = halo->scale_length;
    const real ar = a + r;
    const real c  = a * sqr(halo->vhalo) * ((-ar * rlog1p(r / a)) + r) / (0.2162165954 * cube(r) * ar);

    MULVS(acc, pos, c);
}

/* CHECKME: Seems to have precision related issues for a small number of cases for very small qy */
void triaxialHaloAccel(vectorptr restrict acc, const Halo* h, const vectorptr restrict pos)
{
    /* TODO: More things here can be cached */
    const real qzs      = sqr(h->flattenZ);
    const real rhalosqr = sqr(h->scale_length);
    const real mvsqr    = -sqr(h->vhalo);

    const real xsqr = sqr(X(pos));
    const real ysqr = sqr(Y(pos));
    const real zsqr = sqr(Z(pos));

    const real arst  = rhalosqr + (h->c1 * xsqr) + (h->c3 * X(pos) * Y(pos)) + (h->c2 * ysqr);
    const real arst2 = (zsqr / qzs) + arst;

    X(acc) = mvsqr * (((2.0 * h->c1) * X(pos)) + (h->c3 * Y(pos)) ) / arst2;

    Y(acc) = mvsqr * (((2.0 * h->c2) * Y(pos)) + (h->c3 * X(pos)) ) / arst2;

    Z(acc) = (2.0 * mvsqr * Z(pos)) / ((qzs * arst) + zsqr);
}

