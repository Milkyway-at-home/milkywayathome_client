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
#include "nbody_potential.h"
#include "milkyway_util.h"

static inline mwvector sphericalAccel(const Spherical* sph, const mwvector pos)
{
    const real r   = mw_absv(pos);
    const real tmp = sph->scale + r;

    return mw_mulvs(pos, -sph->mass / (r * sqr(tmp)));
}

/* gets negative of the acceleration vector of this disk component */
static inline mwvector miyamotoNagaiDiskAccel(const Disk* disk, const mwvector pos)
{
    mwvector acc;
    const real a   = disk->scaleLength;
    const real b   = disk->scaleHeight;
    const real zp  = mw_sqrt( sqr(Z(pos)) + sqr(b) );
    const real azp = a + zp;

    const real rp  = sqr(X(pos)) + sqr(Y(pos)) + sqr(azp);
    const real rth = mw_sqrt(cube(rp));  /* rp ^ (3/2) */

    X(acc) = -disk->mass * X(pos) / rth;
    Y(acc) = -disk->mass * Y(pos) / rth;
    Z(acc) = -disk->mass * Z(pos) * azp / (zp * rth);

    return acc;
}

static inline mwvector exponentialDiskAccel(const Disk* disk, const mwvector pos)
{
    const real b = disk->scaleLength;
    const real r = mw_absv(pos);

    const real expPiece = mw_exp(-r / b) * (r + b) / b;
    const real factor   = disk->mass * (expPiece - 1.0) / cube(r);

    return mw_mulvs(pos, factor);
}

static inline mwvector logHaloAccel(const Halo* halo, const mwvector pos)
{
    mwvector acc;

    const real tvsqr = -2.0 * sqr(halo->vhalo);
    const real qsqr  = sqr(halo->flattenZ);
    const real d     = halo->scaleLength;
    const real zsqr  = sqr(Z(pos));

    const real arst  = sqr(d) + sqr(X(pos)) + sqr(Y(pos));
    const real denom = (zsqr / qsqr) +  arst;

    X(acc) = tvsqr * X(pos) / denom;
    Y(acc) = tvsqr * Y(pos) / denom;
    Z(acc) = tvsqr * Z(pos) / ((qsqr * arst) + zsqr);

    return acc;
}

static inline mwvector nfwHaloAccel(const Halo* halo, const mwvector pos)
{
    const real r  = mw_absv(pos);
    const real a  = halo->scaleLength;
    const real ar = a + r;
    const real c  = a * sqr(halo->vhalo) * ((-ar * mw_log1p(r / a)) + r) / (0.2162165954 * cube(r) * ar);

    return mw_mulvs(pos, c);
}

/* CHECKME: Seems to have precision related issues for a small number of cases for very small qy */
static inline mwvector triaxialHaloAccel(const Halo* h, const mwvector pos)
{
    mwvector acc;

    /* TODO: More things here can be cached */
    const real qzs      = sqr(h->flattenZ);
    const real rhalosqr = sqr(h->scaleLength);
    const real mvsqr    = -sqr(h->vhalo);

    const real xsqr = sqr(X(pos));
    const real ysqr = sqr(Y(pos));
    const real zsqr = sqr(Z(pos));

    const real arst  = rhalosqr + (h->c1 * xsqr) + (h->c3 * X(pos) * Y(pos)) + (h->c2 * ysqr);
    const real arst2 = (zsqr / qzs) + arst;

    X(acc) = mvsqr * (((2.0 * h->c1) * X(pos)) + (h->c3 * Y(pos)) ) / arst2;

    Y(acc) = mvsqr * (((2.0 * h->c2) * Y(pos)) + (h->c3 * X(pos)) ) / arst2;

    Z(acc) = (2.0 * mvsqr * Z(pos)) / ((qzs * arst) + zsqr);

    return acc;
}

mwvector acceleration(const Potential* pot, const mwvector pos)
{
    mwvector acc, acctmp;

    /* GCC and clang both turn these into jump tables */
    switch (pot->disk.type)
    {
        case ExponentialDisk:
            acc = exponentialDiskAccel(&pot->disk, pos);
            break;
        case MiyamotoNagaiDisk:
            acc = miyamotoNagaiDiskAccel(&pot->disk, pos);
            break;
        case InvalidDisk:
        default:
            fail("Invalid disk type in acceleration()\n");
    }

    switch (pot->halo.type)
    {
        case LogarithmicHalo:
            acctmp = logHaloAccel(&pot->halo, pos);
            break;
        case NFWHalo:
            acctmp = nfwHaloAccel(&pot->halo, pos);
            break;
        case TriaxialHalo:
            acctmp = triaxialHaloAccel(&pot->halo, pos);
            break;
        case InvalidHalo:
        default:
            fail("Invalid halo type in acceleration()\n");
    }

    mw_incaddv(acc, acctmp);
    acctmp = sphericalAccel(&pot->sphere[0], pos);
    mw_incaddv(acc, acctmp);

    return acc;
}

