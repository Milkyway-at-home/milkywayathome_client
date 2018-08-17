/*
 * Copyright (c) 2010, 2011 Matthew Arsenault
 * Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute.
 * Copyright (c) 2016-2018 Siddhartha Shelton
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

#include "nbody_priv.h"
#include "nbody_potential.h"
#include "milkyway_util.h"
#include "nbody_caustic.h"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

static inline mwvector hernquistSphericalAccel(const Spherical* sph, mwvector pos, real r)
{
    const real tmp = sph->scale + r;

    return mw_mulvs(pos, -sph->mass / (r * sqr(tmp)));
}

static inline mwvector plummerSphericalAccel(const Spherical* sph, mwvector pos, real r)
{
    const real tmp = mw_sqrt(sqr(sph->scale) + sqr(r));

    return mw_mulvs(pos, -sph->mass / cube(tmp));
}

/* gets negative of the acceleration vector of this disk component */
static inline mwvector miyamotoNagaiDiskAccel(const Disk* disk, mwvector pos, real r)
{
    mwvector acc;
    const real a   = disk->scaleLength;
    const real b   = disk->scaleHeight;
    const real zp  = mw_sqrt(sqr(Z(pos)) + sqr(b));
    const real azp = a + zp;

    const real rp  = sqr(X(pos)) + sqr(Y(pos)) + sqr(azp);
    const real rth = mw_sqrt(cube(rp));  /* rp ^ (3/2) */

    X(acc) = -disk->mass * X(pos) / rth;
    Y(acc) = -disk->mass * Y(pos) / rth;
    Z(acc) = -disk->mass * Z(pos) * azp / (zp * rth);

    return acc;
}

static inline mwvector exponentialDiskAccel(const Disk* disk, mwvector pos, real r)
{
    const real b = disk->scaleLength;

    const real expPiece = mw_exp(-r / b) * (r + b) / b;
    const real factor   = disk->mass * (expPiece - 1.0) / cube(r);

    return mw_mulvs(pos, factor);
}

static inline mwvector logHaloAccel(const Halo* halo, mwvector pos, real r)
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

static inline mwvector nfwHaloAccel(const Halo* halo, mwvector pos, real r)
{
    const real a  = halo->scaleLength;
    const real ar = a + r;
//     const real c  = a * sqr(halo->vhalo) * (r - ar * mw_log((r + a) / a)) / (0.2162165954 * cube(r) * ar);
    /* this is done to agree with NEMO. IDK WHY. IDK where 0.2162165954 comes from */
    const real c  = a * sqr(a) * 237.209949228 * (r - ar * mw_log((r + a) / a)) / ( cube(r) * ar);

    return mw_mulvs(pos, c);
}

/* CHECKME: Seems to have precision related issues for a small number of cases for very small qy */
static inline mwvector triaxialHaloAccel(const Halo* h, mwvector pos, real r)
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

static inline mwvector ASHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real gam = h->gamma;
    const real lam = h->lambda;
    const real M = h->mass;
    const real a = h->scaleLength;
    const real scaleR = r/a;
    const real scaleL = lam/a;
    real c;

    if (r<lam)
    {    c = -(M/sqr(a))*mw_pow(scaleR,gam-2)/(1+mw_pow(scaleR,gam-1));
    }
    else
    {    c = -(M/sqr(r))*mw_pow(scaleL,gam)/(1+mw_pow(scaleL,gam-1));
    }

    return mw_mulvs(pos, c/r);
}

static inline mwvector WEHaloAccel(const Halo* h, mwvector pos, real r)
{
    const real a = h->scaleLength;
    const real M = h->mass;
    const real sum2 = sqr(r) + sqr(a);

    const real c = (-M/r)*(a + mw_sqrt(sum2))/(sum2 + a*mw_sqrt(sum2));

    return mw_mulvs(pos, c/r);

}

mwvector nbExtAcceleration(const Potential* pot, mwvector pos)
{
    mwvector acc, acctmp;
    const real r = mw_absv(pos);
    /*Calculate the Disk Accelerations*/
    switch (pot->disk.type)
    {
        case ExponentialDisk:
            acc = exponentialDiskAccel(&pot->disk, pos, r);
            break;
        case MiyamotoNagaiDisk:
            acc = miyamotoNagaiDiskAccel(&pot->disk, pos, r);
            break;
        case InvalidDisk:
        default:
            mw_fail("Invalid disk type in external acceleration\n");
    }
    /*Calculate the Halo Accelerations*/
    switch (pot->halo.type)
    {
        case LogarithmicHalo:
            acctmp = logHaloAccel(&pot->halo, pos, r);
            break;
        case NFWHalo:
            acctmp = nfwHaloAccel(&pot->halo, pos, r);
            break;
        case TriaxialHalo:
            acctmp = triaxialHaloAccel(&pot->halo, pos, r);
            break;
        case CausticHalo:
            acctmp = causticHaloAccel(&pot->halo, pos, r);
            break;
        case AllenSantillanHalo:
            acctmp = ASHaloAccel(&pot->halo, pos, r);
            break;
        case WilkinsonEvansHalo:
            acctmp = WEHaloAccel(&pot->halo, pos, r);
            break;
        case InvalidHalo:
        default:
            mw_fail("Invalid halo type in external acceleration\n");
    }

    mw_incaddv(acc, acctmp);
    /*Calculate the Bulge Accelerations*/
    switch (pot->sphere[0].type)
    {
        case HernquistSpherical:
            acctmp = hernquistSphericalAccel(&pot->sphere[0], pos, r);
            break;
        case PlummerSpherical:
            acctmp = plummerSphericalAccel(&pot->sphere[0], pos, r);
            break;
        case InvalidSpherical:
        default:
            mw_fail("Invalid bulge type in external acceleration\n");
    }

    mw_incaddv(acc, acctmp);

    return acc;
}













