/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "nbody_priv.h"
#include "orbitintegrator.h"

/* CHECKME: order of operations and effect on precision, and where can
 * we share divisions and such */


/* Pure functions are the best ones */

void sphericalAccel(vectorptr restrict acc, const Spherical* sph, const vectorptr restrict pos)
{
    real r;
    ABSV(r, pos);
    MULVS(acc, pos, -sph->mass / (r * sqr(sph->scale + r)));
}

/* gets negative of the acceleration vector of this disk component */
void miyamotoNagaiDiskAccel(vectorptr restrict acc, const Disk* disk, const vectorptr restrict pos)
{
    const real a   = disk->scale_length;
    const real b   = disk->scale_height;
    const real zp  = rsqrt( sqr(Z(pos)) + sqr(b) );
    const real azp = a + zp;
    const real rth = rpow( sqr(X(pos)) + sqr(Y(pos)) + sqr(azp), 1.5);

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
    const real factor   = disk->mass * (expPiece - 1) / cube(r);
    MULVS(acc, pos, factor);
}

void logHaloAccel(vectorptr restrict acc, const Halo* halo, const vectorptr restrict pos)
{
    const real tvsqr = -2.0 * sqr(halo->vhalo);
    const real qsqr  = sqr(halo->flattenZ);
    const real d     = halo->scale_length;
    const real zsqr  = sqr(Z(pos));

    const real arst  = sqr(d) + sqr(X(pos)) + sqr(Y(pos));
    const real denom = arst + zsqr / qsqr;

    X(acc) = tvsqr * X(pos) / denom;
    Y(acc) = tvsqr * Y(pos) / denom;
    Z(acc) = tvsqr * Z(pos) / (qsqr * arst + zsqr);
}

void nfwHaloAccel(vectorptr restrict acc, const Halo* halo, const vectorptr restrict pos)
{
    real r;
    ABSV(r, pos);
    const real a  = halo->scale_length;
    const real ar = a + r;
    const real c  = a * sqr(halo->vhalo) * (r - ar * rlog1p(r / a)) / (0.216 * cube(r) * ar);

    MULVS(acc, pos, c);
}

/* CHECKME: Seems to have precision related issues for a small number of cases for very small qy */
void triaxialHaloAccel(vectorptr restrict acc, const Halo* h, const vectorptr restrict pos)
{
    /* TODO: More things here can be cached */
    const real qzs = sqr(h->flattenZ);
    const real rhalosqr = sqr(h->scale_length);
    const real vsqr = -sqr(h->vhalo);

    const real xsqr = sqr(X(pos));
    const real ysqr = sqr(Y(pos));
    const real zsqr = sqr(Z(pos));

    const real arst = rhalosqr + (h->c1 * xsqr) + (h->c3 * X(pos) * Y(pos)) + (h->c2 * ysqr);

    const real arst2 = arst + zsqr / qzs;

    X(acc) = vsqr * ( 2 * h->c1 * X(pos) + h->c3 * Y(pos) ) / arst2;

    Y(acc) = vsqr * ( 2 * h->c2 * Y(pos) + h->c3 * X(pos) ) / arst2;

    Z(acc) = 2 * vsqr * Z(pos) / (qzs * arst + zsqr);

}

inline void acceleration(vectorptr restrict acc, const NBodyCtx* ctx, const vectorptr restrict pos)
{
    /* lookup table for functions for calculating accelerations */
    static const HaloAccel haloFuncs[] = { [LogarithmicHalo] = logHaloAccel,
                                           [NFWHalo]         = nfwHaloAccel,
                                           [TriaxialHalo]    = triaxialHaloAccel };

    static const DiskAccel diskFuncs[] = { [ExponentialDisk]   = exponentialDiskAccel,
                                           [MiyamotoNagaiDisk] = miyamotoNagaiDiskAccel };

    static const SphericalAccel sphFuncs[] = { [SphericalPotential] = sphericalAccel };

    vector acctmp1 = ZERO_VECTOR;
    vector acctmp2 = ZERO_VECTOR;

    /* Use the type of potential to index into the table, and use the
     * appropriate function */

    diskFuncs[ctx->pot.disk.type](acctmp1, &ctx->pot.disk, pos);
    haloFuncs[ctx->pot.halo.type](acctmp2, &ctx->pot.halo, pos);

    ADDV(acc, acctmp1, acctmp2);

    sphFuncs[ctx->pot.sphere[0].type](acctmp1, &ctx->pot.sphere[0], pos);

    /* add the resulting vectors */
    INCADDV(acc, acctmp1);
}


void reverseOrbit(InitialConditions* fc, const NBodyCtx* ctx, InitialConditions* ic)
{
    vector acc, v, x;
    real t;
    int i;

    const real tstop = ctx->model.time_orbit;
    const real dt    = ctx->model.orbit_timestep;

    // Set the initial conditions
    SETV(x, ic->position);
    SETV(v, ic->velocity);
    INCNEGV(v);

    // Get the initial acceleration
    acceleration(acc, ctx, x);

    // Loop through time
    for (t = 0; t <= tstop; t += dt)
    {
        // Update the velocities and positions
        for (i = 0; i < 3; ++i)
        {
            v[i] += acc[i] * dt;
            x[i] += v[i] * dt;
        }
        // Compute the new acceleration
        acceleration(acc, ctx, x);
    }

    /* Report the final values (don't forget to reverse the velocities) */
    SETV(fc->position, x);
    SETV(fc->velocity, v);
    INCNEGV(fc->velocity);

    fc->useGalC = ic->useGalC;  /* Not actually necessary */
    fc->useRadians = ic->useRadians;

}


