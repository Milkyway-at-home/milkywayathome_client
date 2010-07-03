/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "nbody_priv.h"

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
    const real zp  = rsqrt( sqr(pos[2]) + sqr(b) );
    const real azp = a + zp;
    const real rth = rpow( sqr(pos[0]) + sqr(pos[1]) + sqr(azp), 1.5);

    acc[0] = -disk->mass * pos[0] / rth;
    acc[1] = -disk->mass * pos[1] / rth;
    acc[2] = -disk->mass * pos[2] * azp / (zp * rth);
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
    const real zsqr  = sqr(pos[2]);

    const real arst  = sqr(d) + sqr(pos[0]) + sqr(pos[1]);
    const real denom = arst + zsqr / qsqr;

    acc[0] = tvsqr * pos[0] / denom;
    acc[1] = tvsqr * pos[1] / denom;
    acc[2] = tvsqr * pos[2] / (qsqr * arst + zsqr);
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

    const real xsqr = sqr(pos[0]);
    const real ysqr = sqr(pos[1]);
    const real zsqr = sqr(pos[2]);

    const real arst = rhalosqr + (h->c1 * xsqr) + (h->c3 * pos[0] * pos[1]) + (h->c2 * ysqr);

    const real arst2 = arst + zsqr / qzs;

    acc[0] = vsqr * ( 2 * h->c1 * pos[0] + h->c3 * pos[1] ) / arst2;

    acc[1] = vsqr * ( 2 * h->c2 * pos[1] + h->c3 * pos[0] ) / arst2;

    acc[2] = 2 * vsqr * pos[2] / (qzs * arst + zsqr);

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
    x[0] = ic->position[0];
    x[1] = ic->position[1];
    x[2] = ic->position[2];

    v[0] = -ic->velocity[0];
    v[1] = -ic->velocity[1];
    v[2] = -ic->velocity[2];

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
    fc->position[0] = x[0];
    fc->position[1] = x[1];
    fc->position[2] = x[2];

    fc->velocity[0] = -v[0];
    fc->velocity[1] = -v[1];
    fc->velocity[2] = -v[2];

    fc->useGalC = ic->useGalC;  /* Not actually necessary */
    fc->useRadians = ic->useRadians;

}


