/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "defs.h"
#include "code.h"

/* CHECKME: order of operations and effect on precision, and where can
 * we share divisions and such */

/* gets negative of the acceleration vector of this disk component */
static void miyamotoNagaiAccel(vector acc, const Disk* disk, const vector pos)
{
    const real a    = disk->scale_length;
    const real b    = disk->scale_height;
    const real zp   = sqrt( sqr(pos[2]) + sqr(b) );
    const real azp  = a + zp;
    const real rth  = pow( sqr(pos[0]) + sqr(pos[1]) + sqr(azp), 1.5);

    const real arst = disk->mass / rth;

    acc[0] = disk->mass * pos[0] / rth;
    acc[1] = disk->mass * pos[1] / rth;
    acc[2] = disk->mass * pos[2] * azp / (zp * rth);
}

static void logHaloAccel(vector acc, const Halo* halo, const vector pos)
{
    const real tvsqr = 2.0 * sqr(halo->vhalo);
    const real qsqr  = sqr(halo->flattenZ);
    const real d     = halo->scale_length;
    const real zsqr  = sqr(pos[2]);

    const real arst  = sqr(d) + sqr(pos[0]) + sqr(pos[1]);
    const real denom = arst + zsqr / qsqr;

    acc[0] = tvsqr * pos[0] / denom;
    acc[1] = tvsqr * pos[1] / denom;
    acc[2] = tvsqr * pos[2] / (qsqr * arst + zsqr);
}

static void sphericalAccel(vector acc, const Spherical* sph, const vector pos)
{
    const real r     = sqrt( sqr(pos[0]) + sqr(pos[1]) + sqr(pos[2]) );
    const real denom = r * sqr(sph->scale + r);

    acc[0] = sph->mass * pos[0] / denom;
    acc[1] = sph->mass * pos[1] / denom;
    acc[2] = sph->mass * pos[2] / denom;
}

inline static void acceleration(const NBodyCtx* ctx, const real* pos, real* acc)
{
    vector acc1 = { 0.0, 0.0, 0.0 };
    vector acc2 = { 0.0, 0.0, 0.0 };
    vector acc3 = { 0.0, 0.0, 0.0 };

    sphericalAccel(acc1, ctx->pot.sphere, pos);
    miyamotoNagaiAccel(acc2, &ctx->pot.disk, pos);
    logHaloAccel(acc3, &ctx->pot.halo, pos);

    acc[0] = - (acc1[0] + acc2[0] + acc3[0]);
    acc[1] = - (acc1[1] + acc2[1] + acc3[1]);
    acc[2] = - (acc1[2] + acc2[2] + acc3[2]);
}


void integrate(const NBodyCtx* ctx, InitialConditions* ic)
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
    acceleration(ctx, x, acc);

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
        acceleration(ctx, x, acc);
    }

    /* Report the final values (don't forget to reverse the velocities) */
    ic->position[0] = x[0];
    ic->position[1] = x[1];
    ic->position[2] = x[2];

    ic->velocity[0] = -v[0];
    ic->velocity[1] = -v[1];
    ic->velocity[2] = -v[2];

}


